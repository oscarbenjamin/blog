from collections import OrderedDict
import weakref

#
# A WeakValueDictionary maps keys to values but will discard items as soon as
# there are no other references to the value. We are going to use this to map
# from a tuples of child graphs to graph objects. Then any time we try to
# create an expression that already exists the previously created expression
# will be returned. Note that unlike a cache this will not hold anything in
# memory that is not still being used somewhere.
#
all_expressions = weakref.WeakValueDictionary()

#
# The all_atoms dictionary is similar to all expressions but it maps from the
# internal hashable data of the atomic expression to the atomic expression
# value.
#
all_atoms = weakref.WeakValueDictionary()


class Expression:
    """Class of all expressions.

    This has subclasses for atomic and compound expressions respectively.
    """

    # For now we allow any expression to be a head.
    def __call__(*args):
        return TreeExpression(*args)

    def __neg__(self):
        return Mul(Integer(-1), self)

    def __add__(self, other):
        return Add(self, other)

    def __sub__(self, other):
        return Add(self, Mul(Integer(-1), other))

    def __mul__(self, other):
        return Mul(self, other)

    def __pow__(self, other):
        return Pow(self, other)

    #
    # Here we do not define methods like __eq__ or __hash__. The default
    # methods provided by object are fine because those assume that equality is
    # based on identity. The use of the all_expressions and all_atoms
    # dictionary above ensures that equality is based on identity because each
    # expression is globally unique.
    #


class AtomType:
    """Identifier to distinguish different kinds of atoms.

    An AtomType is not an Expression but is used to construct atomic
    expressions.

    >>> Integer = AtomType('Integer')
    >>> Integer
    Integer
    >>> Integer(1)
    Integer(1)
    """
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name

    def __call__(self, value):
        return AtomExpression(self, value)


class AtomExpression(Expression):
    """The class for all atomic expressions

    An atomic expression has an internal value but no args.
    """
    def __new__(cls, atom_type, value):
        #
        # Include type(value) in the key to avoid 1 == 1.0 in all_atoms.
        #
        key = (atom_type, value, type(value))

        previous = all_atoms.get(key, None)
        if previous is not None:
            return previous

        expression = object.__new__(AtomExpression)

        expression.children = ()
        expression.head = atom_type
        expression.args = ()
        expression.value = value

        all_atoms[key] = expression

        return expression

    def __repr__(self):
        return str(self.value)


class TreeExpression(Expression):
    """The class for all compound expressions

    These expose their args.
    """
    def __new__(cls, *children):
        # Return a previously created expression if possible
        previous = all_expressions.get(children, None)
        if previous is not None:
            return previous

        expression = object.__new__(cls)

        expression.children = children
        expression.head = children[0]
        expression.args = children[1:]
        expression.value = None

        all_expressions[children] = expression
        return expression

    def __repr__(self):
        return f'{self.head}({", ".join(map(repr, self.args))})'

    def __repr__(self):
        return to_str(self)


def topological_sort(expression):
    """List of subexpressions sorted topologically.

    >>> f = Function('f')
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> topological_sort(f(f(x, y), f(f(x))))
    [f, x, y, f(x, y), f(x), f(f(x)), f(f(x, y), f(f(x)))]
    """
    seen = set()
    expressions = []
    stack = [(expression, list(expression.children)[::-1])]

    while stack:
        top, children = stack[-1]
        while children:
            child = children.pop()
            if child not in seen:
                seen.add(child)
                stack.append((child, list(child.children)[::-1]))
                break
        else:
            stack.pop()
            expressions.append(top)

    return expressions


def tree_to_instructions(expression):

    subexpressions = topological_sort(expression)

    atoms = []
    evaluation = []
    compound = []

    indices = {}

    for expr in subexpressions:
        if not expr.children:
            atoms.append(expr)
        else:
            compound.append(expr)

    indices = dict(zip(atoms, range(len(atoms))))

    for index, expr in enumerate(compound, len(atoms)):
        child_indices = [indices[e] for e in expr.children]
        evaluation.append(child_indices)
        indices[expr] = index

    return atoms, evaluation


def instructions_to_tree(atoms, evaluation):

    stack = list(atoms)

    for indices in evaluation:
        args = [stack[i] for i in indices]
        stack.append(TreeExpression(*args))

    expression = stack[-1]
    return expression


def evaluation_form(expression, operators):

    subexpressions = topological_sort(expression)

    atoms = []
    operations = []
    compound = []

    for expr in subexpressions:
        if expr in operators:
            continue
        elif not expr.children:
            atoms.append(expr)
        else:
            compound.append(expr)

    indices = dict(zip(atoms, range(len(atoms))))

    for n, expr in enumerate(atoms):
        operator = operators.get(expr.head, None)
        if operator is not None:
            atoms[n] = operator(expr.value)

    for index, expr in enumerate(compound, len(atoms)):
        head = expr.head
        function = operators.get(head, head)
        arg_indices = [indices[e] for e in expr.args]
        operations.append((function, arg_indices))
        indices[expr] = index

    return atoms, operations


def evaluate(expression, operators, values):

    stack, operations = evaluation_form(expression, operators)

    stack = [values.get(e, e) for e in stack]

    for func, indices in operations:
        args = [stack[i] for i in indices]
        stack.append(func(*args))

    return stack[-1]


Integer = AtomType('Integer')
Symbol = AtomType('Symbol')
Function = AtomType('Function')

Add = Function('Add')
Mul = Function('Mul')
Pow = Function('Pow')

x = Symbol('x')
y = Symbol('y')
sin = Function('sin')
cos = Function('cos')

one = Integer(1)
zero = Integer(0)

class UserExpression:

    def __init__(self, rep):
        self._rep = rep

    def __repr__(self):
        return to_str(self._rep)

    def __neg__(self):
        return Mul(Integer(-1), self._rep)

    def __add__(self, other):
        return Add(self._rep, other._rep)

    def __sub__(self, other):
        return Add(self._rep, Mul(Integer(-1), other._rep))

    def __mul__(self, other):
        return Mul(self._rep, other._rep)

    def __pow__(self, other):
        return Pow(self._rep, other._rep)

    def eval_f64(self, values=None):
        if values is None:
            values = {}
        return eval_f64(self._rep, values)

    def diff(self, symbol, ntimes=1):
        derivative = self._rep
        for n in range(ntimes):
            derivative = diff(derivative, symbol)
        return UserExpression(derivative)

import math

operators_f64 = {
    Integer: float,
    Add: lambda *a: math.fsum(a),
    Mul: lambda *a: math.prod(a),
    Pow: pow,
    sin: math.sin,
    cos: math.cos,
}

def eval_f64(expression, values=None):
    if values is None:
        values = {}
    return evaluate(expression, operators_f64, values)

import sympy

operators_sympy = {
    Integer: sympy.Integer,
    Symbol: sympy.Symbol,
    Add: sympy.Add,
    Mul: sympy.Mul,
    Pow: sympy.Pow,
    sin: sympy.sin,
    cos: sympy.cos,
}

def to_sympy(expression):
    return evaluate(expression, operators_sympy, {})

operators_to_str = {
    Integer: str,
    Symbol: str,
    Add: lambda *a: f'({" + ".join(a)})',
    Mul: lambda *a: f'({"*".join(a)})',
    Pow: lambda b, e: f'{b}^{e}',
    sin: lambda a: f'sin({a})',
    cos: lambda a: f'cos({a})',
}

def to_str(expression):
    return evaluate(expression, operators_to_str, {})

derivatives = {
    (sin, 0): cos,
    (cos, 0): lambda e: -sin(e),
    (Pow, 0): lambda b, e: Pow(b, e - 1),
}

derivatives_applied = {Add}


def diff(expression, sym):
    """Derivative of expression wrt sym.

    Uses forward accumulation algorithm.
    """
    stack, operations = evaluation_form(expression, {})

    diff_stack = [one if expr == sym else zero for expr in stack]

    for func, indices in operations:

        args = [stack[i] for i in indices]
        diff_args = [diff_stack[i] for i in indices]
        expr = func(*args)

        if set(diff_args) == {zero}:
            diff_terms = []
        elif func == Add:
            diff_terms = [da for da in diff_args if da != zero]
        elif func == Mul:
            diff_terms = []
            for n, diff_arg in enumerate(diff_args):
                if diff_arg != zero:
                    term = Mul(*args[:n], diff_arg, *args[n+1:])
                    diff_terms.append(term)
        elif func == Pow:
            assert diff_args[1] == zero
            base, exp = args
            diff_terms = [Mul(exp, Pow(base, exp-Integer(1)))]
        else:
            diff_terms = []
            for n, diff_arg in enumerate(diff_args):
                if diff_arg != zero:
                    pdiff = derivatives[(func, n)]
                    diff_term = pdiff(*args)
                    if diff_arg != one:
                        diff_term = Mul(diff_term, diff_arg)
                    diff_terms.append(diff_term)

        if not diff_terms:
            derivative = zero
        elif len(diff_terms) == 1:
            derivative = diff_terms[0]
        else:
            derivative = Add(*diff_terms)

        stack.append(expr)
        diff_stack.append(derivative)

    return diff_stack[-1]


def make_expression(n):
    e = x**Integer(2) - Integer(1)
    for _ in range(n):
        e = e**Integer(2) + e
    return e

def make_expression2(x, n):
    e = x**2 - 1
    for _ in range(n):
        e = e**2 + e
    return e

e = sin(sin(sin(x)))
for _ in range(10):
    e = diff(e, x)
