from __future__ import annotations

from collections import OrderedDict
import weakref

#
# A WeakValueDictionary maps keys to values but will discard items as soon as
# there are no other references to the value. We are going to use this to map
# from tuples of child graphs to graph objects. Then any time we try to create
# an expression that already exists the previously created expression will be
# returned. Note that unlike a cache this will not hold anything in memory that
# is not still being used somewhere.
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


def rebuild(atoms, evaluation):
    """Inverse of evaluation_form(expression, {})"""
    stack = list(atoms)

    for indices in evaluation:
        args = [stack[i] for i in indices]
        stack.append(TreeExpression(*args))

    expression = stack[-1]
    return expression


def evaluate(expression, operators, values):
    """Evaluate the expression"""
    stack, operations = evaluation_form(expression, operators)

    stack = [values.get(e, e) for e in stack]

    for func, indices in operations:
        args = [stack[i] for i in indices]
        stack.append(func(*args))

    return stack[-1]


# ------------------------------------------------------- #
#   Define some basic expression types 		          #
# ------------------------------------------------------- #


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
negone = Integer(-1)
zero = Integer(0)


# ------------------------------------------------------- #
#   eval_f64: 64-bit real floating point evaluation.      #
# ------------------------------------------------------- #


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


# ------------------------------------------------------- #
#   eval_int: Integer evaluation.                         #
# ------------------------------------------------------- #


import math

operators_int = {
    Integer: int,
    Add: lambda *a: sum(a),
    Mul: lambda *a: math.prod(a),
    Pow: pow,
}

def eval_int(expression, values=None):
    if values is None:
        values = {}
    return evaluate(expression, operators_int, values)


# ------------------------------------------------------- #
#   to_sympy: Convert to SymPy expression.                #
# ------------------------------------------------------- #


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


operators_from_sympy = {
    sympy.Integer: lambda e: Integer(e.p),
    sympy.Symbol: lambda e: Symbol(e.name),
    sympy.Add: Add,
    sympy.Mul: Mul,
    sympy.Pow: Mul,
    sympy.Mul: Mul,
    sympy.sin: sin,
    sympy.cos: cos,
}


def from_sympy(expression):
    if expression.is_Atom:
        if isinstance(expression, sympy.Integer):
            return Integer(expression.p)
        elif isinstance(expression, sympy.Symbol):
            return Symbol(expression.name)
        else:
            assert False
    args = [from_sympy(arg) for arg in expression.args]
    func = operators_from_sympy[type(expression)]
    return func(*args)


# ------------------------------------------------------- #
#   to_str: Printing support				              #
# ------------------------------------------------------- #


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


# ------------------------------------------------------- #
#   count_ops: statistics about expression size           #
# ------------------------------------------------------- #


op1 = lambda e: 1
sum1 = lambda *a: 1 + sum(a)


operators_count_ops = {
    Integer: op1,
    Symbol: op1,
    Add: sum1,
    Mul: sum1,
    Pow: sum1,
    sin: sum1,
    cos: sum1,
}


def count_ops_tree(expression):
    return evaluate(expression, operators_count_ops, {})


def count_ops_graph(expression):
    return len(topological_sort(expression))


# ------------------------------------------------------- #
#   diff: Differentiation                                 #
# ------------------------------------------------------- #


derivatives = {
    (sin, 0): cos,
    (cos, 0): lambda e: -sin(e),
    (Pow, 0): lambda b, e: Pow(b, e - 1),
}


def _diff(expression, sym):
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
            diff_terms = [Mul(exp, Pow(base, Add(exp, negone)))]
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


def _diff_reverse(expr, sym):
    """Derivative of expression wrt sym.

    Uses reverse mode accumulation.
    """
    stack, operations = evaluation_form(expr, {})

    try:
        sym_index = stack.index(sym)
    except ValueError:
        return zero

    # Forward pass for evaluation of the function
    for func, indices in operations:
        args = [stack[i] for i in indices]
        stack.append(func(*args))

    # Reverse pass to collect the derivatives
    nterms = len(stack)
    diff_terms = [[] for _ in range(nterms)]
    diff_terms[-1] = [one]

    for func, indices in reversed(operations):
        args = [stack[i] for i in indices]
        terms = diff_terms.pop()
        if len(terms) == 1:
            ddf = terms[0]
        else:
            ddf = Add(*terms)
        if func == Add:
            for i in indices:
                diff_terms[i].append(ddf)
        elif func == Mul:
            for n, i in enumerate(indices):
                otherargs = args[:n] + args[n+1:] + [ddf]
                otherargs = [arg for arg in otherargs if arg != one]
                if len(otherargs) == 0:
                    ddfi = one
                elif len(otherargs) == 1:
                    ddfi = otherargs[0]
                else:
                    ddfi = Mul(*otherargs)
                diff_terms[i].append(ddfi)
        else:
            for n, i in enumerate(indices):
                pdiff_func = derivatives[(func, n)]
                pdiff = pdiff_func(*args)
                if ddf == one:
                    ddfi = pdiff
                else:
                    ddfi = Mul(ddf, pdiff)
                diff_terms[i].append(ddfi)

    terms_final = diff_terms[sym_index]

    if len(terms_final) == 1:
        dfdx = terms_final[0]
    else:
        dfdx = Add(*terms_final)

    return dfdx


def diff(expr, sym, ntimes=1, reverse=False):
    diff1 = _diff_reverse if reverse else _diff

    deriv = expr
    for _ in range(ntimes):
        deriv = diff1(deriv, sym)
    return deriv


# -------------------------------------------------
# binexpand : replace associative operators with binary operators
# -------------------------------------------------

from functools import reduce

def _binexpand(operator, args):
    if len(args) == 0:
        raise ValueError
    return reduce(operator, args)

operators_binexpand = {
    Add: lambda *a: _binexpand(Add, a),
    Mul: lambda *a: _binexpand(Mul, a),
}

def binexpand(expression):
    return evaluate(expression, operators_binexpand, {})

# -------------------------------------------------
# lambdification with LLVM
# ------------------------------------------------

import struct

def double_to_hex(f):
    return hex(struct.unpack('<Q', struct.pack('<d', f))[0])

operators_llvm = {
    Add: Add,
    Mul: Mul,
    Pow: Pow,
    cos: cos,
    sin: sin,
    Integer: float,
}

_llvm_header = """
; ModuleID = "mod1"
target triple = "unknown-unknown-unknown"
target datalayout = ""

declare double    @llvm.pow.f64(double %Val1, double %Val2)
declare double    @llvm.sin.f64(double %Val)
declare double    @llvm.cos.f64(double %Val)

"""

def to_llvm_f64(symargs, expression):

    expression = binexpand(expression)

    atoms, operations = evaluation_form(expression, operators_llvm)

    argnames = {s: f'%"s"' for s in symargs}
    constants = {}

    identifiers = []
    for a in atoms:
        if a in symargs:
            identifiers.append(argnames[a])
        elif isinstance(a, float):
            identifiers.append(double_to_hex(a))
        else:
            raise ValueError

    args = ', '.join(f'double {argnames[arg]}' for arg in symargs)
    signature = f'define double @"jit_func1"({args})'

    instructions = []
    for func, indices in operations:

        n = len(instructions)
        identifier = f'%".{n}"'
        identifiers.append(identifier)
        argids = [identifiers[i] for i in indices]

        if func == Add:
            instructions.append(f'{identifier} = fadd double ' + ', '.join(argids))
        elif func == Mul:
            instructions.append(f'{identifier} = fmul double ' + ', '.join(argids))
        elif func == Pow:
            instructions.append(f'{identifier} = call double @llvm.pow.f64(double {argids[0]}, double {argids[1]})')
        elif func == sin:
            instructions.append(f'{identifier} = call double @llvm.sin.f64(double {argids[0]})')
        elif func == cos:
            instructions.append(f'{identifier} = call double @llvm.cos.f64(double {argids[0]})')
        else:
            raise ValueError(func)

    instructions.append(f'ret double {identifiers[-1]}')

    function_lines = [signature, '{', *instructions, '}']
    module_code = _llvm_header + '\n'.join(function_lines)
    return module_code


_exe_eng = []


def lambdify(args, expression):
    module_code = to_llvm_f64(args, expression)

    import ctypes
    import llvmlite.binding as llvm

    llvm.initialize()
    llvm.initialize_native_target()
    llvm.initialize_native_asmprinter()

    llmod = llvm.parse_assembly(module_code)

    pmb = llvm.create_pass_manager_builder()
    pmb.opt_level = 2
    pass_manager = llvm.create_module_pass_manager()
    pmb.populate(pass_manager)

    pass_manager.run(llmod)

    target_machine = llvm.Target.from_default_triple().create_target_machine()
    exe_eng = llvm.create_mcjit_compiler(llmod, target_machine)
    exe_eng.finalize_object()
    _exe_eng.append(exe_eng)

    fptr = exe_eng.get_function_address('jit_func1')

    rettype = ctypes.c_double
    argtypes = [ctypes.c_double] * len(args)

    cfunc = ctypes.CFUNCTYPE(rettype, *argtypes)(fptr)
    return cfunc



def make_expression(n):
    e = x**Integer(2) - Integer(1)
    for _ in range(n):
        e = e**Integer(2) + e
    return e

def make_expression2(n, x, intfunc=int):
    e = x**intfunc(2) - intfunc(1)
    for _ in range(n):
        e = e**intfunc(2) + e
    return e

def make_expression3(n, x, intfunc=int):
    #return Poly(range(0, n+1)[::-1], x)
    terms = [intfunc(m)*x**intfunc(m) for m in range(1, n+1)]
    #return sum(terms)
    return Add(*terms)

#e = sin(sin(sin(x)))
#diff(e, x, 10)

#diff(sin(x), x, 1000)

#print(1)
#e = make_expression(100)
#print(2)
#print(eval_f64(e, {x:0}))


#for _ in range(10):
#    e = diff(e, x)
