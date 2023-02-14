.. _data-structures-1:

Data structures for symbolic expressions: part 1
================================================

This will be the first of many blog posts discussing data structures and basic
algorithms for working with symbolic expressions. The reason for this is that I
want to change the way that SymPy's data structures work and I need to make
sure that lots of other people understand the design decisions that I am
considering. Firstly I need to introduce how things currently work in SymPy to
understand the context.

The design of Basic
-------------------

Let's use SymPy to consider what a data structure for a symbolic expression
might look like. In SymPy we can create symbolic expressions by creating
symbols and then calling symbolic functions on those symbols

.. doctest::

    >>> from sympy import symbols, cos
    >>> x = symbols('x')
    >>> expr = (1 + cos(x))**2
    >>> expr
    (cos(x) + 1)**2

Here we have created a symbolic expression ``expr`` which represents the
mathematical expression :math:`(1 + \cos{x})^2`. The expression as a data
structure is represented as a tree where each node in the tree is particular
function and the children are the *arguments* of that function. We can see this
functional representation explicitly using SymPy's ``srepr`` function

.. doctest::

    >>> from sympy import srepr
    >>> print(srepr(expr))
    Pow(Add(cos(Symbol('x')), Integer(1)), Integer(2))

Here the functions are ``Pow``, ``Add`` and ``cos`` and these represent
mathematical operations and functions so e.g. ``Pow(x, y)`` represents
:math:`x^y`, ``Add(x, y)`` represents :math:`x + y` and ``cos`` represents the
:math:`\cos` function. In an expression ``Add(x, y)`` we say that ``Add`` is
the *head* of the expression and that ``x`` and ``y`` are the arguments. The
other parts of the expression are ``Symbol('x')`` and ``Integer(1)`` and
``Integer(2)``. These are the *atoms* of the expression and represent the leaf
nodes of the expression tree. Using SymPy's ``dotprint`` functon and the
graphviz program you can make a visual representation of an expression tree
that looks like this:

.. graphviz:: expr.dot

In SymPy's terminology the expression head is known as the ``func`` and the
arguments as the ``args``.

.. doctest::

    >>> expr.func
    <class 'sympy.core.power.Pow'>
    >>> expr.args
    (cos(x) + 1, 2)

Here the ``func`` is ``Pow`` and the ``args`` are ``cos(x) + 1`` and ``2`` so
this expression represents ``Pow(cos(x) + 1, 2)`` which is :math:`(\cos{x} +
1)^2`.  Of course in turn the expression ``cos(x) + 1`` is represented as an
``Add`` with its ``args`` and so on until we reach the atoms which are the
expressions that do not have ``args``. It might look like ``Integer(2)`` has
the argument ``1`` but this is not considered to be part of the ``args``:

.. doctest::

    >>> from sympy import Integer
    >>> Integer(2).func
    <class 'sympy.core.numbers.Integer'>
    >>> Integer(2).args
    ()

Of course internally the ``Integer(2)`` object holds an ordinary Python integer
representing the number ``2`` to distinguish itself from ``Integer(3)`` but the
expression tree structure considers that ``Integer(2)`` is atomic and does not
recursively have more children so its ``args`` tuple is empty.

Many other symbolic computational systems use thes basic notion of a head (or
function) and arguments as the main way to represent symbolic expressions but
there are lots of subtle details and differences in how exactly that works.
One thing that you might notice in the output above is that the head of a SymPy
expression is a Python ``class`` object. In fact an expression of a given head
is an *instance* of that class:

.. doctest::

    >>> from sympy import Pow
    >>> expr.func == Pow
    True
    >>> type(expr) == Pow
    True

So if ``Pow`` is a class then calling ``Pow(x, y)`` creates an instance of
``Pow``. The instance can then store its arguments ``x`` and ``y`` and expose
them back through the ``args`` attribute. If we want to make new expression
heads we can make new classes. Many operations should behave differently for
different expression heads and with classes we can add *methods* which are
functions whose behavour can be different for each different class. This allows
u to encode the differences between all of the different mathematical functions
and operations when operating on symbolic expressions. Of course we want
expressions to have a relatively uniform interface for end users so we want all
of these different classes to present mostly the same attributes and methods
even if they would behave differently for the different classes. To make this
happen we use superclasses and then have lots of different subclasses that only
change a few of those methods. In SymPy the top-level superclass is called
``Basic`` and it has many levels of subclasses in a hierarchy. A few classes
from the top of this hierarchy are arranged like this:

.. graphviz:: classes.dot

Many of these break down further into more subclasses for example ``Boolean``
have many subclasses in another hierarchy:

.. graphviz:: bool.dot

The idea here is that the top-level ``Basic`` type defines methods that are
common to all symbolic expressions like ``subs`` for substituting values. The
next level subclasses like ``Expr``, ``Boolean`` and ``Set`` define other
methods that only make sense for the kinds of mathematical operations they
represent. For example:

- ``Expr`` represents things that can be added and multiplied so it defines
  methods like ``__add__`` and ``__mul__`` that will create ``Add`` and ``Mul``
  objects when you do e.g. ``expr1 * expr2``.
- ``Set`` represents mathematical sets and defines methods like ``intersect``
  that will create ``Intersection`` and ``Union`` expressions.
- ``Boolean`` represents boolean expressions like ``x > 2`` and so ``Boolean``
  defines methods like ``__and__`` and ``__or__`` which can be used to make
  ``And`` and ``Or`` expressions.

Then each of these high-level ``Basic`` subclasses will have more subclasses
that define more specific behaviour. For example ``Gt`` is for greater-than and
``Gt(x, y)`` represents the boolean statement :math:`x > y`. The `Gt` class has
methods that can work out if the statement is true and can then evaluate to
``true`` or ``false``:

.. doctest::

    >>> from sympy import Gt
    >>> Gt(x, 1)
    x > 1
    >>> Gt(1, 2)
    False
    >>> Gt(2, 1)
    True

The class ``Lt`` represents less-than and so its methods for evaluating to true
or false need to be different from those of ``Gt``. Likewise the ``cos`` and
``sin`` classes will have different methods for numerical evaluation and
differentiation and so on.

Using classes like this to represent expression heads has a number of
advantages because it means that we can leverage the general techniques of
object-oriented programming to build up all the different operations that we
would want in a computer algebra system. However I want to use this document to
explain why this is actually a problematic design decision. In fact in a new
design of how SymPy's symbolic expressions are represented we should *not* use
classes for expression heads.


Expression heads should not be classes
--------------------------------------

SymPy is not the only computer algebra system or symbolic manipulation system
to use classes for expression heads. For example ``SymEngine`` is a C++ library
based on essentially the same design as SymPy but using C++ classes instead.
The core symbolic library in SageMath is based on Pynac and is based on GiNaC
which is a C++ library. Each of these systems follows much of the design that I
referred to above where there are classes like ``Add``, ``Mul`` etc having
different methods to give different behaviour. Exactly how the methods work and
how the internal data structures work is different in each case but the basic
principles are the same: each different expression head is a different class
and the meaning of the different expression heads is encoded in the fact that
they have methods that do different things in different subclasses.

There are more examples of class-based symbolic systems that I could list but
broadly they are similar for the purposes of our discussion here. What I do
want to say though is that these different systems have different scopes and
are used in different ways and so the arguments that I will make for not using
classes do not necessarily apply to them as well. The system that I know best
is SymPy and in that context I can see that the use of classes for expression
heads works out badly.

At the same time there are other symbolic computing engines that are not based
on the object-oriented design. For example Mathematica allows any variable to
be either a symbol or a function and in Mathematica the internal representation
for something like ``cos(x)`` (spelled ``Cos[x]`` in Wolfram language) is
something more like a tuple::

    (Cos, x)

where both ``Cos`` and ``x`` are just atomic symbols. There is no basic
distinction between what kind of objects ``Cos`` and ``x`` are here::

    In[1]:= expr1 = Cos[x]
    Out[1]:= Cos[x]

    In[2]:= expr2 = x[Cos]
    Out[2]:= x[Cos]

Indexing into the expression gives the ``args`` but also the head is simply
the first of the ``args``::

    In[3]:= expr2[[0]]
    Out[3]:= x

    In[4]:= expr2[[1]]
    Out[4]:= Cos

Not distinguishing the head from the ``args`` makes it possible for the head
itself to be a compound symbolic expression. Here we differentiate a function
``f(x)`` and we can see that the head is now an expression::

    In[7]:= df = D[f[x], x]
    Out[7]:= f'[x]

    In[8]:= df[[0]]
    Out[8]:= f'

    In[9]:= df[[0]] // FullForm
    Out[9]:= Derivative[1][f]

This is something that just is not directly possible in SymPy's class-based
system. The reason is that in SymPy a head is a class and an expression is an
instance. The ``args`` of an expression also need to be instances. So now if
``f(x)`` makes sense as an expression ``f`` must be a class. But if ``f`` is a
class then it cannot be an expression and so it cannot be in the ``args`` of an
expression like ``Derivative``: only its instances can. Hence in SymPy we can
represent ``f'(x)`` only as ``Derivative(f(x), x)`` but there is no way to
represent ``f'`` as simply the derivative of the function ``f`` without that
function having been called with an argument. An operation that is natural in
Mathematica is to differentiate an undefined function and then substitute a
different value for its argument::

    In[10]:= D[f[x], x] /. x -> y^2
    Out[10]:= f'[y ^ 2]

    In[11]:= f'[y^2]
    Out[11]:= f'[y ^ 2]

    In[12]:= f'[y^2] // FullForm
    Out[12]:= Derivative[1][f][Power[y, 2]]

In SymPy we just cannot represent the object ``Derivative[1][f]`` because it
requires ``f`` to be an expression but instead ``f`` needs to be a head.
Instead we have awkward workarounds using ``Subs``:

.. doctest::

    >>> from sympy import Function, symbols, pprint
    >>> x, y = symbols('x, y')
    >>> f = Function('f')
    >>> f(x).diff(x).subs(x, y**2)
    Subs(Derivative(f(x), x), x, y**2)
    >>> pprint(_)  # doctest: +NORMALIZE_WHITESPACE
    /d       \|
    |--(f(x))||   2
    \dx      /|x=y

..
    --|

There are ways that this limitation could be worked around but it is not easy
and does not simply arise naturally. Likewise many users want to be able to
treat an undefined function as if it was a symbolic expression and do things
like:

.. doctest::

    >>> f(x).subs(f, cos)
    cos(x)

This barely works in SymPy and in fact really should just give an error. The
fact that it works at all is because well meaning people have put fudges into
the codebase that really just should not be there. The ``Function`` class in
fact needs to be a metaclass to be able to make things like this work which
really just shows that the problem is that ``f`` should not be a class. When we
call ``f = Function('f')`` above in fact what happens is that a class will be
created dynamically behind the scenes and then ``f`` will be that class. The
class itself would not behave in the way that people would expect of a symbol
and so it needs to have a metaclass to pretend that it is like a symbolic
expression even though it cannot be one. The problem is that ``f`` is a
*subclass* of ``Basic`` rather than an *instance* of ``Basic``.

Another problem with using classes for expression heads is that the expression
head does not necessarily represent the type of object that we have. As an
example consider an integral of a matrix:

.. doctest::

    >>> from sympy import Matrix, Integral
    >>> M = Matrix([[1, 2], [x, 4]])
    >>> M
    Matrix([
    [1, 2],
    [x, 4]])
    >>> expr = Integral(M, (x, 0, 1))
    >>> print(expr)
    Integral(Matrix([
    [1, 2],
    [x, 4]]), (x, 0, 1))

The expression tree for this expression looks like:

.. graphviz:: integral.dot

We see then that the head here is ``Integral``. However conceptually this
object is really a matrix. It should have the attributes and properties of a
matrix such as ``.shape``:

.. doctest::

    >>> M.shape
    (2, 2)
    >>> expr.shape
    Traceback (most recent call last):
      ...
    AttributeError: 'Integral' object has no attribute 'shape'
    >>> expr.doit().shape
    (2, 2)

Many other operations involving ``expr`` will not work properly like
``det(expr)``, ``expr.det()``, ``expr * M`` and so on. This basic problem
permeates much of SymPy. Objects are distinguished based on their Python type
which is their expression head. However the expression head does not really
represent the "type" of mathematical object we have or what the expression
represents. There are ways to work around this but it shows quite clearly that
there is a mismatch between identifying an expression head with the kind of
object and therefore that identifying expression heads with classes muddles
things up.


Object-oriented symbolic computing
----------------------------------

Many operations in SymPy are slow. Sometimes this slowness is inevitable just
because many symbolic algorithms are inherently slow. Most of the time though
it would be possible for SymPy to be much faster and this is an important
motivation for wanting to have something different from the current ``Basic``
design.

One suggestion for why SymPy can be slow is that it is written in Python rather
than say C++. This was part of the motivation for SymEngine which essentially
translates the basic design of SymPy's core into C++. The intention has always
been that it might be possible to use SymEngine as a replacement for the core
of SymPy while still using SymPy. While SymEngine is a lot faster than SymPy
the reasons for that are more complicated than just the difference between
Python and C++. SymEngine has a more limited scope than SymPy and has been
better optimised for particular common operations. SymEngine also makes
different choices for how to do things like numerical evaluation. Many of those
different choices would also make SymPy faster without rewriting anything in
C++. There are still problems with SymEngine's use of classes though in
particular if we eventually want SymEngine to be able to replace the core of
SymPy. Here the problem is extensibility: the fact that SymEngine uses C++
classes to define all of its behaviour means that it cannot be extended from
Python. That means that the only way SymEngine could replace the core of SymPy
is if pretty much *all* methods of *every* SymPy class were translated into
C++. And then after doing that it would not be possible to extend the behaviour
of any of those classes from within SymPy's Python codebase and to do all of
things that many downstream Python projects do.

The principle that SymPy could have a core written in something like C++
does make sense because this is how most things in Python work. The idea is
usually that Python is a nice language to work in for end users but you would
not usually write the computationally intensive parts of a widely used codebase
in Python itself. This is why SymPy is unusual in being a widely used
CPU-bound pure Python library. What would be better though than the SymEngine
approach of translating SymPy's classes into C++ is having a system based on
rules and pattern-matching where the engine that applies the rules can be
implemented as efficiently as possible but the rules themselves are still
controllable from the Python level. To make this work we would need to *not*
use classes for expression heads and *not* encode all behaviour in methods
defined on classes.

In any case the point that I really want to make here is that the idea that
SymPy is slow because it is written in Python is also something of a red
herring when it comes to thinking about how to make SymPy faster. It is
definitely possible to make something *much* faster than SymPy itself while
still working in pure Python. The design of ``Basic`` is poorly optimised
and in fact makes any attempt at optimisation very difficult. If SymPy had a
better optimised design in Python then it would be much easier to speed that
design up by rewriting some parts in another language if needed.

Currently an operation in SymPy such as ``solve`` involves executing the
methods of potentially hundreds of different classes and any one of those can
choose to do some slow operation if it wants. The authors of any one class will
have little idea of the implications that any operation can have for the
performance of higher-level algorithms such as differentiation and integration
etc. This is because in SymPy all of the logic that encodes the meanings of the
different kinds of expressions is encoded as executable code in methods. It is
not possible to call these methods without creating an instance of the class
and even just creating an instance can be slow because of all of the code that
executes during the constructor. This scattering of all of the important code
through hundreds of classes means that *no one* can tell you what the big-O
behaviour is for any significant operation that uses SymPy expressions.

In principle the advantage of using classes is that you can optimise each
individual class so that it uses the best data structure internally and then
users of the class do not need to care about how those internals work. The
problem with this is that then there can be only one possible choice of data
structure for each expression head. It is not possible for example to use one
kind of data structure for numerical evaluation and then a different kind of
data structure for differentiation because the choice of data structure is
encoded in the classes. What we need is a design that separates the
specification of expression heads from the choice of data structure so that we
can switch to different data structures for different operations. That means
not using classes for expression heads.


Different kinds of data structure
---------------------------------

If we were going to use different data structures then what would that look
like? Here we will consider some examples to motivate the idea of using
different data structures in different situations.

One of the models that is most commonly used in symbolic computing
systems is to replace the *tree* representation that I mentioned above with an
expression *graph* instead. Here the basic idea is that most large expressions
will have many repeating subexpressions. This happens because many operations
that create large expressions do so precisely because the operation ends up
repeating subexpressions in many places. A good example of this is the result
of applying the chain rule in differentiation:

.. doctest::

    >>> from sympy import cos, sin, exp, tan, diff
    >>> e = x*cos(sin(exp(tan(x))))
    >>> print(diff(e, x))
    -x*(tan(x)**2 + 1)*exp(tan(x))*sin(sin(exp(tan(x))))*cos(exp(tan(x))) + cos(sin(exp(tan(x))))

Note how many times ``tan(x)`` appears and how many times ``exp(tan(x))``
appears and so on. Differentiating a large expression almost always results in
many repeating expressions. As a tree this derivative expression looks like this:

.. graphviz:: diff.dot

The idea in representing the expression as a graph is that we never duplicate
any repeating expression. Every place in the expression that features the same
expression like ``tan(x)`` for example will just reuse the previously created
expression. Now our expression becomes a *directed acyclic graph* (DAG) and
its representation as a data structure is more like this:

.. graphviz:: diff_graph_2.dot

What we can see here is that this DAG has fewer nodes than the original tree
which will obviously save on memory. In practice a lot of the time this is what
happens in SymPy already because of the use of the cache. If an operation
repeatedly attempts to construct some expression then the cache will just
return new references to the already created expression e.g.

.. doctest::

    >>> Pow(x, 2) is Pow(x, 2)
    True

Here in Python the ``is`` operator verifies that these are the exact same
object in memory. The problem though is that the cache is just a cache and
cannot be relied upon. The cache will often reduce the memory usage of
repeating subexpressions but any operation in SymPy will still need to traverse
the expression as if it were the fully expanded tree. In SymPy if we wanted to
substitute a value like say ``x -> pi`` we would need to recurse down the whole
tree like this::

    def subs(expression, old_value, new_value):
        """Replace old_value with new_value in expression"""
        # Handle substitution
        if expression == old_value:
            return new_value

        # Base case for recursion is an atom
        old_args = expression.args
        if not old_args:
            return expression

        # Call subs recursively on all args
        new_args = [subs(arg, old_value, new_value) for arg in old_args]

        # Make a new expression with the new args
        return expression.func(*new_args)

Here the operation will traverse the entire tree processing every repeating
subexpression many times. Again we can use the cache to save some of the
repeated work. It would be much better though if we could actually evaluate
each subexpression exactly once rather than depending on the cache. One of the
problems with depending on the cache is that when we get to large expressions
that the cache itself gets full and suddenly the cost of all of the duplicated
work processing subexpressions becomes enormous. A better approach is to switch
from recursing down the tree to building up from the bottom. For that we want a
*topological sort* of the DAG which we will represent using a different data
structure.


Demonstration
-------------

Here now is some simple code that demonstrates all of the ideas above
