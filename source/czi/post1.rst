Towards a new SymPy: part 1 - Outline
=====================================

This is the first part of a series of posts on big changes for SymPy with
particular focus on speed. What I am going to describe in these posts are:

- A description of the of the issues around speed with SymPy currently.
- Work that I (and others) have done to speed up SymPy after the past few years.
- Work in particular that I have done over the last year to make SymPy faster.
- What I think are the next steps to make SymPy faster.

The other posts in this series can be found at :ref:`new sympy`.

Over the last year in particular I have been working as a part of a CZI funded
project that has three strands. One of those three strands is for me to work on
speeding up SymPy. Now that we come to the end of that year, I want to describe
what has been done and spell out a vision for the future.

I will be writing this in a series of blog posts. This first post will outline
the structure of the foundations of a computer algebra system (CAS) like SymPy,
describe some problems SymPy currently has and what can be done to address
them. Then subsequent posts will focus in more detail on particular components
and the work that has been done and what should be done in the future.

I am writing this with the intention that it should be accessible to someone
who is not a SymPy developer although other SymPy developers are the intended
audience for many of the points that I will make. Many of the things that I
will describe here are not well understood even by many SymPy developers though
and a major goal of this series of posts is to help to try to change that.

What is a CAS?
--------------

There are many different kinds of computer algebra system (CAS) and they have
different intended uses that are driven by their different kinds of authors,
contributors and users. Some like SageMath are directed more at pure
mathematicians, while others like SymPy are directed more at scientists and
engineers.

Most scientific computing is not done with computer algebra systems. Instead
numerical libraries (e.g. NumPy, SciPy, etc in the context of Python) are the
usual backbone of computing in science and engineering. What those libraries do
is to provide fast implementations of numerical computing based on machine
precision floating point numbers. When I say machine precision floating point I
mean something like a 64-bit binary floating point number like the ``float``
type in Python:

    >>> from math import pi
    >>> pi
    3.141592653589793
    >>> type(pi)
    <class 'float'>

Many things can be done with machine precision floating point numbers, but
there are situations where they are not good enough. Programmers learn fairly
early on that floating point calculations are not exact e.g.:

    >>> (1 / 49) * 49
    0.9999999999999999

We can see here the limits of the precision of 64-bit floating point because we
should have a result of 1 but in fact we have something that is quite close to
1 but still not exactly 1. Since the answer is quite close most of the time
this is usually not a huge problem even if it does cause some confusion for
beginners.

However, there are many situations where this is not acceptable:

- Small errors can accumulate and become large errors.
- To solve the error accumulation problem some calculations should be done with
  many more digits of precision than machine precision floating point numbers
  can allow.
- Some calculations cannot be done reliably at all unless they are done
  *exactly*.
- Sometimes we want to calculate not in terms of explicit numbers but in terms
  of symbols to obtain e.g. "analytic" rather than "numerical" results.

Traditionally these other situations are handled by scientists and engineers
turning to a CAS and so in the context of scientific computing the usual
distinction is between working "numerically" with floating point numbers or
otherwise using a CAS to work "symbolically" or "analytically". There is a lot
more to a CAS than this but in simple terms this is the role that SymPy plays
within the Python scientific computing ecosystem. Of course there are many
people also using SymPy for applications that are more like pure mathematics
but most SymPy users and developers are from a scientific background and that
has a big influence on the way that SymPy is developed.

The "core" of a CAS
-------------------

Different people will have different ideas about what the "core" of a CAS like
SymPy is. For example users might think of functions like ``solve`` or
``integrate`` as being the core functions but these are really high level
functions and not really the foundation of a CAS. Most SymPy developers would
probably think of the ``sympy.core`` module as being the core of SymPy but this
actually misses a major component of SymPy which is the ``sympy.polys`` module.
In my mind these three components will be the core of any CAS like SymPy:

- A robust numerical evaluation subsystem (``sympy.core.evalf``).
- A symbolic manipulation subsystem (``sympy.core.basic``).
- A computational algebra subsystem (``sympy.polys``).

I will explain each of these in turn because it is easy to miss the importance
of each of them and I think that many SymPy developers do not understand the
role each of these components has and how they relate to one another. Of course
there are many other components of SymPy that are important but these three are
the foundation of almost all of the higher-level things that SymPy does.
Importantly when we want to think about the speed of SymPy we need to
understand that almost everything depends on the speed of these three
components. Higher level algorithms are absolutely important but at least in
SymPy those are all built on top of these three pieces as the foundation for
everything else.

A robust numerical evaluation subsystem
---------------------------------------

Above I use the word "numerical" to mean calculations that are done with
machine precision floating point numbers. This is the usual intended meaning of
the word but this use of the word is tricky because in the context of a CAS we
can also do numerical calculations with arbitrary precision approximate or
exact numbers and those are also "numerical" calculations. I will define
several levels of numerical calculation ranging from machine precision floating
point to exact and symbolic calculations. As we go down this list the
calculations become slower but more accurate:

- Machine precision floating point: This is the usual 64-bit binary floating
  point numbers that are used in most numerical libraries. They are fast and
  have a wide range of support in hardware and software but are limited in
  both range and precision. Most scientific Python libraries like NumPy, SciPy
  etc work (almost) exclusively with machine precision floating point numbers.
- Arbitrary *fixed* precision floating point: This means using floating point
  numbers that can have as many digits of precision as desired. These
  calculations are much slower than machine precision floating point
  calculations but by increasing the precision we can make the calculations as
  accurate as we want. The ``mpmath`` library that is used by SymPy is an
  example of this. With ``mpmath`` you can choose to use say 100 digits and
  then all calculations will be done with 100 digit numbers.
- "Robust" floating point calculations: These are calculations that are done
  with arbitrary precision floating point numbers but are done in such a way
  that the result is guaranteed (perhaps formally or perhaps heuristically) to
  be correct to within a specified tolerance. Making this work requires using
  *adaptive* precision so if the end result should be accurate to 100 digits
  then intermediate calculations might need to use many more than 100 digits to
  control the growth of errors. SymPy's ``evalf`` is an example of *heuristic*
  robust numerics and the Arb library that I will discuss later is an example
  of a library that does *formal* robust numerics.
- Exact numerical calculations: these are numerical calculations that are done
  using e.g. exact numbers like rational numbers or more complicated types of
  exact numbers.
- Symbolic calculations: these are calculations that are done with symbols
  rather than just numbers. The result of a symbolic calculation is a symbolic
  expression that can be manipulated further. This is the kind of calculation
  that probably most people think of when they think of a CAS.

Let us just quickly demonstrate each of these:

- Machine precision floating point:

    >>> from math import sqrt
    >>> sqrt(2)  # Use machine precision floating point
    1.4142135623730951

- Arbitrary precision floating point:

    >>> from mpmath import mp
    >>> mp.dps = 50 # Use 50 digits for the calculation
    >>> mp.sqrt(2)
    mpf('1.4142135623730950488016887242096980785696718753769')

- Robust floating point calculations:

    >>> from sympy import sqrt
    >>> sqrt(2).evalf(50) # Compute a result that is correct to 50 digits
    1.4142135623730950488016887242096980785696718753769

- Exact numerical calculations:

    >>> from sympy import Rational
    >>> Rational(1, 2) + Rational(1, 3) # Compute exactly
    5/6

- Symbolic calculations:

    >>> from sympy import Symbol, sqrt
    >>> x = Symbol('x')
    >>> sqrt(x)**2   # Compute with symbols
    x

The distinction between these different kinds of calculations can be a bit
fuzzy but the first point to note is that the vast majority of scientific
computing is done with machine precision floating point numbers as mentioned at
the top of the list. Everything below that is what a CAS like SymPy is
typically used for. It is also possible to do symbolic calculations with
e.g. machine precision floating point numbers so perhaps including "symbolic"
in this list does not make sense but I think that from a user's perspective
"symbolic" is in some way the next level after exact numerical calculations.

Many SymPy users do not understand the distinctions between these different
kinds of numeric calculations but also many SymPy developers do not understand
them either. For example I am not sure how many SymPy developers would
automatically realise ``mpmath`` and ``evalf`` are at different levels in this
scheme despite the fact that ``evalf`` works entirely by using ``mpmath``. As
an example to demonstrate the differences we will compute :math:`\sin(1)^2 +
\cos(1)^2 - 1` using SymPy's ``evalf``, mpmath and also SymEngine which is a
C++ library that recreates the "core" of SymPy in C++:

    >>> import sympy as sym
    >>> e = sym.cos(1)**2 + sym.sin(1)**2 - 1
    >>> e.evalf()
    -0.e-124

    >>> from mpmath import mp
    >>> mp.dps = 100
    >>> mp.cos(1)**2 + mp.sin(1)**2 - 1
    mpf('0.0')

    >>> import symengine as se
    >>> (se.cos(1)**2 + se.sin(1)**2 - 1).evalf()
    1.11022302462516e-16

The correct answer here is zero because :math:`\sin^2(x) + \cos^2(x) = 1` for
all :math:`x`. SymPy's ``evalf`` method returns the strange looking result
``-0.e-124``. This is ``evalf``'s way of saying that the result is smaller than
:math:`10^{-124}` and is possibly zero but not proven to be zero. The robust
numerics in ``evalf`` should ensure that we know that the result is very close
to zero but can still never prove that it is or is not exactly zero. We do not
get an exact result of zero because ``evalf`` is careful not to claim that the
result is exactly zero if it cannot prove that. Usually ``evalf`` will prove
that a nonzero expression is nonzero but it can never prove that a zero
expression is exactly zero. This is a fundamental limitation of robust (rather
than exact) numerics: we can never prove that something is *exactly* zero. In
order to compute this result ``evalf`` will have used mpmath but with
higher and higher precision until it reached the maximum allowed precision (124
digits) and gave up.

The mpmath result is exactly zero which is exactly correct but it is only
*exactly* correct by "luck". Here mpmath is using 100 decimal digits but it
still does not *prove* that the result is zero. It does however return
something that is indistinguishable from an exact zero unlike ``evalf`` which
made it clear that it could not prove that the result is exactly zero. This is
a fundamental limitation of arbitrary *fixed* precision numerics as compared to
robust numerics: we can use as many digits as we like and the result will
likely be accurate but without very careful analysis we generally do not know
*anything* about how accurate the result is. It is always possible that had we
used more digits we would have found that the result was not zero after all but
mpmath does not indicate any level of uncertainty about the result it returns.

By contrast SymEngine's ``evalf`` has computed the result here using
machine precision floating point numbers and so it gives a fast result but the
result is not zero and there is no way to know if that is due to a rounding
error or not. We can ask SymEngine's ``evalf`` to use more digits but it
will still only use arbitrary *fixed* precision (like ``mpmath``) and not
*robust* numerics (like SymPy's ``evalf``):

    >>> (se.cos(1)**2 + se.sin(1)**2 - 1).evalf(100)
    -3.9443045261050590270586428264e-31

There is nothing inherently wrong with SymEngine's ``evalf`` and there are good
reasons to use all of the different levels of numeric calculation that I have
described above. Many SymPy users would actually be happier if SymPy's
``evalf`` was faster even if less accurate and would prefer the behaviour of
SymEngine's ``evalf``. However these levels of numeric calculation are not
interchangeable and SymPy's other two core systems (the expression system and
the computational algebra subsystem) absolutely *do* depend on SymPy's
``evalf`` giving *robust* numeric evaluation and not just arbitrary *fixed*
precision numeric evaluation: many things would break if SymPy's ``evalf`` was
changed to behave like SymEngine's ``evalf``.

I will propose later a plan for how to improve SymPy's numeric evaluation
subsystem both in terms of speed and accuracy. For now I just want to note that
the original author of both mpmath and SymPy's ``evalf`` is Fredrik Johansson
who subsequently went on to create Arb which is a library for doing *formal*
robust numerics. SymPy's ``evalf`` should change to using Arb-like formal
robust numerics and should ultimately provide the option to use the Arb library
as the basis for this subsystem.

The symbolic expression system
------------------------------

SymPy's symbolic expression system is the second component of what I call the
"core" of SymPy and is located in the ``sympy.core`` package. This is what most
SymPy users and contributors are used to working with. This system defines
expressions in a symbolic tree representation e.g.:

    >>> import sympy as sym
    >>> x = sym.Symbol("x")
    >>> e = x**2 + 2*x + 1
    >>> e
    x**2 + 2*x + 1
    >>> sym.srepr(e)
    "Add(Pow(Symbol('x'), Integer(2)), Mul(Integer(2), Symbol('x')), Integer(1))"
    >>> sym.print_tree(e, assumptions=False)
    Add: x**2 + 2*x + 1
    +-One: 1
    +-Pow: x**2
    | +-Symbol: x
    | +-Integer: 2
    +-Mul: 2*x
      +-Integer: 2
      +-Symbol: x

Most work that goes on in SymPy is done either on the internals of this
expression system or on the other code that operates with these expressions. In
many ways this system is nice but there are also many problems with it.
Essentially it is designed with an emphasis on what would be nice for end users
who are doing simple things and as a result is not well suited as a foundation
for building more complex algorithms. What we do not have is any alternative
that can be used instead of this system for the internals of SymPy.

Many of the problems with SymPy and in particular its performance stem from
overuse of this expression system. Unfortunately the prominent exposure of the
symbolic subsystem in the SymPy API makes it very difficult to change and so
realistically the best path forward is to reduce the usage of this system at
least in the internals of SymPy. That is difficult though because we usually do
not have a clear alternative to use instead and most SymPy developers do not
understand how to use SymPy apart from by using this system.

When I say that the symbolic subsystem is overused I should be clear about
what the alternatives would be:

- Use a symbolic subsystem that has a very different design.
- Use the computational algebra subsystem instead.

Much of the work that I have done recently in SymPy has been to try to expand
the computational algebra subsystem, to make it faster and to make more use of
it for heavier algorithms in the internals of things like ``solve``,
``integrate`` etc. Many things in SymPy (e.g. matrices) can almost immediately
be made a lot faster simply by having them use the computational algebra
subsystem instead of the symbolic subsystem. I will talk more about this later.

Of course being a CAS that is primarly intended for symbolics there are many
things in SymPy that do need to use a symbolic subsystem but the current
design of core symbolics in SymPy is not suitable for most of those things.
The problems manifest both in terms of:

- speed: many things are much slower than they could be.
- behaviour: many things being more difficult to do.
- features: many users want to do things that the design cannot really support.
- bugs: the system is hard to use robustly and this leads to bugs.
- maintainability: making changes to any part of the expression system
  (including just fixing obvious bugs) is difficult because the effects of any
  change are far reaching and unpredictable.

.. _symengine:

What about SymEngine?
~~~~~~~~~~~~~~~~~~~~~

In terms of speed one approach that has been mooted for making the symbolic
subsystem faster is to rewrite the "core" in e.g. C++ and this is essentially
what SymEngine is. So one possibility to make SymPy faster would be to use
SymEngine instead of SymPy's symbolic subsystem. However SymEngine is not a
drop-in replacement for SymPy's symbolic subsystem and it could never be made
to be so. There are two basic problems with attempting to use SymEngine for the
core of SymPy:

- On the one hand SymEngine is not sufficiently similar to SymPy's existing
  symbolic subsystem so simply switching to use it inside SymPy would break all
  sorts of things.
- On the other hand SymEngine is too much like SymPy's symbolic subsystem so
  using it would not solve many of the problems that SymPy has and in fact
  would make it *much harder* to solve those problems.

The part of SymPy that SymEngine aims to emulate is literally made up hundreds
of thousands of lines of code and has poorly defined semantics in many cases.
It is basically impossible to make any drop-in replacement for this system that
would not differ in ways that would break things. At the same time some of the
SymEngine interface is *deliberately* different from SymPy (e.g. ``evalf``
above) so that changing SymPy's behaviour to be like SymEngine would break
compatibility for users and downstream projects.

It is also not possible to adapt the behaviour of SymEngine to be more like
what SymPy needs because it is not extensible from Python. Previously Julia
used both SymPy and SymEngine for symbolics but subsequently they decided to
move away from both and create ``SymbolicUtils.jl`` and ``Symbolics.jl`` to be
the foundation for symbolics in Julia. The reasons given at the time that
neither SymPy nor SymEngine were suitable were that:

- SymPy worked for what they needed but was too slow.
- SymEngine was fast enough but too inflexible and not extensible from Julia.

These same two considerations actually apply to SymPy itself: SymPy's symbolic
subsystem is too slow for SymPy itself and SymEngine is too inflexible for
SymPy to use it as a replacement.

There is a third problem which the Julia people seemed to overlook which is
that *both* SymPy's symbolic subsystem *and* SymEngine's reimplementation of it
are based on a design that is not suitable for the kinds of things that SymPy
needs either for its internals or also for what many users and downstream
libraries would like. There are many aspects to this design problem but the
most fundamental one is the problem of automatic evaluation. What I mean by
automatic evaluation is basically this:

    >>> import sympy as sym
    >>> sym.cos(sym.pi/4)
    sqrt(2)/2

Quite simply you asked to create the expression :math:`\cos(\pi/4)` and SymPy
instead gave you the expression :math:`\sqrt{2}/2`. Maybe that is what you
wanted but maybe it is not. The problem is that you do not have any real
control over this. There is ``evaluate=False`` but it does not work in general
and it cannot be made to work in general.

In terms of speed the problem with automatic evaluation is that it makes it
impossible to control the performance of higher-level algorithms because every
time any expression is created a huge amount of computational work is done in
order to try to "evaluate" the expression. Every now and again someone will
ask a question like "what is the computational complexity of <some SymPy
function>" but if this function uses the symbolic subsystem then it is
entirely impossible to answer this question: unbounded computation can occur
just during the creation of an expression.

We can try to reduce the cost of the computational work during automatic
evaluation but actually we really need to be able to reduce it to *zero* which
is what it would be if there were *no* automatic evaluation. Creating an
expression needs to be so cheap that we think of it as free compared to doing
any actual computation (many SymPy contributors *do* think of it as free which
leads them to write very inefficient code).

Many other problems with the design of SymPy's symbolic subsystem could in
principle be fixed without most users really noticing the internal changes.
Almost all code that uses the symbolic subsystem relies on automatic evaluation
though and so simply making a change to disable evaluation would break almost
everything.

The reason this makes it difficult for SymPy to use SymEngine is that SymEngine
is *also* based on automatic evaluation and in such a way that is impossible to
control from the outside. In the case of SymEngine ``evaluate=False`` is not
an option and I believe there are even plans for SymEngine's internal data
structures to change so that it would be *impossible* to implement something
like ``evaluate=False``. SymEngine broadly does "less" automatic evaluation
than SymPy which is one reason why it is faster but there is still some. Even
just swapping the order of the terms here is unacceptable if there is no way to
disable it:

    >>> import symengine as se
    >>> x = se.symbols('x')
    >>> se.exp(x) + x
    x + exp(x)

A new symbolic subsystem
~~~~~~~~~~~~~~~~~~~~~~~~

What SymPy needs is a symbolic subsystem that is *not* based on automatic
evaluation. Of course many users would still want automatic evaluation and
there should be a way to provide that. It is *easy* though to build a system
that has automatic evaluation on top of a system that does not whereas it is
*impossible* to build a system that does not have automatic evaluation on top
of a system that does. The current design of ``Basic`` and ``Expr`` in SymPy
and all of their hundreds of subclasses builds automatic evaluation into the
core data structures of the symbolic subsystem. This is a fundamental design
flaw that causes all kinds of problems for speed, behaviour and extensibility.
Using SymEngine (in its current form) can make this subsystem faster but it
would then make all of the *other* problems unsolvable.

The solution to this problem is to build a new symbolic subsystem but it really
needs to be built from the ground up: there is no viable way to get there
through incremental changes to the existing system. This can be done for the
internals of SymPy to get a lot of the benefit in terms of speed and behaviour
without breaking compatibility for users and downstream projects (which could
still use the existing system).

At some point though it would be better to switch the user facing symbolic
system to an implementation that would be based internally on that new
subsystem. We can try to minimise the noticeable impact of this change but it
would definitely be a breaking change. There has been some discussion of a
SymPy 2.0 that would make some backwards incompatible changes and this is the
number one thing that should be done. In my opinion we are not ready for SymPy
2.0 until we are ready to switch out the internals of the symbolic subsystem. I
don't see any other potential change that is important enough to warrant any
major compatibility break.

Of course if SymPy were redesigned to have a hypothetical new symbolic
subsystem with a different design then either SymEngine or something similar
(let's call it SymEngineX) could be designed to provide a faster implementation
of that subsystem. The new system though should be designed so that it is
possible for a faster implementation to be created without all of the problems
that would apply to SymPy using SymEngine as it is now. In particular:

- The parts that would be provided by SymEngineX must be a *small* part of the
  code that makes up the whole system (rather than a rewrite of *everything*).
- Those parts would need to have well defined semantics such that it is
  actually possible to make multiple compatible implementations (at least the
  SymEngineX implementation would need to be compatible with the pure Python
  SymPy implementation).
- The behaviour of the overall system still needs to be controllable from
  Python so *none* of the evaluation etc rules of SymEngineX could be
  hard-coded (in the way that SymPy and SymEngine both hard-code everything
  right now).

If a new symbolic subsystem was designed in this way then it would become much
easier to make an alternate core that would play the role that SymEngine was
intended to play in relation to SymPy. It could be SymEngine itself that still
provides that core but the parts that SymPy would need from it would look very
different to what SymEngine provides right now.

I want to be clear here that I do not want to criticise the excellent work done
by many contributors to both SymPy and SymEngine working on these symbolic
systems. There is a *lot* of great code in both projects that make up these
systems but there are also problems that are just not fixable without a
*complete* rebuild of the internals starting from the foundations. The existing
code and contributions would not be wasted but much of that code would need to
be adapted over time to a new framework.

I will talk about my proposals for the symbolic subsystem more concretely in a
separate blog post.

The computational algebra subsystem
-----------------------------------

SymPy's computational algebra subsystem is the third component of what I call
the "core" of SymPy and is located in the ``sympy.polys`` package. This
subsystem is underused and underappreciated but is absolutely critical to
much SymPy's ability to do anything nontrivial. In a CAS that was aimed more at
pure mathematics this subsystem would be the most prominent feature but in
SymPy it is mostly hidden away and not well understood. Few SymPy users or
developers know how to use the lower levels of this system and few are aware of
the features that it provides or why you might want to use them instead of the
"higher level" symbolic subsystem.

To give a clear example of how this subsystem is used I will show how to
convert an expression from the symbolic subsystem into the computational
algebra subsystem using the high-level ``Poly`` representation:

    >>> import sympy as sym
    >>> x = sym.Symbol("x")
    >>> e = (x + 1)**2
    >>> p = sym.Poly(e, x)
    >>> p
    Poly(x**2 + 2*x + 1, x, domain='ZZ')
    >>> p.as_expr()
    x**2 + 2*x + 1

The ``Poly`` representation is the highest level representation in the
computational algebra subsystem but its internals are still very different from
the symbolic subsystem. We can see this by looking at the internal
representation:

    >>> p.rep.rep
    [1, 2, 1]

This is a representation of the polynomial as a list of coefficients. Within
the computational algebra subsystem this is referred to as the dense univariate
polynomial (DUP) representation. In this representation we cannot represent
the unexpanded power ``(x + 1)**2`` but instead we have to expand it out to
``x**2 + 2*x + 1`` because we can only store the coefficients of the expanded
form. This is a very simple example but this computational algebra subsystem
can represent much more complex expressions. I have described this subsystem in
more detail in this page of the SymPy documentation:

https://docs.sympy.org/latest/modules/polys/domainsintro.html

Many important SymPy functions like ``factor``, ``cancel`` etc work by
converting from the symbolic subsystem to the computational algebra subsystem,
doing some computation and then converting back to the symbolic subsystem.
These conversions back and forth between the two systems are quite expensive
and often just converting back to the symbolic subsystem can be more expensive
than the actual calculation that is done in the computational algebra system.
This is because just *creating* expressions is slow in the symbolic system
(mainly because of automatic evaluation). It is much more efficient to perform
an entire calculation in the computational algebra subsystem without converting
back and forth but any code that does this needs to be written by someone who
has some understanding of how to use the computational algebra subsystem. At
least with the conversions someone can write code that seems like it only uses
the symbolic subsystem without needing to know anything about the computational
algebra subsystem. The price for that convenience is that it makes various
things in SymPy slower than they need to be.

So many things could be made faster if SymPy just used the computational
algebra subsystem more. The problem is that a certain mathematical background
is needed just to understand how to use it and most of the people who use and
contribute to SymPy do not have that background. This is also not helped by the
fact that the system is not well documented: there is a lot of documentation
but mostly it has not been written so that someone who does not have a
background in computational algebra could understand even what the basic
functions do (I wrote the doc page linked above to try to address this
partially).

Speaking of documentation I think it is quite telling about both SymPy users
and developers that while here I am describing the ``sympy.polys`` module as
one of the three pillars in the foundation of all of SymPy the docs do not even
bother to mention it at top level (as of SymPy 1.12):

https://docs.sympy.org/latest/reference/index.html#reference

The main front-page of the API docs consider that "polynomials" should just be
listed under "topics" unlike ``ntheory`` (number theory) which is a whole
section to itself. The tutorial does not mention the ``polys`` module directly,
there is nothing in the "guides" section about it etc. To give some context
here is a very rough count of how many lines of code are in each of SymPy's
top-level modules (this includes tests, comments, docstrings etc)::

    polys            100651
    physics           80805
    core              58583
    printing          48638
    solvers           43966
    functions         42605
    matrices          32942
    utilities         32248
    combinatorics     26144
    integrals         25588
    parsing           25213
    stats             22356
    tensor            19663
    simplify          18296
    geometry          15534
    assumptions       11524
    ntheory           11464
    series            11258
    plotting          11034
    sets              10400
    concrete           7093
    logic              6929
    vector             6772
    codegen            6489
    categories         4730
    testing            4610
    holonomic          4299
    crypto             3958
    calculus           3665
    diffgeom           3046
    liealgebras        2189
    external           2182
    algebras           2086
    discrete           1813
    interactive        1416
    strategies         1406
    multipledispatch   1239

The other things listed as "topics" alongside ``polys`` are ``geometry``,
``holonomic``, ``liealgebras``, ``categories``, ``crypto``, ``diffgeom``,
``plotting``, and ``stats``. Some of these are barely used by anyone while
others like ``plotting`` are widely used and should also not be listed here.
The "polynomials" module is the largest module in SymPy making up over 10% of
the all of the code and is used by *everything* else but is sorely undervalued
by the SymPy community.

The basic problem here is that SymPy is a project that is mostly used and
developed by people who are not pure mathematicians and the ``polys`` module
does not look like something that would be relevant to them. We need to
understand though that a CAS needs to have a computational algebra subsystem in
order to do the things that scientists and engineers want it to do as well.
Higher level functions that like ``solve``, ``integrate`` etc that are used by
scientists and engineers need to be built on top of a computational algebra
subsystem in order to be fast and reliable.

Maintenance of the computational algebra subsystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When I started working on SymPy I did not have a background in computational
algebra and I also did not understand how to use this subsystem (the ``polys``
module). Over time I have learned more about it and I have also learned how to
use it and I have made various improvements to it. To reach this point I had to
read a lot of textbooks and papers and look at what other (more
algebra-oriented) CASs do.

To begin with what I found was that there was this large almost unmaintained
subsystem that was critical to SymPy but yet did not have any active maintainer
working on it. Kalevi Suominen was the only SymPy developer who seemed to
really understand this part of the codebase and was critical to ensuring that
at least pull requests could be reviewed but was also not actively working on
major improvements or any redesign. Kalevi's guidance made it possible for me
to make improvements to this subsystem while in the process of learning about
it.

The code in the ``polys`` module is well designed. Of course there can always
be improvements but the fundamentals of the design are good. Looking back to
before my own involvement in SymPy it looks like it was originally developed by
Mateusz Paprocki and others. It is clear from the code that it was at some
point well maintained and loved by people who understood it well and cared
about it and had some vision for how it should work and what the plan was for
the future. Then it looks though as if those people at some point just moved on
leaving the code with various things in a sort of work-in-progress state or in
the middle of a redesign that was never finished. This is probably the hardest
part about trying to work on it because there are lots of pieces there but no
one remembers what the plan was supposed to be for the future.

More recently though there has been renewed interest in this subsystem and
there are several maintainers who can understand the code and make
improvements. The most significant improvements recently are

- The addition of the ``DomainMatrix`` class which has been mostly work done by
  myself (I will talk about this more later).
- The significant redesign and expansion of the ``sympy.polys.numberfields``
  module by Steve Keiffer who has also added some significant algorithms and
  improved others in parts of the code that had been untouched for years.

There is another contribution that I want to draw attention to that is very
relevant which is this PR by another maintainer Sangyub Lee:

    https://github.com/sympy/sympy/pull/20614

The significance of this PR is that it transfers the calculation of a
particular operation (here the ``Matrix.eigenvects`` method) from the symbolic
subsystem to the computational algebra subsystem. This is a very good example
of exactly what needs to be done more in SymPy and how things that many SymPy
users could relate to (eigenvectors in this example) can be computed much
faster and better simply by using the computational algebra subsystem more.

The main point here though is just that the ``polys`` module is not
(almost-)unmaintained any more. There are multiple maintainers who can
understand it and make improvements to it and we are now in a better position
to improve its capabilities and expand its use within SymPy.

Improvements to the computational algebra subsystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Above I casually suggested that the numeric subsystem should be completely
rewritten to use Arb-like semantics. I then went on to say that the symbolic
subsystem should be rewritten (which is a major undertaking for SymPy!). I
would not say the same thing about the computational algebra subsystem:

- The fundamentals of the design are good.
- Much of the code is written exactly how it should be written and is well
  micro-optimized.
- The more widely used parts are mature and well tested.

Of course there can be improvements like adding more features, improving
algorithms and so on, but this system can be improved incrementally. There is
however a major limitation here. While I say that much of the code is written
exactly how it should be written, it is all Python code and the lower level
parts of this subsystem are precisely the sort of thing that should be written
in a lower level language like C. It is clear from the code that this was
always part of the development plan but then it just didn't happen. Happily we
do not need to rewrite everything in C ourselves because there is now already a
C library that does exactly what we need: Flint. I have been working on making
it possible for SymPy to leverage Flint which is a library for doing fast
computations with polynomials, matrices and other things in the same way that
SymPy's computational algebra subsystem does.

As for the other topics I mentioned above, I will talk about this more in a
subsequent post.

Summary
-------

Above I have described what I consider to be the subsystems that make up the
core of SymPy and what the current state of each of them is, as well as what I
think needs to be done to improve them. In subsequent posts I will talk about
each of these in more detail to explain what has been done and what should be
done going forward. Briefly though:

- The numeric subsystem should be rewritten to use Arb-like semantics and
  should be able to use Arb directly at least as an optional backend.
- The symbolic subsystem should be rebuilt from the ground up based on a
  non-evaluating representation of expressions.
- The computational algebra subsystem should be made faster both by improving
  its algorithms and also by leveraging Flint.
- The computational algebra subsystem should be used more by the rest of SymPy.
