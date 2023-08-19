Towards a new SymPy: part 2 - Computational Algebra
===================================================

This is the second part of a series of posts on big changes for SymPy with
particular focus on speed. The other posts in this series can be found at
:ref:`new sympy`. The first post outlined what I consider to be the three
core subsystems of SymPy:

1. Robust numerical evaluation (``sympy.core.evalf``)
2. The symbolic subsystem (``sympy.core.basic``)
3. The computational algebra subsystem (``sympy.polys``)

In relation to the computational algebra subsystem, there are three ways in
which SymPy can be made faster:

1. Improve/extend some of the algorithms and features.
2. Make more use of the computational algebra subsystem in the rest of SymPy.
3. Make use of Flint which is a fast C library for computational algebra.

This post will focus on the first two points. I will write a separate post
about Flint and python-flint because I know that some people in other projects
will be more interested in using python-flint than SymPy and I hope to
encourage them to contribute to python-flint.

As before I am writing this with the intention that it should be understandable
to non SymPy developers. The primary intended audience though is other SymPy
developers because I want them to understand the significance of the work done
so far and the changes that I think are needed for the future.

Integers and rationals
----------------------

First, I will give a brief overview of the computational algebra subsystem. In
SymPy the computational algebra subsystem is located in the ``sympy.polys``
module. The lowest level of this system are the "domains" which are objects
that represent rings and fields. Briefly a ring is a set of "elements" that can
be added and multiplied and a field is a ring in which it is also possible to
divide any two elements (excepting 0). The most well known examples of each are
the ring of integers :math:`\mathbf{Z}` and the field of rational numbers
:math:`\mathbf{Q}`.

Engineers and scientists (rather than pure mathematicians) might wonder why we
need to distinguish between e.g. integers and rationals and worry about whether
or not division is "possible". In some sense we can always divide two integers
(excepting 0) and just get a rational number. Certainly in SymPy's *symbolic*
subsystem there is no difficulty in dividing two integers although we do need
to be careful to divide SymPy's ``Integer`` objects and not Python's ``int``
objects (whose division operator gives ``float``)::

    >>> from sympy import Integer
    >>> Integer(1)/Integer(2)  # dividing Integer gives Rational
    1/2
    >>> 1/2  # dividing ints gives floats
    0.5

The basic reasons that we need to distinguish carefully between
:math:`\mathbf{Z}` and :math:`\mathbf{Q}` are not particularly deep or
mathematical but rather purely *computational*. Efficient computing means
having well defined types, representations and operations. On a very basic
level the representation of an integer is a sequence of bits and the
representation of a rational number is a pair of integers. Arithmetic
operations manipulate these representations and so an efficient arithmetic
heavy algorithm needs to use a well defined representation for their elementary
types. If we choose our representation to be :math:`\mathbf{Z}` then we cannot
represent :math:`1/2` without converting to a different representation and so
division is not allowed. Of course we can always convert any element of
:math:`\mathbf{Z}` to an element of :math:`\mathbf{Q}` but this conversion is
not a free operation and the representation as an element of :math:`\mathbf{Q}`
would be less efficient if actually all of our elements are integers and we do
not need to divide them.

The other reason that we need to distinguish between e.g. :math:`\mathbf{Z}`
and :math:`\mathbf{Q}` in particular or more generally between different
"domains" is that we usually want to use different algorithms for different
domains. For example, the most efficient algorithm for inverting a matrix or
solving a linear system of equations over :math:`\mathbf{Z}` is different from
the most efficient algorithm for doing the same over :math:`\mathbf{Q}`. Often
the algorithms will convert to different domains for example to invert a matrix
of integers we might convert to :math:`\mathbf{Q}` and use a division-based
algorithm. Alternatively to invert a matrix of rational numbers we might factor
out denominators and convert to a matrix over :math:`\mathbf{Z}` before using a
fraction-free algorithm for inverting a matrix over :math:`\mathbf{Z}`.

In SymPy's domain system, there are domain objects like ``ZZ`` and ``QQ`` that
represent :math:`\mathbf{Z}` and :math:`\mathbf{Q}` respectively. These domain
objects are used to construct elements of the domain and to convert between
domains e.g.::

    >>> from sympy import ZZ, QQ
    >>> QQ(2,3)  # construct an element of QQ
    MPQ(2,3)
    >>> q = QQ(2)  # construct an integer-valued element of QQ
    >>> q
    MPQ(2,1)
    >>> z = ZZ.convert_from(q, QQ)  # convert from QQ to ZZ
    >>> z
    2
    >>> type(z)
    <class 'int'>

In this example we construct the rational number :math:`2/3`. Its
representation is the ``MPQ`` type which is a pure Python implementation of
rational numbers based on Python's ``int`` type (analogous to the Python
standard library ``fractions.Fraction`` type). The representation for ``ZZ``
uses Python's ``int`` type.

GMP and gmpy2
-------------

I mentioned the two domains ``ZZ`` and ``QQ`` above and I will describe the
other domains below. Since these two domains are absolutely foundational it is
worth discussing their implementation a bit first because it is important to
understand the role that gmpy2 and the underlying GMP C library play in this.
If gmpy2 is installed then SymPy will use its ``mpz`` and ``mpq`` types for
``ZZ`` and ``QQ`` respectively. These are based on the underlying GMP C
library::

    >>> from sympy import ZZ, QQ
    >>> ZZ(2)  # pip install gmpy2 to see this
    mpz(2)
    >>> QQ(2,3)
    mpq(2,3)

For small rational numbers gmpy2's ``mpq`` type is about 20x faster than
SymPy's pure Python ``MPQ`` type (when both are used from Python code). For
small integers ``mpz`` is about the same speed as Python's ``int`` type.
Potentially ``mpz`` is actually a little slower for small integers just because
CPython's ``int`` is heavily micro-optimised for very small integers. In Python
all integers use a single arbitrary precision type ``int`` which is unusual
among programming languages. Since most integers are small, CPython tries to
optimise for this case. Possibly also CPython can do some optimisations with
``int`` that it cannot do with "third party types" like gmpy2's ``mpz``.

For larger integers (bigger than machine precision) both ``mpz`` and ``int``
will represent an integer using multiple "digits" although the digits are
referred to as "limbs" in GMP. CPython uses 30-bit limbs because its
implementation of integer arithmetic is designed to be portable "generic C"
that can be compiled by any compiler. By contrast GMP's implementation is
designed to be as fast as possible at all costs and so it uses handcrafted
assembly code for many different CPU architectures. In CPU specific assembly it
is possible to access instructions like ``mulx`` or to read the CPU's carry
flag after an addition etc. These things are essential for being able to use
64-bit limbs but are not available in generic C. GMP's ``mpz`` type therefore
uses 64-bit limbs (as do all widely used fast big integer implementations).

As the bit size increases say to 1000 bits then gmpy2's ``mpz`` becomes a lot
faster than CPython's ``int``. At this bit size the difference in speed to
multiply two integers is something like 4x and is due to the fact that ``mpz``
uses 64-bit limbs whereas ``int`` uses 30-bit limbs. This smaller limb-size
essentially means a 4x increase in the number of CPU-level operations needed
for a big integer multiplication. At these bit sizes the fact that SymPy's
``MPQ`` is implemented in pure Python becomes less significant and what matters
more is just that SymPy's ``MPQ`` is built on top of CPython's slower ``int``
whereas the gmpy2/GMP ``mpq`` is built over GMP's faster ``mpz``.

For very large integers gmpy2's ``mpz`` is *much* faster than ``int``. This is
because GMP has a whole hierarchy of algorithms for different bit sizes.
CPython's ``int`` type tops out at the Karatsuba algorithm which is also used
by GMP for intermediate bit-sizes but GMP has more complex algorithms that are
used for larger bit sizes. At these large bit sizes the difference in speed
between ``mpz`` and ``int`` can be enormous. In fact the slowness of some
algorithms used by CPython's ``int`` type was even recently considered to be a
security vulnerability to the extent that certain operations were disabled:

https://discuss.python.org/t/int-str-conversions-broken-in-latest-python-bugfix-releases/18889

That did lead to some work on improving CPython's ``int`` algorithms but the
disabled operations remain disabled even though the algorithms have been
improved to some extent::

    >>> import math
    >>> x = math.factorial(1559)
    >>> print(x)
    ...
    ValueError: Exceeds the limit (4300 digits) for integer string conversion;
    use sys.set_int_max_str_digits() to increase the limit

What this all means is that for some operations SymPy is a lot faster when
gmpy2 is installed. Some people reading this might think that it seems absurd
to worry about the performance of megabyte sized integers but many symbolic
algorithms will generate much larger integers than might be expected. Also when
gmpy2 is installed it will be used by mpmath and so it speeds up SymPy's
numeric subsystem as well as the computational algebra subsystem.

The symbolic subsystem is the only one of the three core subsystems of SymPy
that does not use gmpy2. I did recently look at changing SymPy to have the
symbolic ``Integer`` and ``Rational`` types use ``mpz`` and ``mpq`` when gmpy2
is installed but hit a stumbling block that the symbolic subsystem allows
"unevaluated rationals"::

    >>> from sympy import Rational
    >>> Rational(3, 6, 1)
    3/6

That is apparently documented behaviour. This is not really a useful feature
though because the symbolic subsystem could also represent the same thing using
``Mul`` and ``Pow``. Having the symbolic ``Rational`` type use gmpy2's ``mpq``
for its internal representation would break these "unevaluated rationals" (that
should be done anyway though).

One of the nice things about the computational algebra subsystem is that it is
possible to swap out the implementation of the domain objects like this. So on
the one hand SymPy and its only hard dependency mpmath can be used as entirely
pure Python code without gmpy2. This is useful for many people who do not need
the performance of gmpy2 and do not want to install it or cannot install it for
example if they use a different Python implementation like PyPy. On the other
hand if gmpy2 is installed then SymPy will use it and will be a *lot* faster
for some operations but with *no other observable change in behaviour*.

(Contrast this last point with my previous comments about the difficulty of
SymPy using SymEngine to speed up the symbolic subsystem in :ref:`symengine`)

The domain system
-----------------

I have talked a lot about the two domains ``ZZ`` and ``QQ`` but there are many
more domains in SymPy. You can read more about them here:

https://docs.sympy.org/latest/modules/polys/domainsintro.html

People who are familiar with computational algebra will recognise these
domains:

- ``GF(n)``: integers mod ``n`` (the name ``GF`` is misleading)
- ``ZZ``: the integers
- ``QQ``: the rational numbers
- ``ZZ_I``: Gaussian integers
- ``QQ_I``: Gaussian rationals
- ``QQ(a)``: algebraic number field generated by ``a``
- ``RR``: the real numbers (floats with fixed precision provided by mpmath)
- ``CC``: the complex numbers (complex floats with fixed precision)
- ``K[x]``: polynomials in e.g. ``x`` with coefficients in another domain ``K``
- ``K[x,y]``: multivariate polynomials in ``x`` and ``y``.
- ``K(x,y)``: rational functions (ratios of polynomials) in ``x`` and ``y``.
- ``EX``: The expression domain (basically the symbolic subsystem)

There are more domains but these are the most important ones. Perhaps the
easiest way to see what the domains are for is by using the
``construct_domain`` function. This function is used internally by SymPy to
choose a domain that could represent an expression from the symbolic subsystem
in the domain system::

    >>> from sympy import *
    >>> x, y = symbols('x, y')
    >>> construct_domain([1, 2])
    (ZZ, [mpz(1), mpz(2)])

Here we asked for a domain that could represent both ``1`` and ``2``. What was
returned was the domain ``ZZ`` (meaning the integers) and a list of two
elements representing ``1`` and ``2`` in that domain (as gmpy2 ``mpz`` objects
in this case). We can try more examples::

    >>> construct_domain([x**2, 1])
    (ZZ[x], [x**2, 1])
    >>> construct_domain([x**2, y])
    (ZZ[x,y], [x**2, y])
    >>> construct_domain([x**2, y/x])
    (ZZ(x,y), [x**2, y/x])
    >>> construct_domain([sin(x), y])
    (ZZ[y,sin(x)], [(sin(x)), y])
    >>> construct_domain([x, 2.0])
    (RR[x], [x, 2.0])
    >>> construct_domain([x/2, 1])
    (QQ[x], [1/2*x, 1])

Sparse and dense polynomials
----------------------------

Importantly the polynomial domains are always implemented as "sparse"
polynomials. This means that only nonzero terms are stored. This example
contrasts the sparse and dense representations of polynomials::

    >>> from sympy import QQ, symbols
    >>> x, y = symbols('x, y')
    >>> e = x**10 + y
    >>> p_sparse = QQ[x,y].convert(e)
    >>> p_dense = e.as_poly()

This is how the sparse ``PolyElement`` and the dense ``Poly`` usually look::

    >>> p_sparse
    x**10 + y
    >>> p_dense
    Poly(x**10 + y, x, y, domain='ZZ')

This is what their internal representations look like::

    >>> dict(p_sparse) # internal sparse representation
    {(0, 1): mpq(1,1), (10, 0): mpq(1,1)}
    >>> e.as_poly().rep.rep  # internal dense representation
    [[mpz(1)], [], [], [], [], [], [], [], [], [], [mpz(1), mpz(0)]]

Here the sparse representation is a dictionary mapping exponent tuples to
coefficients. The dense representation is a list of lists of coefficients. The
empty lists in the dense representation represent zero terms. The dense
representation is described as "dense" because it needs to store zero terms
explicitly.

Integers mod ``n``
------------------

Some domains that people might expect to be find in the domain system are
missing like finite fields of non-prime order e.g. ``GF(2**3)``. Essentially
that is because most SymPy users are not interested in such things and the rest
of the codebase does not need them. Some people have expressed interest in
these and contributions are certainly welcome but I guess it has not happened
because it is not considered high priority and most users who are interested in
such things are more likely to use something like Sage rather than SymPy.

Probably most SymPy users are not interested in ``GF(n)`` (integers mod ``n``)
either but that is there because it is needed for the algorithms in the other
domains. Let me give an example to show how all of this is used::

    >>> from sympy import symbols, factor
    >>> x, y = symbols('x, y')
    >>> e = x**4 - y**4/16
    >>> e
    x**4 - y**4/16
    >>> factor(e)
    (2*x - y)*(2*x + y)*(4*x**2 + y**2)/16

So how does this work? First ``factor`` converts the expression to a polynomial
with coefficients in some domain::

    >>> p = e.as_poly()
    >>> p
    Poly(x**4 - 1/16*y**4, x, y, domain='QQ')
    >>> p.domain
    QQ

Here the ``Poly`` identifies that we have two variables ``x`` and ``y`` and
that the coefficients are in the domain ``QQ`` (the rational numbers). We are
now ready to call the factorisation algorithm. The factorisation algorithm for
polynomials with coefficients in ``QQ`` will first factor out the denominator
``16`` giving a polynomial ``16*x**4 - y**4`` with coefficients in ``ZZ``. We
now want to factorise this but then the algorithm for factorising polynomials
over ``ZZ`` will convert the problem to factorising polynomials over ``GF(p)``
for some prime ``p``. Then we compute the factorisation over ``GF(p)`` and
convert the result back to a factorisation over ``ZZ`` and so on. So the steps
in the computation are (with ``EX`` representing the ordinary symbolic
expressions)::

    EX -> QQ[x,y] -> ZZ[x,y] -> GF(p)[x,y] -> factored -> ... -> EX

The ``...`` here obscures a bunch of complexity that I don't want to get into.
For those familiar with these things the ``Zassenhaus`` algorithm is used by
default but the main weakness is that LLL-based techniques are not implemented
(so worst case is not polynomial time). For everyone else the algorithms used
here are usually good but for certain inputs ``factor`` can be very slow when a
different algorithm would be a lot faster.

The main point of this factorisation example is just to show the significance
of ``GF(p)``. Most SymPy users do not care about ``GF(p)`` but it is crucial
for things that they do care about because it is used by e.g. ``factor`` which
is in turn used by ``solve`` and ``simplify`` and so on.

Algebraic number fields
-----------------------

In this example we get ``EX`` which is what ``construct_domain`` returns when
it gives up::

    >>> construct_domain([sqrt(2), 1])
    (EX, [EX(sqrt(2)), EX(1)])

There is a domain for this but it will not be used by default (we have to pass
``extension=True``)::

    >>> construct_domain([sqrt(2), 1], extension=True)
    (QQ<sqrt(2)>, [ANP([mpq(1,1), mpq(0,1)], [mpq(1,1), mpq(0,1), mpq(-2,1)], QQ), ANP([mpq(1,1)], [mpq(1,1), mpq(0,1), mpq(-2,1)], QQ)])

It might not look nice but that is the domain for the algebraic number field
:math:`\mathbb{Q}(\sqrt{2})`. The ``ANP`` stands for "algebraic number
polynomial". The representation of algebraic number fields always uses a
*primitive element*. This representation is efficient for arithmetic but
computing the primitive element can be expensive which means that it can be
slow to construct the domain. Timings are (on a slow computer)::

    In [10]: %time ok = construct_domain([sqrt(2)], extension=True)
    CPU times: user 26 ms, sys: 0 ns, total: 26 ms
    Wall time: 25.2 ms

    In [11]: %time ok = construct_domain([sqrt(2), sqrt(3)], extension=True)
    CPU times: user 47.4 ms, sys: 0 ns, total: 47.4 ms
    Wall time: 46.3 ms

    In [12]: %time ok = construct_domain([sqrt(2), sqrt(3), sqrt(5)], extension=True)
    CPU times: user 55.2 ms, sys: 0 ns, total: 55.2 ms
    Wall time: 53.9 ms

    In [13]: %time ok = construct_domain([sqrt(2), sqrt(3), sqrt(5), sqrt(7)], extension=True)
    CPU times: user 120 ms, sys: 0 ns, total: 120 ms
    Wall time: 118 ms

    In [14]: %time ok = construct_domain([sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11)], extension=True)
    CPU times: user 688 ms, sys: 0 ns, total: 688 ms
    Wall time: 686 ms

    In [15]: %time ok = construct_domain([sqrt(2), sqrt(3), sqrt(5), sqrt(7), sqrt(11), sqrt(13)], extens
        ...: ion=True)
    ^C^C
    KeyboardInterrupt

I don't know how long that last command would take but I interrupted it after
about 5 minutes. I have not investigated why it is so slow but I expect that it
can be made faster. I think what it really shows though is that it is a bad
idea to even try to compute the primitive element and that it is better to
represent algebraic number fields differently in the case of having many
algebraic generators.

