.. _data-structures-1:

Data structures for symbolic expressions: part 1
================================================

This will be the first of many blog posts discussing data structures and basic
algorithms for working with symbolic expressions. The reason for this is that I
want to change the way that SymPy's data structures work and I need to make
sure that lots of other people understand the design decisions that I am
considering.

Let us first use SymPy to consider what a data structure for a symbolic
expression might look like. In SymPy we can create symbolic expressions by
creating symbols and then calling symbolic functions on those symbols::

    >>> x = symbols('x')
