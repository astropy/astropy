*********************************************
Decorating functions to accept NDData objects
*********************************************

.. important:: The functionality described here is still experimental and will
               likely evolve over time as more packages make use of it.

Introduction
============

The `astropy.nddata` module includes a decorator
:func:`~astropy.nddata.support_nddata` that makes it easy for developers and
users to write functions that can accept either :class:`~astropy.nddata.NDData`
objects and also separate arguments.

Getting started
===============

Let's consider the following function::

    def test(data, wcs=None, unit=None, n_iterations=3):
        ...

Now let's say that we want to be able to call the function as ``test(nd)``
where ``nd`` is a :class:`~astropy.nddata.NDData` instance. We can decorate
this function using :func:`~astropy.nddata.support_nddata`::

    from astropy.nddata import support_nddata

    @support_nddata
    def test(data, wcs=None, unit=None, n_iterations=3):
        ...

which makes it so that when the user calls ``test(nd)``, the function would
automatically be called with::

    test(nd.data, wcs=nd.wcs, unit=nd.unit)

That is, the decorator looks at the signature of the function and checks if any
of the arguments are also properties of the ``NDData`` object, and passes them
as individual arguments. The function can also be called with separate
arguments as if it wasn't decorated.

An exception is raised if an ``NDData`` property is set but the function does
not accept it - for example, if ``wcs`` is set, but the function cannot support
WCS objects, an error would be raised. On the other hand, if an argument in the
function does not exist in the ``NDData`` object or is not set, it is simply
left to its default value.

If the function call succeeds, then the decorator returns the values from the
function unmodified by default. However, in some cases we may want to return
separate ``data``, ``wcs``, etc. if these were passed in separately, and a new
:class:`~astropy.nddata.NDData` instance otherwise. To do this, you can specify
``repack=True`` in the decorator and provide a list of the names of the output
arguments from the function::

    @support_nddata(repack=True, returns=['data', 'wcs'])
    def test(data, wcs=None, unit=None, n_iterations=3):
        ...

With this, the function will return separate values if ``test`` is called with
separate arguments, and an object with the same class type as the input if the
input is an :class:`~astropy.nddata.NDData` or subclass instance.

Finally, the decorator can be made to restrict input to specific ``NDData``
sub-classes (and sub-classes of those) using the ``accepts`` option::

    @support_nddata(accepts=CCDImage)
    def test(data, wcs=None, unit=None, n_iterations=3):
        ...

