.. _nddata_subclassing:

Subclassing
===========

`~astropy.nddata.NDData`
------------------------

This class serves as the base for subclasses that use a `numpy.ndarray` (or
something that presents a numpy-like interface) as the ``data`` attribute.

.. note::
  Each attribute is saved as attribute with one leading underscore. For example
  the ``data`` is saved as ``_data`` and the ``mask`` as ``_mask``, and so on.

Adding another property
^^^^^^^^^^^^^^^^^^^^^^^

    >>> from astropy.nddata import NDData

    >>> class NDDataWithFlags(NDData):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super(NDDataWithFlags, self).__init__(*args, **kwargs)
    ...
    ...     @property
    ...     def flags(self):
    ...         return self._flags
    ...
    ...     @flags.setter
    ...     def flags(self, value):
    ...         self._flags = value

    >>> ndd = NDDataWithFlags([1,2,3])
    >>> ndd.flags is None
    True

    >>> ndd = NDDataWithFlags([1,2,3], flags=[0, 0.2, 0.3])
    >>> ndd.flags
    [0, 0.2, 0.3]

.. note::
  To simplify subclassing each setter (except for ``data``) is called during
  ``__init__`` so putting restrictions on any attribute can be done inside
  the setter and will also apply duing instance creation.

Customize the setter for a property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    >>> import numpy as np

    >>> class NDDataMaskBoolNumpy(NDData):
    ...
    ...     @NDData.mask.setter
    ...     def mask(self, value):
    ...         # Convert mask to boolean numpy array.
    ...         self._mask = np.array(value, dtype=np.bool_)

    >>> ndd = NDDataMaskBoolNumpy([1,2,3], mask=True)
    >>> ndd.mask
    array(True, dtype=bool)

Extend the setter for a property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``unit``, ``meta`` and ``uncertainty`` implement some additional logic in their
setter so subclasses might define a call to the superclass and let the
super property set the attribute afterwards::

    >>> import numpy as np

    >>> class NDDataUncertaintyShapeChecker(NDData):
    ...
    ...     @NDData.uncertainty.setter
    ...     def uncertainty(self, value):
    ...         value = np.asarray(value)
    ...         if value.shape != self.data.shape:
    ...             raise ValueError('uncertainty must have the same shape as the data.')
    ...         # Call the setter of the super class in case it might contain some
    ...         # important logic (only True for meta, unit and uncertainty)
    ...         super(NDDataUncertaintyShapeChecker, self.__class__).uncertainty.fset(self, value)

    >>> ndd = NDDataUncertaintyShapeChecker([1,2,3], uncertainty=[2,3,4])
    INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
    >>> ndd.uncertainty
    UnknownUncertainty([2, 3, 4])

Having a setter for the data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    >>> class NDDataWithDataSetter(NDData):
    ...
    ...     @NDData.data.setter
    ...     def data(self, value):
    ...         # Convert mask to numpy array
    ...         self._data = np.asarray(value)

    >>> ndd = NDDataWithDataSetter([1,2,3])
    >>> ndd.data = [3,2,1]
    >>> ndd.data
    array([3, 2, 1])

`~astropy.nddata.NDDataRef`
---------------------------

`~astropy.nddata.NDDataRef` itself inherits from `~astropy.nddata.NDData` so
any of the possibilities there also apply to NDDataRef. But NDDataRef also
inherits from the Mixins:

- `~astropy.nddata.NDSlicingMixin`
- `~astropy.nddata.NDArithmeticMixin`
- `~astropy.nddata.NDIOMixin`

which allow additional operations.

Slicing an existing property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose you have a class expecting a 2 dimensional ``data`` but the mask is
only 1D. This would lead to problems if one were to slice in two dimensions.

    >>> from astropy.nddata import NDDataRef
    >>> import numpy as np

    >>> class NDDataMask1D(NDDataRef):
    ...     def _slice_mask(self, item):
    ...         # Multidimensional slices are represented by tuples:
    ...         if isinstance(item, tuple):
    ...             # only use the first dimension of the slice
    ...             return self.mask[item[0]]
    ...         # Let the superclass deal with the other cases
    ...         return super(NDDataMask1D, self)._slice_mask(item)

    >>> ndd = NDDataMask1D(np.ones((3,3)), mask=np.ones(3, dtype=bool))
    >>> nddsliced = ndd[1:3,1:3]
    >>> nddsliced.mask
    array([ True,  True], dtype=bool)

.. note::
  The methods doing the slicing of the attributes are prefixed by a
  ``_slice_*`` where ``*`` can be ``mask``, ``uncertainty`` or ``wcs``. So
  simply overriding them is the easiest way to customize how the are sliced.

.. note::
  If slicing should affect the ``unit`` or ``meta`` see the next example.


Slicing an additional property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building on the added property ``flags`` we want them to be sliceable:

    >>> class NDDataWithFlags(NDDataRef):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super(NDDataWithFlags, self).__init__(*args, **kwargs)
    ...
    ...     @property
    ...     def flags(self):
    ...         return self._flags
    ...
    ...     @flags.setter
    ...     def flags(self, value):
    ...         self._flags = value
    ...
    ...     def _slice(self, item):
    ...         # slice all normal attributes
    ...         kwargs = super(NDDataWithFlags, self)._slice(item)
    ...         # The arguments for creating a new instance are saved in kwargs
    ...         # so we need to add another keyword "flags" and add the sliced flags
    ...         kwargs['flags'] = self.flags[item]
    ...         return kwargs # these must be returned

    >>> ndd = NDDataWithFlags([1,2,3], flags=[0, 0.2, 0.3])
    >>> ndd2 = ndd[1:3]
    >>> ndd2.flags
    [0.2, 0.3]

If you wanted to keep just the original ``flags`` instead of the sliced ones
you could use ``kwargs['flags'] = self.flags`` and omit the ``[item]``.


Arithmetic on an existing property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Customizing how an existing property is handled during arithmetics is possible
with some arguments to the function calls like
:meth:`~astropy.nddata.NDArithmeticMixin.add` but if it's possible to hardcode
behaviour too. The actual operation on the attribute (except for ``unit``) is
done in a method ``_arithmetic_*`` where ``*`` is the name of the property.

For example to customize how the ``meta`` will be affected during arithmetics::

    >>> from copy import deepcopy
    >>> class NDDataWithMetaArithmetics(NDDataRef):
    ...
    ...     def _arithmetic_meta(self, operation, operand, handle_mask, **kwds):
    ...         # the function must take the arguments:
    ...         # operation (numpy-ufunc like np.add, np.subtract, ...)
    ...         # operand (the other NDData-like object, already wrapped as NDData)
    ...         # handle_mask (see description for "add")
    ...
    ...         # The meta is dict like but we want the keywords exposure to change
    ...         # Anticipate that one or both might have no meta and take the first one that has
    ...         result_meta = deepcopy(self.meta) if self.meta else deepcopy(operand.meta)
    ...         # Do the operation on the keyword if the keyword exists
    ...         if result_meta and 'exposure' in result_meta:
    ...             result_meta['exposure'] = operation(result_meta['exposure'], operand.data)
    ...         return result_meta # return it

To trigger this method the ``handle_meta`` argument must not be ``None``,
``"ff"`` or ``"first_found"``::

    >>> ndd = NDDataWithMetaArithmetics([1,2,3], meta={'exposure': 10})
    >>> ndd2 = ndd.add(10, handle_meta='')
    >>> ndd2.meta
    {'exposure': 20}

    >>> ndd3 = ndd.multiply(0.5, handle_meta='')
    >>> ndd3.meta
    {'exposure': 5.0}

.. warning::
  To use these internal `_arithmetic_*` methods there are some restrictions on
  the attributes when calling the operation:

  - ``mask``: ``handle_mask`` must not be ``None``, ``"ff"`` or ``"first_found"``.
  - ``wcs``: ``compare_wcs`` argument with the same restrictions as mask.
  - ``meta``: ``handle_meta`` argument with the same restrictions as mask.
  - ``uncertainty``: ``propagate_uncertainties`` must be ``None`` or evaluate
    to ``False``. ``arithmetic_uncertainty`` must also accepts different
    arguments: ``operation, operand, result, correlation, **kwargs``


Arithmetic with an additional property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This also requires overriding the ``_arithmetic`` method. Suppose we have a
``flags`` attribute again::

    >>> from copy import deepcopy
    >>> import numpy as np

    >>> class NDDataWithFlags(NDDataRef):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super(NDDataWithFlags, self).__init__(*args, **kwargs)
    ...
    ...     @property
    ...     def flags(self):
    ...         return self._flags
    ...
    ...     @flags.setter
    ...     def flags(self, value):
    ...         self._flags = value
    ...
    ...     def _arithmetic(self, operation, operand, *args, **kwargs):
    ...         # take all args and kwargs to allow arithmetic on the other properties
    ...         # to work like before.
    ...
    ...         # do the arithmetics on the flags (pop the relevant kwargs, if any!!!)
    ...         if self.flags is not None and operand.flags is not None:
    ...             result_flags = np.logical_or(self.flags, operand.flags)
    ...             # np.logical_or is just a suggestion you can do what you want
    ...         else:
    ...             if self.flags is not None:
    ...                 result_flags = deepcopy(self.flags)
    ...             else:
    ...                 result_flags = deepcopy(operand.flags)
    ...
    ...         # Let the superclass do all the other attributes note that
    ...         # this returns the result and a dictionary containing other attributes
    ...         result, kwargs = super(NDDataWithFlags, self)._arithmetic(operation, operand, *args, **kwargs)
    ...         # The arguments for creating a new instance are saved in kwargs
    ...         # so we need to add another keyword "flags" and add the processed flags
    ...         kwargs['flags'] = result_flags
    ...         return result, kwargs # these must be returned

    >>> ndd1 = NDDataWithFlags([1,2,3], flags=np.array([1,0,1], dtype=bool))
    >>> ndd2 = NDDataWithFlags([1,2,3], flags=np.array([0,0,1], dtype=bool))
    >>> ndd3 = ndd1.add(ndd2)
    >>> ndd3.flags
    array([ True, False,  True], dtype=bool)

Another arithmetic operation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
    Add example!)


`~astropy.nddata.NDDataBase`
----------------------------

The class `~astropy.nddata.NDDataBase` is a metaclass -- when subclassing it,
all properties of `~astropy.nddata.NDDataBase` *must* be overriden in the
subclass.

Subclassing from `~astropy.nddata.NDDataBase` gives you complete flexibility
in how you implement data storage and the other properties. If your data is
stored in a numpy array (or something that behaves like a numpy array), it may
be more straightforward to subclass `~astropy.nddata.NDData` instead of
`~astropy.nddata.NDDataBase`.

Implementing the NDDataBase interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example to create a readonly container.

    >>> from astropy.nddata import NDDataBase

    >>> class NDDataReadOnlyNoRestrictions(NDDataBase):
    ...     def __init__(self, data, unit, mask, uncertainty, meta, wcs):
    ...         self._data = data
    ...         self._unit = unit
    ...         self._mask = mask
    ...         self._uncertainty = uncertainty
    ...         self._meta = meta
    ...         self._wcs = wcs
    ...
    ...     @property
    ...     def data(self):
    ...         return self._data
    ...
    ...     @property
    ...     def unit(self):
    ...         return self._unit
    ...
    ...     @property
    ...     def mask(self):
    ...         return self._mask
    ...
    ...     @property
    ...     def uncertainty(self):
    ...         return self._uncertainty
    ...
    ...     @property
    ...     def meta(self):
    ...         return self._meta
    ...
    ...     @property
    ...     def wcs(self):
    ...         return self._wcs

    >>> # A meaningless test to show that creating this class is possible:
    >>> NDDataReadOnlyNoRestrictions(1,2,3,4,5,6) is not None
    True

..note::
  Actually defining an ``__init__`` is not necessary and the properties could
  return arbitary values but the properties **must** be defined.

Subclassing `~astropy.nddata.NDUncertainty`
-------------------------------------------

.. warning::
    This section is old and probably outdated. Will be provided soon... maybe :-)

This is done by using classes to represent the uncertainties of a given type.
For example, to set standard deviation uncertainties on the pixel values, you
can do::

    >>> import numpy as np
    >>> from astropy.nddata import NDData, StdDevUncertainty
    >>> array = np.zeros((12, 12, 12))  # a 3-dimensional array with all zeros
    >>> ndd = NDData(array)
    >>> uncertainty = StdDevUncertainty(np.ones((12, 12, 12)) * 0.1)
    >>> ndd_uncertainty = NDData(ndd, uncertainty=uncertainty)

New error classes should sub-class from `~astropy.nddata.NDUncertainty`, and
should provide methods with the following API::

   class MyUncertainty(NDUncertainty):

       def propagate_add(self, other_nddata, result_data):
           ...
           result_uncertainty = MyUncertainty(...)
           return result_uncertainty

       def propagate_subtract(self, other_nddata, result_data):
           ...
           result_uncertainty = MyUncertainty(...)
           return result_uncertainty

       def propagate_multiply(self, other_nddata, result_data):
           ...
           result_uncertainty = MyUncertainty(...)
           return result_uncertainty

       def propagate_divide(self, other_nddata, result_data):
           ...
           result_uncertainty = MyUncertainty(...)
           return result_uncertainty

All error sub-classes inherit an attribute ``self.parent_nddata`` that is
automatically set to the parent `~astropy.nddata.NDData` object that they
are attached to. The arguments passed to the error propagation methods are
``other_nddata``, which is the `~astropy.nddata.NDData` object that is being
combined with ``self.parent_nddata``, and ``result_data``, which is a Numpy
array that contains the data array after the arithmetic operation. All these
methods should return an error instance ``result_uncertainty``, and should not
modify ``parent_nddata`` directly. For subtraction and division, the order of
the operations is ``parent_nddata - other_nddata`` and ``parent_nddata /
other_nddata``.

To make it easier and clearer to code up the error propagation, you can use
variables with more explicit names, e.g::

   class MyUncertainty(NDUncertainty):

       def propogate_add(self, other_nddata, result_data):

           left_uncertainty = self.parent.uncertainty.array
           right_uncertainty = other_nddata.uncertainty.array

           ...

Note that the above example assumes that the errors are stored in an ``array``
attribute, but this does not have to be the case.

For an example of a complete implementation, see `~astropy.nddata.StdDevUncertainty`.
