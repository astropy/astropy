.. _nddata_subclassing:

Subclassing
***********

`~astropy.nddata.NDData`
========================

This class serves as the base for subclasses that use a `numpy.ndarray` (or
something that presents a ``numpy``-like interface) as the ``data`` attribute.

.. note::
  Each attribute is saved as an attribute with one leading underscore. For
  example, the ``data`` is saved as ``_data`` and the ``mask`` as ``_mask``,
  and so on.

Adding Another Property
-----------------------

    >>> from astropy.nddata import NDData

    >>> class NDDataWithFlags(NDData):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super().__init__(*args, **kwargs)
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
  To simplify subclassing, each setter (except for ``data``) is called during
  ``__init__`` so putting restrictions on any attribute can be done inside
  the setter and will also apply during instance creation.

Customize the Setter for a Property
-----------------------------------

    >>> import numpy as np

    >>> class NDDataMaskBoolNumpy(NDData):
    ...
    ...     @NDData.mask.setter
    ...     def mask(self, value):
    ...         # Convert mask to boolean numpy array.
    ...         self._mask = np.array(value, dtype=np.bool_)

    >>> ndd = NDDataMaskBoolNumpy([1,2,3])
    >>> ndd.mask = [True, False, True]
    >>> ndd.mask
    array([ True, False,  True]...)

Extend the Setter for a Property
--------------------------------

``unit``, ``meta``, and ``uncertainty`` implement some additional logic in their
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
    ...         # Unlike "super(cls_name, cls_name).uncertainty.fset" or
    ...         # or "NDData.uncertainty.fset" this will respect Pythons method
    ...         # resolution order.

    >>> ndd = NDDataUncertaintyShapeChecker([1,2,3], uncertainty=[2,3,4])
    INFO: uncertainty should have attribute uncertainty_type. [astropy.nddata.nddata]
    >>> ndd.uncertainty
    UnknownUncertainty([2, 3, 4])

Having a Setter for the Data
----------------------------

    >>> class NDDataWithDataSetter(NDData):
    ...
    ...     @NDData.data.setter
    ...     def data(self, value):
    ...         self._data = np.asarray(value)

    >>> ndd = NDDataWithDataSetter([1,2,3])
    >>> ndd.data = [3,2,1]
    >>> ndd.data
    array([3, 2, 1])

.. _NDDataRef:

`~astropy.nddata.NDDataRef`
===========================

`~astropy.nddata.NDDataRef` itself inherits from `~astropy.nddata.NDData` so
any of the possibilities there also apply to NDDataRef. But NDDataRef also
inherits from the Mixins:

- `~astropy.nddata.NDSlicingMixin`
- `~astropy.nddata.NDArithmeticMixin`
- `~astropy.nddata.NDIOMixin`

Which allow additional operations.

Add Another Arithmetic Operation
--------------------------------

Adding another operation is possible provided the ``data`` and ``unit`` allow
it within the framework of `~astropy.units.Quantity`.

Examples
^^^^^^^^

..
  EXAMPLE START
  Adding Operations When Working with NDDataRef

To add a power function::

    >>> from astropy.nddata import NDDataRef
    >>> import numpy as np
    >>> from astropy.utils import sharedmethod

    >>> class NDDataPower(NDDataRef):
    ...     @sharedmethod # sharedmethod to allow it also as classmethod
    ...     def pow(self, operand, operand2=None, **kwargs):
    ...         # the uncertainty doesn't allow propagation so set it to None
    ...         kwargs['propagate_uncertainties'] = None
    ...         # Call the _prepare_then_do_arithmetic function with the
    ...         # numpy.power ufunc.
    ...         return self._prepare_then_do_arithmetic(np.power, operand,
    ...                                                 operand2, **kwargs)

This can be used like the other arithmetic methods similar to
:meth:`~astropy.nddata.NDArithmeticMixin.add`. So it works when calling it
on the class or the instance::

    >>> ndd = NDDataPower([1,2,3])

    >>> # using it on the instance with one operand
    >>> ndd.pow(3)
    NDDataPower([ 1,  8, 27])

    >>> # using it on the instance with two operands
    >>> ndd.pow([1,2,3], [3,4,5])
    NDDataPower([  1,  16, 243])

    >>> # or using it as classmethod
    >>> NDDataPower.pow(6, [1,2,3])
    NDDataPower([  6,  36, 216])

To allow propagation also with ``uncertainty`` see subclassing
`~astropy.nddata.NDUncertainty`.

..
  EXAMPLE END

The ``_prepare_then_do_arithmetic`` implements the relevant checks if it was
called on the class or the instance, and if one or two operands were given,
converts the operands, if necessary, to the appropriate classes. Overriding
``_prepare_then_do_arithmetic`` in subclasses should be avoided if
possible.

Arithmetic on an Existing Property
----------------------------------

Customizing how an existing property is handled during arithmetic is possible
with some arguments to the function calls such as
:meth:`~astropy.nddata.NDArithmeticMixin.add`, but it is possible to hardcode
behavior too. The actual operation on the attribute (except for ``unit``) is
done in a method ``_arithmetic_*`` where ``*`` is the name of the property.

Examples
^^^^^^^^

..
  EXAMPLE START
  Customizing Existing Properties During Arithmetic in NDData

To customize how the ``meta`` will be affected during arithmetics::

    >>> from astropy.nddata import NDDataRef

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

To trigger this method, the ``handle_meta`` argument to arithmetic methods can
be anything except ``None`` or ``"first_found"``::

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

  - ``mask``: ``handle_mask`` must not be ``None``, ``"ff"``, or
    ``"first_found"``.
  - ``wcs``: ``compare_wcs`` argument with the same restrictions as mask.
  - ``meta``: ``handle_meta`` argument with the same restrictions as mask.
  - ``uncertainty``: ``propagate_uncertainties`` must be ``None`` or evaluate
    to ``False``. ``arithmetic_uncertainty`` must also accept different
    arguments: ``operation``, ``operand``, ``result``, ``correlation``,
    ``**kwargs``.

..
  EXAMPLE END

Changing the Default Argument for Arithmetic Operations
-------------------------------------------------------

If the goal is to change the default value of an existing parameter for
arithmetic methods, such as when explicitly specifying the parameter each
time you call an arithmetic operation is too much effort, you can change the
default value of existing parameters by changing it in the method signature of
``_arithmetic``.

Example
^^^^^^^

..
  EXAMPLE START
  Changing the Default Argument for Arithmetic Operations in NDData

To change the default value of an existing parameter for arithmetic methods::

    >>> from astropy.nddata import NDDataRef
    >>> import numpy as np

    >>> class NDDDiffAritDefaults(NDDataRef):
    ...     def _arithmetic(self, *args, **kwargs):
    ...         # Changing the default of handle_mask to None
    ...         if 'handle_mask' not in kwargs:
    ...             kwargs['handle_mask'] = None
    ...         # Call the original with the updated kwargs
    ...         return super()._arithmetic(*args, **kwargs)

    >>> ndd1 = NDDDiffAritDefaults(1, mask=False)
    >>> ndd2 = NDDDiffAritDefaults(1, mask=True)
    >>> ndd1.add(ndd2).mask is None  # it will be None
    True

    >>> # But giving other values is still possible:
    >>> ndd1.add(ndd2, handle_mask=np.logical_or).mask
    True

    >>> ndd1.add(ndd2, handle_mask="ff").mask
    False

The parameter controlling how properties are handled are all keyword-only
so using the ``*args``, ``**kwargs`` approach allows you to only alter one
default without needing to care about the positional order of arguments.

..
  EXAMPLE END

Arithmetic with an Additional Property
--------------------------------------

This also requires overriding the ``_arithmetic`` method. Suppose we have a
``flags`` attribute again::

    >>> from copy import deepcopy
    >>> import numpy as np

    >>> class NDDataWithFlags(NDDataRef):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super().__init__(*args, **kwargs)
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
    ...         result, kwargs = super()._arithmetic(operation, operand, *args, **kwargs)
    ...         # The arguments for creating a new instance are saved in kwargs
    ...         # so we need to add another keyword "flags" and add the processed flags
    ...         kwargs['flags'] = result_flags
    ...         return result, kwargs # these must be returned

    >>> ndd1 = NDDataWithFlags([1,2,3], flags=np.array([1,0,1], dtype=bool))
    >>> ndd2 = NDDataWithFlags([1,2,3], flags=np.array([0,0,1], dtype=bool))
    >>> ndd3 = ndd1.add(ndd2)
    >>> ndd3.flags
    array([ True, False,  True]...)

Slicing an Existing Property
----------------------------

Suppose you have a class expecting a 2D ``data`` but the mask is
only 1D. This would lead to problems if you were to slice in two dimensions.

    >>> from astropy.nddata import NDDataRef
    >>> import numpy as np

    >>> class NDDataMask1D(NDDataRef):
    ...     def _slice_mask(self, item):
    ...         # Multidimensional slices are represented by tuples:
    ...         if isinstance(item, tuple):
    ...             # only use the first dimension of the slice
    ...             return self.mask[item[0]]
    ...         # Let the superclass deal with the other cases
    ...         return super()._slice_mask(item)

    >>> ndd = NDDataMask1D(np.ones((3,3)), mask=np.ones(3, dtype=bool))
    >>> nddsliced = ndd[1:3,1:3]
    >>> nddsliced.mask
    array([ True,  True]...)

.. note::
  The methods slicing the attributes are prefixed by a ``_slice_*`` where ``*``
  can be ``mask``, ``uncertainty``, or ``wcs``. So overriding them is the
  most convenient way to customize how the attributes are sliced.

.. note::
  If slicing should affect the ``unit`` or ``meta`` see the next example.


Slicing an Additional Property
------------------------------

Building on the added property ``flags``, we want them to be sliceable:

    >>> class NDDataWithFlags(NDDataRef):
    ...     def __init__(self, *args, **kwargs):
    ...         # Remove flags attribute if given and pass it to the setter.
    ...         self.flags = kwargs.pop('flags') if 'flags' in kwargs else None
    ...         super().__init__(*args, **kwargs)
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
    ...         kwargs = super()._slice(item)
    ...         # The arguments for creating a new instance are saved in kwargs
    ...         # so we need to add another keyword "flags" and add the sliced flags
    ...         kwargs['flags'] = self.flags[item]
    ...         return kwargs # these must be returned

    >>> ndd = NDDataWithFlags([1,2,3], flags=[0, 0.2, 0.3])
    >>> ndd2 = ndd[1:3]
    >>> ndd2.flags
    [0.2, 0.3]

If you wanted to keep just the original ``flags`` instead of the sliced ones,
you could use ``kwargs['flags'] = self.flags`` and omit the ``[item]``.

`~astropy.nddata.NDDataBase`
============================

The class `~astropy.nddata.NDDataBase` is a metaclass — when subclassing it,
all properties of `~astropy.nddata.NDDataBase` *must* be overridden in the
subclass.

Subclassing from `~astropy.nddata.NDDataBase` gives you complete flexibility
in how you implement data storage and the other properties. If your data is
stored in a ``numpy`` array (or something that behaves like a ``numpy`` array),
it may be more convenient to subclass `~astropy.nddata.NDData` instead of
`~astropy.nddata.NDDataBase`.

Example
-------

..
  EXAMPLE START
  Implementing the NDDataBase Interface

To implement the NDDataBase interface by creating a read-only container::

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

.. note::
  Actually defining an ``__init__`` is not necessary and the properties could
  return arbitrary values but the properties **must** be defined.

..
  EXAMPLE END

Subclassing `~astropy.nddata.NDUncertainty`
===========================================

.. warning::
    The internal interface of NDUncertainty and subclasses is experimental and
    might change in future versions.

Subclasses deriving from `~astropy.nddata.NDUncertainty` need in order to
implement:

- Property ``uncertainty_type`` should return a string describing the
  uncertainty, for example, ``"ivar"`` for inverse variance.
- Methods for propagation: `_propagate_*` where ``*`` is the name of the
  universal function (ufunc) that is used on the ``NDData`` parent.

Creating an Uncertainty without Propagation
-------------------------------------------

`~astropy.nddata.UnknownUncertainty` is a minimal working implementation
without error propagation. We can create an uncertainty by storing
systematic uncertainties::

    >>> from astropy.nddata import NDUncertainty

    >>> class SystematicUncertainty(NDUncertainty):
    ...     @property
    ...     def uncertainty_type(self):
    ...         return 'systematic'
    ...
    ...     def _data_unit_to_uncertainty_unit(self, value):
    ...         return None
    ...
    ...     def _propagate_add(self, other_uncert, *args, **kwargs):
    ...         return None
    ...
    ...     def _propagate_subtract(self, other_uncert, *args, **kwargs):
    ...         return None
    ...
    ...     def _propagate_multiply(self, other_uncert, *args, **kwargs):
    ...         return None
    ...
    ...     def _propagate_divide(self, other_uncert, *args, **kwargs):
    ...         return None

    >>> SystematicUncertainty([10])
    SystematicUncertainty([10])
