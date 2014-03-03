Subclassing `~astropy.nddata.NDData` and `~astropy.nddata.NDUncertainty`
=======================================================================================

Subclassing `~astropy.nddata.NDUncertainty`
---------------------------------------------------

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
