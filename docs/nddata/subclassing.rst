Subclassing `~astropy.nddata.nddata.NDData` and `~astropy.nddata.nderror.NDError`
=================================================================================

Subclassing `~astropy.nddata.nderror.NDError`
---------------------------------------------

New error classes should sub-class from `~astropy.nddata.nderror.NDError`, and
should provide methods with the following API::

   class MyError(NDError):

       def propagate_add(self, other_nddata, result_data):
           ...
           result_error = MyError(...)
           return result_error

       def propagate_subtract(self, other_nddata, result_data):
           ...
           result_error = MyError(...)
           return result_error

       def propagate_multiply(self, other_nddata, result_data):
           ...
           result_error = MyError(...)
           return result_error

       def propagate_divide(self, other_nddata, result_data):
           ...
           result_error = MyError(...)
           return result_error

The error class will have an attribute ``self.parent_nddata`` that is set to
the parent `~astropy.nddata.nderror.NDError` object. The arguments passed to
the arithmetic routines are ``other_nddata``, which is the dataset that is
being combined with ``self.parent_nddata``. All these methods should return an
error instance ``result_error``, and should not modify ``parent_nddata``
directly. For subtraction and division, the order of the operations is
``parent_nddata - other_nddata`` and ``parent_nddata / other_nddata``.

To make it easier and clearer to code up the error propagation, you can use
variables with more explicit names, e.g::

   class MyError(NDError):

       def propogate_add(self, other_nddata, result_data):

           left_error = self.parent.error.array
           right_error = other_nddata.error.array

           ...
           
Note that the above example assumes that the errors are stored in an ``array``
attribute, but this does not have to be the case.