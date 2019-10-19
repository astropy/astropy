.. _add_units:

Adding support for units in a model (Advanced)
==============================================

Evaluation
----------

To make it so that your models can accept parameters with units and be evaluated
using inputs with units, you need to make sure that the
:meth:`~astropy.modeling.Model.evaluate` method works correctly with
input values and parameters with units. For simple arithmetic, this may work
out of the box since :class:`~astropy.units.Quantity` objects are understood by
a number of Numpy functions.

If users of your models provide input during evaluation that is not compatible
with the parameter units, they may get cryptic errors such as::

    UnitsError : Can only apply 'subtract' function to dimensionless quantities
    when other argument is not a quantity (unless the latter is all
    zero/infinity/nan)

There are several attributes or properties that can be set on models that adjust
the behavior of models with units. These attributes can be changed from the
defaults in the class definition, e.g.::

    class MyModel(Model):
        input_units = {'x': u.deg}
        ...

Note that these are all optional.

.. _models_input_units:

``input_units``
^^^^^^^^^^^^^^^

You can easily add checking of the input units by adding an ``input_units``
property or attribute on your model class. This should return either `None` (to
indicate no constraints) or a dictionary where the keys are the input names
(e.g. ``x`` for many 1D models) and the values are the units expected, which can
be a function of the parameter units::

    @property
    def input_units(self):
        if self.mean.unit is None:
            return None
        else:
            return {'x': self.mean.unit}

If the user then gives values with incorrect input units, a clear error will be
displayed::

    UnitsError: Units of input 'x', (dimensionless), could not be converted to
    required input units of m (length)

Note that the input units don't have to match exactly those returned by
``input_units``, but be convertible to them. In addition, ``input_units`` can
also be specified as an attribute rather than a property in simple cases::

    input_units = {'x': u.deg}

``return_units``
^^^^^^^^^^^^^^^^

Similarly to :ref:`models_input_units`, this should be dictionary that maps the return
values of a model to units. If :meth:`~astropy.modeling.Model.evaluate` was called
with quantities but returns unitless values, the units are added to the output.
If the return values are quantities in different units, they are converted to
``return_units``.

``input_units_strict``
^^^^^^^^^^^^^^^^^^^^^^

If set to `True`, values that are passed in compatible units will be converted
to the exact units specified in ``input_units``.

This attribute can also be a
dictionary that maps input names to a Boolean to enable converting of that input
to the specified unit.

``input_units_equivalencies``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be set to a dictionary that maps the input names to a list of
equivalencies, for example::

    input_units_equivalencies = {'nu': u.spectral()}

``_input_units_allow_dimensionless``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If set to `True`, values that are plain scalars or Numpy arrays can be passed to
evaluate even if ``input_units`` specifies that the input should have units. It
is up to the :meth:`~astropy.modeling.Model.evaluate` to then decide how to
handle these dimensionless values. This attribute can also be a dictionary that
maps input names to a Boolean to enable passing dimensionless values to
:meth:`~astropy.modeling.Model.evaluate` for that input.


Fitting
-------

To allow models with parameters that have units to be fitted to data with units,
you will need to add a method called ``_parameter_units_for_data_units`` to your
model class. This should take two arguments ``input_units`` and
``output_units`` - ``input_units`` will be set to a dictionary with
the units of the independent variables in the data, while ``output_units`` will
be set to a dictionary with the units the dependent variables in the data (for
example, for a simple 1D model, ``input_units`` will have one key, ``x``, and
``output_units`` will have one key, ``y``). This method should then return
a dictionary giving for each parameter the units the parameter should be
converted to so that the model could be used on the data if units were removed
from both the models and the data. The following example shows the
implementation for the 1D Gaussian::

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {'mean': inputs_unit['x'],
                'stddev': inputs_unit['x'],
                'amplitude': outputs_unit['y']}

With this method in place, the model can then be fit to data that has units.
