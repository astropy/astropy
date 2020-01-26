.. include:: links.inc

.. _modeling-major-changes-for-4.0:

Changes to Modeling in v4.0
***************************

In order to make the internal code less complex, improve performance, and
make the behavior of parameters in compound models more intuitive,
many changes have been made internally to the modeling code,
particularly to compound models and parameters. This page
summarizes the important changes. More technical details are
given at the end, but it is generally not necessary to read those
unless you want to understand why some changes were necessary.

- Support for expressions of modeling classes has been removed.
  Expressions of model class instances are still fully supported.
  This was done to streamline the implementation, improve performance,
  and support the new parameter semantics. For example:

  No longer works::

        NewCompoundClass = Gaussian1D + Const1D

  Still works::

        newmodel = Gaussian1D(3.2, 7.1, 2.1) + Const1D(3.)

- Previous to v4.0, parameters were class descriptors, which meant
  that they could not hold values for the models. Instead, the
  values were held inside the models. This resulted in confusion
  when compound models were used since this necessitated that
  the compound models make copies of the values. As a result,
  changing the value in the compound model did not change the
  constituent model's parameter value and vice versa. Now
  parameters are distinct instances for each use and they
  do hold the value of the parameter, so compound models
  now share the same values as the constituent models.

- Previously when model sets were used, the parameter shape
  did not show the corresponding dimension for the number
  of models. Now it does. For example:

  Old::

        In [1]: g = Gaussian1D([1,1], [1,2], [2,4], n_models=2)
        In [2]: g.amplitude
        Out[2]: Parameter('amplitude', value=[1. 1.])
        In [3]: g.amplitude.shape
        Out[3]: ()

  New::

        In [1]: g = Gaussian1D([1,1], [1,2], [2,4], n_models=2)
        In [2]: g.amplitude
        Out[2]: Parameter('amplitude', value=[1. 1.])
        In [3]: g.amplitude.shape
        Out[3]: (2,)

- Previously the values were held in an array within the model
  instance and it was possible to assign values to slices of
  that array. Reassigning the array did update the parameters,
  but assigning slices does not. The new approach is to either
  replace the whole array by assigning to the parameters property
  or assign directly to the parameter value.

- The use of 'imputed' units, i.e., supplying input/output units
  to a compound model without them but where the component models
  support the ``_parameter_units_for_data_units()`` method is much more
  restricted in its applicability, which will only work when the
  compound expression uses the ``+`` or ``-`` operators. Past behavior
  led to sometimes arbitrary assignments of units, and sometimes
  incorrect units to the parameters.

- Slicing is more restrictive now. Previously a model defined as
  such::

        m = m1 * m2 + m3

  permitted this slice::

        m[1:] # results in m2 + m3

  Now, only submodels in the expression tree (think of it as
  the sequence of operations as performed) are permitted as
  slices. This means some slices that make sense do not work
  now. The following code illustrates what is permitted
  and what isn't::

        m = m1 + m2 + m3 + m4
        m[:2] # Results in m1 + m2, works
        m[1:] # Should result in m2 + m3 + m4, does not work
               # since m1 is part of all subexpressions.



- Generally, all public methods have remained unchanged. The exception
  is ``n_submodels`` that used to be a method but now is a property.

- Many of the non-public methods have changed, particularly for
  compound models.

- The ``_CompoundModelMeta`` metaclass no longer exists.


Technical Details
=================

Parameter-related changes
-------------------------

Previously Python descriptors were used to define parameters. The
drawback of descriptors is that all instances of the models that
they are defined in share the same parameter descriptor instance.
This means the parameter cannot hold different values for different
model instances. The way that this was worked around was to have
the model hold the actual parameter value (and any other information
relevant to that particular instance of the model). The actual
values of the parameters were held in an array that the model
held, and that array held all values for all parameters in one
array. The drawback of this approach was that compound models
had to create their own array holding values for all the parameters
of the compound model. As a result, the parameter values in a
compound model were completely disassociated with those in
the constituent model that made up the compound model. Changes
to one were not reflected in the other, often leading to confusion.

The new implementation still uses attributes defined at the
class level for the parameters, but only
as an intermediate solution. In creating an instance of the model
the class instance is used to create a local instance of a
parameter that is not shared between model instances. So now
parameters hold all the information directly, and can be
shared between compound and constituent models.

One consequence is that the interface to fitting engines still
use the model parameter array, but this array is constructed
from the parameters on the fly. One may set a completely new
array; the model is defined to detect such assignments and
back propagate the values to the parameters. For example, when
doing a fit, the fit results are propagated to the parameters
in such a way. But if one assigns a slice to the model parameter
array, there is no simple way to detect that assignment and
make the necessary parameter updates. So this mode no longer
works unless an explicit method call is made to back propagate
the values. It is not recommended to use this kind of interface
now.

Likewise, parameter-specific attributes also are kept by the
parameter, such as the constraints on the parameter minimum or
maximum, and whether it is fixed as far as fitting is concerned.

Since the parameter no longer is required to link to a specific
model, it holds any dimension corresponding to the model_set
size; when it previously did not.

Compound Model Implementation Changes
-------------------------------------

Compound models previously were implemented using a metaclass
for the compound model while also inheriting from the Model
class (which itself has a metaclass). The primary reason for
this approach was to support expressions of Model classes.
However, this leads to a confusing implementation,
some decrease in performance, and some odd results when
expressions of model classes also include model instances.

The new implementation does away with the metaclass for
compound models and correspondingly no longer supports
expressions of model classes, but only expressions of
model instances. Previously the expression tree was a private
attribute. Now the compound class is itself an expression
tree.

Many of the private methods of Compound Models have changed.
