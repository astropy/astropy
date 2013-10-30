**********
Parameters
**********

Parameters are used in three different contexts within this package:
communicating with fitters, model evaluation and getting values to/from users.

Models maintain a list of parameter names,
`~astropy.modeling.core.Model.param_names`.  Single parameters are instances of
`~astropy.modeling.parameters.Parameter` which provide a proxy for the actual
parameter values.  Simple mathematical operations can be performed with them,
but they also contain additional attributes specific to model parameters, such
as any constraints on their values. The preferred way for users to interact
with models is through individual parameters.

The goal of this package is, when possible, to allow simultaneous model
evaluation and fitting with multiple parameter sets. Because of this, all
models have a `~astropy.modeling.core.Model.param_sets` attribute, an array of
shape ``(len(param_names), param_dim)``, where
`~astropy.modeling.core.Model.param_dim` is the number of parameter sets.
Typically the array is of type float but can become an object array in some
cases. `~astropy.modeling.core.Model.param_sets` is used for model evaluation.

In addition, fittable models maintain an attribute,
`~astropy.modeling.core.ParametricModel.parameters`, which is a flattened 1D
array of parameter values. It serves as the primary storage of the raw values
of fittable models' parameters, and is used directly by fitters as an efficient
means of reading and updating a model's parameters.


Parameter examples
------------------

- Polynomial models are created by default with all coefficients set to 0::

    >>> from astropy.modeling import *
    >>> p1 = models.Polynomial1D(degree=4)
    >>> p1.param_names
    ['c0', 'c1', 'c2', 'c3', 'c4']
    >>> p1.parameters
    array([ 0.,  0.,  0.,  0.,  0.])

- Coefficients can be set using the
  `~astropy.modeling.core.ParametricModel.parameters` attribute::

    >>> p1.parameters = [0, 1, 2, 3, 4]
    >>> p1.parameters
    array([ 0.,  1.,  2.,  3.,  4.])

- It is possible to set the coefficients passing the parameters in a
  dictionary::

    >>> ch2 = models.Chebyshev2D(x_degree=2, y_degree=3, param_dim=2)
    >>> coeff = dict((name, [idx, idx + 10])
    ...              for idx, name in enumerate(ch2.param_names))
    >>> ch2 = models.Chebyshev2D(x_degree=2, y_degree=3, **coeff)
    INFO: Inferred 2 dimensions when creating a Chebyshev2D model. Resetting param_dim to 2 [astropy.modeling.polynomial] 
    >>> ch2.param_sets
    array([[  0.,  10.],
           [  1.,  11.],
           [  2.,  12.],
           [  3.,  13.],
           [  4.,  14.],
           [  5.,  15.],
           [  6.,  16.],
           [  7.,  17.],
           [  8.,  18.],
           [  9.,  19.],
           [ 10.,  20.],
           [ 11.,  21.]])

- or directly, using keyword arguments::

    >>> ch2 = models.Chebyshev2D(x_degree=2, y_degree=3,
    ...                               c0_0=[0, 10], c0_1=[3, 13],
    ...                               c0_2=[6, 16], c0_3=[9, 19],
    ...                               c1_0=[1, 11], c1_1=[4, 14],
    ...                               c1_2=[7, 17], c1_3=[10, 20,],
    ...                               c2_0=[2, 12], c2_1=[5, 15],
    ...                               c2_2=[8, 18], c2_3=[11, 21])
    INFO: Inferred 2 dimensions when creating a Chebyshev2D model. Resetting param_dim to 2 [astropy.modeling.polynomial]

- It is possible to change a single parameter::

    >>> ch2.c0_0
    Parameter('c0_0', value=array([  0.,  10.]))
    >>> ch2.c0_0[0] = -34.2
    >>> ch2.c0_0
    Parameter('c0_0', value=array([-34.2,  10. ]))

- The number of parameter sets is stored in an attribute
  `~astropy.modeling.core.Model.param_dim`::

    >>> ch2.param_dim
    2
