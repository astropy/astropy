******************************
Creating and Evaluating Models
******************************

The base class of all models is `~astropy.modeling.core.Model`, however
fittable models should subclass `~astropy.modeling.core.ParametricModel`.
Parametric models can be linear or nonlinear in a regression analysis sense.

Model instances are callabale, that is, o evaluate a model it is called like a
function. When possible the transformation is done using multiple
:attr:`param_sets <astropy.modeling.core.Model.param_sets>`.  The number of
parameter sets is stored in an attribute
`~astropy.modeling.core.Model.param_dim`.

Parametric models also store a flat array of all parameter values.  When
fitting, this array is directly modified by a subclass of
`~astropy.modeling.fitting.Fitter`, in turn updating all of the model's
parameter values simultaneously.  When fitting nonlinear models, the values of
the parameters are used as initial guesses by the fitting class.

Models have an `~astropy.modeling.core.Model.n_inputs` attribute, which shows
how many coordinates the model expects as an input. All models expect
coordinates as separate arguments.  For example a 2D model expects x and y to
be passed separately, e.g. as two arrays or two lists. When a model has
multiple parameter sets and x, y are 2D arrays, the model is evaluated with
each of the parameter sets and the same x, y as input. The shape of the output
array is ``(param_dim, x_shape, y_shape)`` where param_dim is the number of
parameter sets and ``x_shape, y_shape`` is the shape of the input array.  In
all other cases the shape of the output array is the same as the shape of the
input arrays.

Models also have an attribute `~astropy.modeling.core.Model.n_outputs`, which
shows the number of output coordinates. The
`~astropy.modeling.core.Model.n_inputs` and
`~astropy.modeling.core.Model.n_outputs` attributes are used to chain
transforms by adding models in :class:`series
<astropy.modeling.core.SerialCompositeModel>` or in :class:`parallel
<astropy.modeling.core.SummedCompositeModel>`. Because composite models can
be nested within other composite models, creating theoretically infinitely
complex models, a mechanism to map input data to models is needed. In this case
the input may be wrapped in a `~astropy.modeling.core.LabeledInput` object-- a
dict-like object whose items are ``{label: data}`` pairs.


Model examples
--------------

The examples here assume this import statement was executed::

    >>> from astropy.modeling import *
    >>> import numpy as np

- Create a 1D Gaussian with 2 parameter sets::

    >>> x = np.arange(1, 10, .1)
    >>> models.Gaussian1D.param_names
    ['amplitude', 'mean', 'stddev']
    >>> g1 = models.Gaussian1D(amplitude=[10, 9], mean=[2, 3],
    ...                        stddev=[0.15, .1])
    >>> g1.param_sets
    array([[ 10.  ,   9.  ],
           [  2.  ,   3.  ],
           [  0.15,   0.1 ]])

  Evaluate the model on one data set::

      >>> y = g1(x)
      >>> print(y.shape)
      (90, 2)

  or two data sets (any other number would be an error)::

      >>> y = g1(np.array([x, x]).T)
      >>> print(y.shape)
      (90, 2)

.. plot::

   import matplotlib.pyplot as plt
   import numpy as np
   from astropy.modeling import models, fitting
   x = np.arange(1, 10, .1)
   g1 = models.Gaussian1D(amplitude=[10, 9], mean=[2,3], stddev=[.15,.1])
   y = g1(x)
   plt.plot(x, y)
   plt.title('Evaluate a Gaussian1D model with 2 parameter sets and 1 set of '
             'input data')
   plt.show()

.. plot::

   import matplotlib.pyplot as plt
   import numpy as np
   from astropy.modeling import models, fitting
   x = np.arange(1, 10, .1)
   g1 = models.Gaussian1D(amplitude=[10, 9], mean=[2,3], stddev=[.15,.1])
   y = g1(np.array([x, x]).T)
   plt.plot(x, y)
   plt.title('Evaluating a Gaussian1D model with 2 parameter sets and 2 sets '
             'of input data')
   plt.show()


- Evaluating polynomial models with multiple parameter sets with one input data
  set creates multiple output data sets::

    >>> len(p1.parameters)  # doctest: +SKIP
    10
    >>> p1.c1 = [0, 1, 2, 3, 4]  # doctest: +SKIP
    >>> p1.param_sets  # doctest: +SKIP
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  2.,  3.,  4.]])
    >>> y = p1(x)  # doctest: +SKIP


.. plot::

   import matplotlib.pyplot as plt
   import numpy as np
   from astropy.modeling import models, fitting
   x = np.arange(1, 10, .1)
   p1 = models.Polynomial1D(1, param_dim=5)
   p1.c1 = [0, 1, 2, 3, 4]
   y = p1(x)
   plt.plot(x, y)
   plt.title("Polynomial1D model with 5 parameter sets")
   plt.show()

- When passed a 2D array, the same polynomial will map parameter sets to array
  columns::

    >>> x = np.ones((10,5))
    >>> y = p1(x)  # doctest: +SKIP
    >>> print(y)  # doctest: +SKIP
    array([[ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.],
           [ 0.,  1.,  2.,  3.,  4.]])
    >>> print(y.shape)  # doctest: +SKIP
    (10,5)

- Create and evaluate a parallel composite model::

    >>> x = np.arange(1,10,.1)
    >>> p1 = models.Polynomial1D(1)
    >>> g1 = models.Gaussian1D(amplitude=10., stddev=2.1, mean=4.2)
    >>> sum_of_models = SummedCompositeModel([g1, p1])
    >>> y = sum_of_models(x)

  This is equivalent to applying the two models in parallel::

      >>> y = x + g1(x) + p1(x)

In more complex cases the input and output may be mapped to transformations::

    >>> x, y = np.mgrid[:5, :5]
    >>> off = models.Shift(-3.2)
    >>> poly2 = models.Polynomial2D(2)
    >>> serial_composite_model = SerialCompositeModel(
    ...     [off, poly2], inmap=[['x'], ['x', 'y']], outmap=[['x'], ['z']])

The above composite transform will apply an inplace shift to x, followed by a
2D polynomial and will save the result in an array, labeled 'z'.  To evaluate
this model use a `~astropy.modeling.core.LabeledInput` object::

    >>> labeled_data = LabeledInput([x, y], ['x', 'y'])
    >>> result = serial_composite_model(labeled_data)

The output is also a `~astropy.modeling.core.LabeledInput` object and the
result is stored in label 'z'::

    >>> print(result)  # doctest: +SKIP
    {'x': array([[-3.2, -3.2, -3.2, -3.2, -3.2],
           [-2.2, -2.2, -2.2, -2.2, -2.2],
           [-1.2, -1.2, -1.2, -1.2, -1.2],
           [-0.2, -0.2, -0.2, -0.2, -0.2],
           [ 0.8,  0.8,  0.8,  0.8,  0.8]]),
     'y': array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]]),
     'z': array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])}
