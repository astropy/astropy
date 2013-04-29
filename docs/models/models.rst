******************************
Creating and Evaluating Models
******************************

The base class of all models is `~astropy.models.core.Model`, however fittable
models should subclass `~astropy.models.core.ParametricModel`. Parametric 
models can be linear or nonlinear in a regression analysis sense.

To evaluate a model, it is called like a function. When possible the 
transformation is done using multiple parameter sets,
`~astropy.models.core.Model.param_sets`.
The number of parameter sets is stored in an attribute
`~astropy.models.core.Model.param_dim`. 

Parametric models also store a flat list of all parameters as an instance of 
`~astropy.models.parameters.Parameters`. When fitting, this list-like object is
modified by a subclass of `~astropy.models.fitting.Fitter`. When fitting nonlinear models,
the values of the parameters are used as initial guesses by the fitting class.

Models have an `~astropy.models.core.Model.ndim` attribute, which shows
how many coordinates the 
model expects as an input. All models expect coordinates as separate arguments.
For example a 2D model expects x and y to be passed separately, 
e.g. as two arrays or two lists. When a model has multiple parameter sets and x, y are 
2D arrays, the model is evaluated with each of the parameter sets and the same x, y as 
input. The shape of the output array is (param_dim, x_shape, y_shape) where param_dim is the number 
of parameter sets and x_shape, y_shape is the shape of the input array.
In all other cases the shape of the output array is the same as the shape of the 
input arrays. 

Models also have an attribute `~astropy.models.core.Model.outdim`, which shows
the number of output coordinates. The `~astropy.models.core.Model.ndim` and
`~astropy.models.core.Model.outdim` attributes are used to chain transforms by
adding models in series, `~astropy.models.core.SCompositeModel`, or in parallel,
`~astropy.models.core.PCompositeModel`. Because composite models can 
be nested within other composite models, creating 
theoretically infinetely complex models, a mechanism to map input data to models 
is needed. In this case the input may be wrapped in a
`~astropy.models.core.LabeledInput` object - a dict like object whose items are {label: data} pairs.

Models Examples
---------------

The examples here assume this import statement was executed:

>>> from astropy.models import *

- Create a 1D Gaussian with 2 parameter sets

>>> x = np.arange(1, 10, .1)
>>> builtin_models.Gaussian1DModel.param_names
['amplitude', 'mean', 'xsigma']
>>> g1 = builtin_models.Gaussian1DModel(amplitude=[10, 9], mean=[2,3], fwhm=[.3,.2])
>>> g1.param_sets
array([[ 10.      ,   9.      ],
       [  2.      ,   3.      ],
       [  0.127398,   0.084932]])

Evaluate the model on one data set

>>> y = g1(x)
>>> print y.shape
(90, 2)

or two data sets (any other number would be an error)

>>> y = g1(np.array([x, x]).T)
>>> print y.shape
(90, 2)

.. plot::

  import matplotlib.pyplot as plt
  import numpy as np
  from astropy.models import builtin_models, fitting
  x = np.arange(1, 10, .1)
  g1 = builtin_models.Gaussian1DModel(amplitude=[10, 9], mean=[2,3], fwhm=[.3,.2])
  y = g1(x)
  plt.plot(x, y)
  plt.title('Evaluate a Gaussian1DModel with 2 parameter sets and 1 set of input data')
  plt.show()
  
.. plot::

  import matplotlib.pyplot as plt
  import numpy as np
  from astropy.models import builtin_models, fitting
  x = np.arange(1, 10, .1)
  g1 = builtin_models.Gaussian1DModel(amplitude=[10, 9], mean=[2,3], fwhm=[.3,.2])
  y = g1(np.array([x, x]).T)
  plt.plot(x, y)
  plt.title('Evaluating a Gaussian1DModel with 2 parameter sets and 2 sets of input data')
  plt.show()
  
  
- Evaluating polynomial models with multiple parameter sets with one input data set creates multiple output data sets

>>> p1 = builtin_models.Poly1DModel(1, param_dim=5)
>>> len(p1.parameters)
10
>>> p1.c1 = [0, 1, 2, 3, 4]
>>> p1.param_sets
array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  2.,  3.,  4.]])
>>> y = p1(x)


.. plot::

  import matplotlib.pyplot as plt
  import numpy as np
  from astropy.models import builtin_models, fitting
  x = np.arange(1, 10, .1)
  p1 = builtin_models.Poly1DModel(1, param_dim=5)
  p1.c1 = [0, 1, 2, 3, 4]
  y = p1(x)
  plt.plot(x, y)
  plt.title("Poly1DModel with 5 parameter sets")
  plt.show()
  
- When passed a 2D array, the same polynomial will map parameter sets to array columns

>>> x = np.ones((10,5))
>>> y = p1(x)
>>> print y
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
>>> print y.shape
(10,5)

- Create and evaluate a parallel composite model

>>> x = np.arange(1,10,.1)
>>> p1 = builtin_models.Poly1DModel(1)
>>> g1 = builtin_models.Gaussian1DModel(10., xsigma=2.1, mean=4.2)
>>> parallel_composite_model = PCompositeModel([g1, p1])
>>> y = parallel_composite_model(x)

This is equivalent to applying the two models in parallel:

>>> y = x + (g1(x) - x) + (p1(x) - x)

In more complex cases the input and output may be mapped to transformations:

>>> x, y = np.mgrid[:5, :5]
>>> off = builtin_models.ShiftModel(-3.2)
>>> poly2 = builtin_models.Poly2DModel(2)
>>> serial_composite_model = SCompositeModel([off, poly2], inmap=[['x'], ['x', 'y']], outmap=[['x'], ['z']])

The above composite transform will apply an inplace shift to x, followed by a 2D 
polynomial and will save the result in an array, labeled 'z'.
To evaluate this model use a LabeledInput object

>>> labeled_data = LabeledInput([x, y], ['x', 'y'])
>>> result = serial_composite_model(labeled_data)

The output is also a LabeledInput object and the result is stored in label 'z'.

>>> print result
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


