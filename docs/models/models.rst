.. _models:

******
Models
******

The base class of all models is `~astropy.models.models.Model`, however fittable
models should subclass `~astropy.models.models.ParametricModel`. Parametric 
models can be linear or nonlinear in a regression analysis sense.

All models provide a `__call__` method which performs the transformation in a 
purely mathematical way, i.e. the models are unitless. In addition, when possible the 
transformation is done using multiple parameter sets, `psets`.
The number of parameter sets is stored in an attribute `paramdim`. 

Parametric models also store a flat list of all parameters as an instance of 
`~astropy.models.parameters.Parameters`. When fitting, this list-like object is
modified by a subclass of `~astropy.fitting.fitting.Fitter`. When fitting nonlinear models,
the values of the parameters are used as initial guesses by the fitting class.

Models have dimensions, `ndim` attribute, which show how many coordinates the 
model expects as an input. All models expect coordinates as separate arguments.
For example a 2D model expects x and y to be passed separately, 
e.g. as two arrays or two lists. When a model has multiple parameter sets and x, y are 
2D arrays, the model is evaluated with each of the parameter sets and the same x, y as 
input. The shape of the  output array is (paramdim, xsh, ysh) where paramdim is the number 
of parameter sets and xsh, ysh is the shape of the input array.
In all other cases the shape of the output array is the same as the shape of the 
input arrays. 

Models also have an attribute  `outdim`, which shows the number of output 
coordinates. The `ndim` and `outdim` attributes are used to chain transforms by
adding models (in series or in  parallel), to form an 
instance of `~astropy.models.models._CompositeModel`.  Because composite models can 
be nested within other composite models, creating 
theoretically infinetely complex models, a mechanism to map input data to models 
is needed. In this case the input may be wrapped in a
`~astropy.models.models.LabeledInput` object - a dict like object whose items are {label: data} pairs.

Models Examples
---------------

- Create a 1D Gaussian with 2 parameter sets

>>> x=np.arange(1,10,.1)
>>>models.Gauss1DModel.parnames
['amplitude', 'xcen', 'xsigma']
>>> g1=models.Gauss1DModel(amplitude=[10, 9], xcen=[2,3], fwhm=[.3,.2])
>>> g1.psets
array([[ 10.      ,   9.      ],
       [  2.      ,   3.      ],
       [  0.127398,   0.084932]])

Evaluate the model on one data set

>>> y = g1(x)
>>> print y.shape
(90, 2)

or two data sets (any other number would be an error)

>>> y=g1(np.array([x,x]).T)
>>> print y.shape
(90, 2)

.. figure:: images/gauss1D_eval_psets.png
   :scale: 25 %

- Evaluating polynomial models with multiple parameter sets with one input data set creates multiple output data sets

>>> p1 = models.Poly1DModel(1, paramdim=5)
>>> len(p1.parameters)
10
>>> p1.c1=[0, 1, 2, 3, 4]
>>> p1.psets
array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  1.,  2.,  3.,  4.]])
>>> y = p1(x)

.. image:: images/p1d_5psets.png
   :scale: 25 %

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
>>> p1 = models.Poly1DModel(1)
>>> g1 = models.Gauss1DModel(10., xsigma=2.1, xcen=4.2)
>>> pcomptr = models.PCompositeModel([g1, p1])
>>> y = pcomptr(x)

This is equivalent to applying the two models in parallel:

>>> y = x + (g1(x) - x) + (p1(x) - x)

In more complex cases the input and output may be mapped to transformations:

>>> x, y = np.mgrid[:5, :5]
>>> off = models.ShiftModel(-3.2)
>>> poly2 = models.Poly2DModel(2)
>>> scomptr = models.SCompositeModel([off, poly2], inmap=[['x'], ['x', 'y']], outmap=[['x'], ['z']])

The above composite transform will apply an inplace shift to x, followed by a 2D 
polynomial and will save the result in an array, labeled 'z'.
To evaluate this model use a LabeledInput object

>>> ado = models.LabeledInput([x, y], ['x', 'y'])
>>> result = scomptr(ado)

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


