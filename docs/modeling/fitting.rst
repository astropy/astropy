**********************
Fitting Models to Data
**********************

This module provides wrappers, called Fitters, around some Numpy and Scipy
fitting functions. All Fitters can be called as functions. They take an
instance of `~astropy.modeling.FittableModel` as input and modify its
``parameters`` attribute. The idea is to make this extensible and allow
users to easily add other fitters.

Linear fitting is done using Numpy's `numpy.linalg.lstsq` function.  There are
currently two non-linear fitters which use `scipy.optimize.leastsq` and
`scipy.optimize.fmin_slsqp`.

The rules for passing input to fitters are:

* Non-linear fitters currently work only with single models (not model sets).

* The linear fitter can fit a single input to multiple model sets creating
  multiple fitted models.  This may require specifying the ``model_set_axis``
  argument just as used when evaluating models; this may be required for the
  fitter to know how to broadcast the input data.

* The `~astropy.modeling.fitting.LinearLSQFitter` currently works only with
  simple (not compound) models.

* The current fitters work only with models that have a single output
  (including bivariate functions such as
  `~astropy.modeling.polynomial.Chebyshev2D` but not compound models that map
  ``x, y -> x', y'``).

Plugin Fitters
==============


Fitters defined outside of astropy's core can be inserted into the
`astropy.modeling.fitting` namespace through the use of entry points.
Entry points are references to importable objects. A tutorial on
defining entry points can be found in `setuptools' documentation
<http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`_.
Plugin fitters are required to extend from the `~astropy.modeling.fitting.Fitter`
base class. For the fitter to be discovered and inserted into
`astropy.modeling.fitting` the entry points must be inserted into
the `astropy.modeling` entry point group

.. doctest-skip::

    setup(
          # ...
          entry_points = {'astropy.modeling': 'PluginFitterName = fitter_module:PlugFitterClass'}
    )

This would allow users to import the ``PlugFitterName`` through `astropy.modeling.fitting` by

.. doctest-skip::

    from astropy.modeling.fitting import PlugFitterName

One project which uses this functionality is `Saba <https://saba.readthedocs.io/>`_,
which insert its `SherpaFitter <http://saba.readthedocs.io/en/stable/api.html#saba.SherpaFitter>`_
class and thus allows astropy users to use `Sherpa's <http://cxc.cfa.harvard.edu/contrib/sherpa/>`_
fitting routine.
