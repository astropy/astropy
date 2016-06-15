# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tabular models."""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import abc
import numpy as np
from .core import Model
from ..utils import minversion

try:
    import scipy
    from scipy.interpolate import interpn
    has_scipy = True
except ImportError:
    has_scipy = False

if has_scipy and not minversion(scipy, "0.14"):
    has_scipy = False


__all__ = ['tabular_model']

__doctest_requires__ = {('Tabular', 'tabular_model'): ['scipy']}


class Tabular(Model):
    """
    Returns an interpolated lookup table value.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, ), optional
        The points defining the regular grid in n dimensions.
    lookup_table : array_like, shape (m1, ..., mn, ...)
        The data on a regular grid in n dimensions.
    method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data. Default is "linear".
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.
    fill_value : float, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.  Extrapolation is not supported by method
        "splinef2d".

    Returns
    -------
    value : ndarray
        Interpolated values at input coordinates.

    Raises
    ------
    ImportError
        Scipy is not installed.

    Notes
    -----
    Uses ``scipy.interpolate.interpn``

    """

    linear = False
    fittable = False

    standard_broadcasting = False
    outputs = ('y',)

    lookup_table = abc.abstractproperty()

    _is_dynamic = True

    _id = 0

    def __init__(self, points=None, lookup_table=None, method='linear',
                 bounds_error=True, fill_value=np.nan, **kwargs):

        super(Tabular, self).__init__(**kwargs)
        if lookup_table is not None:
            lookup_table = np.asarray(lookup_table)
            if self.lookup_table.ndim != lookup_table.ndim:
                raise ValueError("lookup_table should be an array with "
                                 "{0} dimensions".format(self.lookup_table.ndim))
            self.lookup_table = lookup_table
        if points is None:
            self._points = tuple(np.arange(x, dtype=np.float)
                                 for x in self.lookup_table.shape)
        else:
            if self.lookup_table.ndim == 1:
                self._points = (points,)
            else:
                self._points = points
            if len(self._points) != self.lookup_table.ndim:
                raise ValueError("Expected grid points in "
                                 "{0} directions, got {1}".format(self.lookup_table.ndim,
                                                                  len(self._points)))

        self.bounds_error = bounds_error
        self.method = method
        self.fill_value = fill_value

    @property
    def points(self):
        return self._points

    def evaluate(self, *inputs):
        """
        Return the interpolated values at the input coordinates.

        Parameters
        ----------
        inputs : list of scalars or ndarrays
            Input coordinates. The number of inputs must be equal
            to the dimensions of the lookup table.
        """
        inputs = np.array(inputs[: self.n_inputs]).T
        if not has_scipy:
            raise ImportError("This model requires scipy >= v0.14")
        return interpn(self.points, self.lookup_table, inputs,
                       method=self.method, bounds_error=self.bounds_error,
                       fill_value=self.fill_value)


def tabular_model(lookup_table, name=None):
    """
    Make a ``Tabular`` model where ``n_inputs`` is
    based on the dimension of the lookup_table.

    This model has to be further initialized and when evaluated
    returns the interpolated values.

    Parameters
    ----------
    lookup_table : array_like, shape (m1, ..., mn)
        The data on a regular grid in n dimensions.
    name : str
        Name for the class.

    Examples
    --------
    >>> table = np.array([[ 3.,  0.,  0.],
    ...                  [ 0.,  2.,  0.],
    ...                  [ 0.,  0.,  0.]])

    >>> tab = tabular_model(table, name='Tabular2D')
    >>> print(tab)
        <class 'abc.Tabular2D'>
        Name: Tabular2D
        Inputs: (u'x0', u'x1')
        Outputs: (u'y',)

    >>> points = ([1, 2, 3], [1, 2, 3])

    Setting fill_value to None, allows extrapolation.
    >>> m = tab(points, name='my_table', bounds_error=False, fill_value=None, method='nearest')

    >>> xinterp = [0, 1, 1.5, 2.72, 3.14]
    >>> m(xinterp, xinterp)
        array([ 3., 3., 3., 0., 0.])

    """
    table = np.asanyarray(lookup_table)
    inputs = tuple('x{0}'.format(idx) for idx in range(table.ndim))
    members = {'lookup_table': table, 'inputs': inputs}
    if name is None:
        model_id = Tabular._id
        Tabular._id += 1
        name = 'Tabular{0}'.format(model_id)

    return type(str(name), (Tabular,), members)
