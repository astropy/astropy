# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tabular models.

Tabular models of any dimension can be created using `tabular_model`.
For convenience `Tabular1D` and `Tabular2D` are provided.

Examples
--------
>>> table = np.array([[ 3.,  0.,  0.],
...                  [ 0.,  2.,  0.],
...                  [ 0.,  0.,  0.]])
>>> points = ([1, 2, 3], [1, 2, 3])
>>> t2 = Tabular2D(points, lookup_table=table, bounds_error=False, fill_value=None, method='nearest')

"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import abc
import numpy as np
from .core import Model
from ..utils import minversion
from ..extern.six.moves import range

try:
    import scipy
    from scipy.interpolate import interpn
    has_scipy = True
except ImportError:
    has_scipy = False

has_scipy = has_scipy and minversion(scipy, "0.14")


__all__ = ['tabular_model', 'Tabular1D', 'Tabular2D']


__doctest_requires__ = {('tabular_model'): ['scipy']}


class _Tabular(Model):
    """
    Returns an interpolated lookup table value.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, ), optional
        The points defining the regular grid in n dimensions.
    lookup_table : array-like, shape (m1, ..., mn, ...)
        The data on a regular grid in n dimensions.
    method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data. Default is "linear".
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then ``fill_value`` is used.
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
    Uses `scipy.interpolate.interpn`.

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

        n_models = kwargs.get('n_models', 1)
        if n_models > 1:
            raise NotImplementedError('Only n_models=1 is supported.')
        super(_Tabular, self).__init__(**kwargs)

        if lookup_table is not None:
            lookup_table = np.asarray(lookup_table)
            if self.lookup_table.ndim != lookup_table.ndim:
                raise ValueError("lookup_table should be an array with "
                                 "{0} dimensions".format(self.lookup_table.ndim))
            self.lookup_table = lookup_table
        if points is None:
            self.points = tuple(np.arange(x, dtype=np.float)
                                for x in self.lookup_table.shape)
        else:
            if self.lookup_table.ndim == 1 and not isinstance(points, tuple):
                self.points = (points,)
            else:
                self.points = points
            if len(self.points) != self.lookup_table.ndim:
                raise ValueError("Expected grid points in "
                                 "{0} directions, got {1}".format(self.lookup_table.ndim,
                                                                  len(self.points)))

        self.bounds_error = bounds_error
        self.method = method
        self.fill_value = fill_value

    def __repr__(self):
        fmt = "<{0}(points={1}, lookup_table={2})>".format(self.__class__.__name__,
                                                           self.points, self.lookup_table)
        return fmt

    def __str__(self):
        default_keywords = [
            ('Model', self.__class__.__name__),
            ('Name', self.name),
            ('Inputs', self.inputs),
            ('Outputs', self.outputs),
            ('Parameters', ""),
            ('  points', self.points),
            ('  lookup_table', self.lookup_table),
            ('  method', self.method),
            ('  fill_value', self.fill_value),
            ('  bounds_error', self.bounds_error)
        ]

        parts = ['{0}: {1}'.format(keyword, value)
                 for keyword, value in default_keywords
                 if value is not None]

        return '\n'.join(parts)


    @property
    def bounding_box(self):
        """
        Tuple defining the default ``bounding_box`` limits,
        ``(points_low, points_high)``.

        Examples
        --------
        >>> from astropy.modeling.models import Tabular1D, Tabular2D
        >>> t1 = Tabular1D(points=[1,2,3], lookup_table=[10, 20, 30])
        >>> t1.bounding_box
        (1, 3)
        >>> t2 = Tabular2D(points=[[1,2,3],[2,3,4]], lookup_table=[[10,20,30],[20,30,40]])
        >>> t2.bounding_box
        ((2, 4), (1, 3))

        """
        bbox = [(min(p), max(p)) for p in self.points][::-1]
        if len(bbox) == 1:
            bbox = bbox[0]
        return tuple(bbox)

    def evaluate(self, *inputs):
        """
        Return the interpolated values at the input coordinates.

        Parameters
        ----------
        inputs : list of scalars or ndarrays
            Input coordinates. The number of inputs must be equal
            to the dimensions of the lookup table.
        """
        inputs = [inp.flatten() for inp in inputs[: self.n_inputs]]
        inputs = np.array(inputs).T
        if not has_scipy:
            raise ImportError("This model requires scipy >= v0.14")
        return interpn(self.points, self.lookup_table, inputs,
                       method=self.method, bounds_error=self.bounds_error,
                       fill_value=self.fill_value)


def tabular_model(dim, name=None):
    """
    Make a ``Tabular`` model where ``n_inputs`` is
    based on the dimension of the lookup_table.

    This model has to be further initialized and when evaluated
    returns the interpolated values.

    Parameters
    ----------
    dim : int
        Dimensions of the lookup table.
    name : str
        Name for the class.

    Examples
    --------
    >>> table = np.array([[ 3.,  0.,  0.],
    ...                  [ 0.,  2.,  0.],
    ...                  [ 0.,  0.,  0.]])

    >>> tab = tabular_model(2, name='Tabular2D')
    >>> print(tab)
        <class 'abc.Tabular2D'>
        Name: Tabular2D
        Inputs: (u'x0', u'x1')
        Outputs: (u'y',)

    >>> points = ([1, 2, 3], [1, 2, 3])

    Setting fill_value to None, allows extrapolation.
    >>> m = tab(points, lookup_table=table, name='my_table', bounds_error=False, fill_value=None, method='nearest')

    >>> xinterp = [0, 1, 1.5, 2.72, 3.14]
    >>> m(xinterp, xinterp)
        array([ 3., 3., 3., 0., 0.])

    """
    table = np.zeros([2] * dim)
    inputs = tuple('x{0}'.format(idx) for idx in range(table.ndim))
    members = {'lookup_table': table, 'inputs': inputs}
    if name is None:
        model_id = _Tabular._id
        _Tabular._id += 1
        name = 'Tabular{0}'.format(model_id)

    return type(str(name), (_Tabular,), members)


Tabular1D = tabular_model(1, name='Tabular1D')


Tabular2D = tabular_model(2, name='Tabular2D')

_tab_docs = """
    method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data. Default is "linear".
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then ``fill_value`` is used.
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
    Uses `scipy.interpolate.interpn`.
"""


Tabular1D.__doc__  =  """
    Tabular model in 1D.
    Returns an interpolated lookup table value.

    Parameters
    ----------
    points : array-like of float of ndim=1.
        The points defining the regular grid in n dimensions.
    lookup_table : array-like, of ndim=1.
        The data in one dimensions.
""" + _tab_docs


Tabular2D.__doc__  =  """
    Tabular model in 2D.
    Returns an interpolated lookup table value.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, m2), optional
        The points defining the regular grid in n dimensions.
    lookup_table : array-like, shape (m1, m2)
        The data on a regular grid in 2 dimensions.

""" + _tab_docs
