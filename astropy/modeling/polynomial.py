# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains predefined polynomial models.
"""
from __future__ import division
import collections
import numpy as np
from ..logger import log
from . import parameters
from .core import ParametricModel, Model, SCompositeModel
from .core import _convert_input, _convert_output
from .core import LabeledInput
from .utils import poly_map_domain, comb
from .functional_models import ShiftModel

__all__ = ['Chebyshev1DModel', 'Chebyshev2DModel', 'Legendre2DModel',
           'Legendre1DModel', 'Poly1DModel', 'Poly2DModel', 'SIPModel',
           'OrthogPolyBase', 'PolynomialModel']


class PolynomialModel(ParametricModel):

    """
    Base class for all polynomial models.

    Its main purpose is to determine how many coefficients are needed
    based on the polynomial order and dimension and to provide their
    default values, names and ordering.

    """
    def __init__(self, degree, n_inputs=1, n_outputs=1, param_dim=1, **pars):
        self.deg = degree
        self._order = self.get_numcoeff(n_inputs)
        self.param_names = self._generate_coeff_names(n_inputs)
        if not pars:
            self.set_coeff(pardim=param_dim)
        else:
            p = pars.get('c0', pars.get('c0_0'))
            if isinstance(p, collections.Sequence):
                lenpars = len(p)
            else:
                lenpars = 1
            if param_dim != lenpars:
                if param_dim == 1:
                    log.info("Inferred {0} dimensions when creating a {1} model. "
                             "Resetting param_dim to {2}".format(lenpars,
                                                                self.__class__.__name__,
                                                                lenpars))
                    param_dim = lenpars
                else:
                    raise ValueError("Number of coefficient sets ({0}) does not match number "
                                     "of parameter sets ({1}).".format(lenpars, param_dim))
            self._validate_pars(**pars)
            self.set_coeff(pardim=param_dim, **pars)
        super(PolynomialModel, self).__init__(self.param_names, n_inputs=n_inputs,
                                              n_outputs=n_outputs, param_dim=param_dim)

    def _invlex(self):
        c = []
        lencoeff = self.deg + 1
        for i in range(lencoeff):
            for j in range(lencoeff):
                if i + j <= self.deg:
                    c.append((j, i))
        return c[::-1]

    def _generate_coeff_names(self, ndim):
        ncoeff = self._order
        names = []
        if ndim == 1:
            for n in range(ncoeff):
                names.append('c{0}'.format(n))
        else:
            for i in range(self.deg + 1):
                names.append('c{0}_{1}'.format(i, 0))
            for i in range(1, self.deg + 1):
                names.append('c{0}_{1}'.format(0, i))
            for i in range(1, self.deg):
                for j in range(1, self.deg):
                    if i + j < self.deg + 1:
                        names.append('c{0}_{1}'.format(i, j))
        return names

    def _validate_pars(self, **pars):
        numcoeff = self._order
        assert(len(pars) == numcoeff)

    def set_coeff(self, pardim=1, **pars):
        """
        Set default values for coefficients
        """
        if not pars:
            for name in self.param_names:
                uname = '_' + name
                if pardim == 1:
                    self.__setattr__(uname, parameters.Parameter(
                        name, 0., self, pardim))
                else:
                    self.__setattr__(uname, parameters.Parameter(
                        name, [0.] * pardim, self, pardim))
        else:
            for name in self.param_names:
                uname = '_' + name
                self.__setattr__(uname, parameters.Parameter(
                                 name, pars[name], self, pardim))

    def get_numcoeff(self, ndim):
        """
        Return the number of coefficients in one parameter set
        """
        if self.deg < 1 or self.deg > 16:
            raise ValueError("Degree of polynomial must be 1< deg < 16")
        # deg+1 is used to account for the difference between iraf using
        # degree and numpy using exact degree
        if ndim != 1:
            nmixed = comb(self.deg, ndim)
        else:
            nmixed = 0
        numc = self.deg * ndim + nmixed + 1
        return numc

class OrthogPolyBase(ParametricModel):

    """

    This is a base class for the 2D Chebyshev and Legendre models.

    The polynomials implemented here require a maximum degree in x and y.

    Parameters
    ----------

    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : list or None
        domain of the x independent variable
    y_domain : list or None
        domain of the y independent variable
    x_window : list or None
        range of the x independent variable
    y_window : list or None
        range of the y independent variable
    param_dim : int
        number of parameter sets
    **pars : dict
        {keyword: value} pairs, representing {parameter_name: value}
    """
    def __init__(self, x_degree, y_degree, x_domain=None, x_window=None, y_domain=None,
                 y_window=None, param_dim=1, **pars):
        self.xdeg = x_degree
        self.ydeg = y_degree
        self._order = self.get_numcoeff()
        self.x_domain = x_domain
        self.y_domain = y_domain
        self.x_window = x_window
        self.y_window = y_window
        self.param_names = self._generate_coeff_names()

        if not pars:
            self.set_coeff(pardim=param_dim)

        else:
            p = pars.get('c0_0')
            if isinstance(p, collections.Sequence):
                lenpars = len(p)
            else:
                lenpars = 1
            if param_dim != lenpars:
                if param_dim == 1:
                    log.info("Inferred {0} dimensions when creating a {1} model. "
                             "Resetting param_dim to {2}".format(lenpars,
                                                                self.__class__.__name__,
                                                                lenpars))
                    param_dim = lenpars
                else:
                    raise ValueError("Number of coefficient sets {0} does not match number "
                                     "of parameter sets {1}".format(lenpars, param_dim))
            self._validate_pars(**pars)
            self.set_coeff(pardim=param_dim, **pars)
        super(OrthogPolyBase, self).__init__(self.param_names, n_inputs=2, n_outputs=1,
                                             param_dim=param_dim)

    def _generate_coeff_names(self):
        names = []
        for j in range(self.ydeg + 1):
            for i in range(self.xdeg + 1):
                names.append('c{0}_{1}'.format(i, j))
        return names

    def set_coeff(self, pardim=1, **pars):
        if not pars:
            for name in self.param_names:
                uname = '_' + name
                self.__setattr__(uname, parameters.Parameter(
                                 name, [0.] * pardim, self, pardim))
        else:
            for name in self.param_names:
                uname = '_' + name
                self.__setattr__(uname, parameters.Parameter(
                    name, pars[name], self, pardim))

    def get_numcoeff(self):
        """
        Determine how many coefficients are needed

        Returns
        -------
        numc : int
            number of coefficients

        """
        numc = (self.xdeg + 1) * (self.ydeg + 1)
        return numc

    def _validate_pars(self, **pars):
        numcoeff = self.get_numcoeff()
        assert(len(pars) == numcoeff)

    def _invlex(self):
        c = []
        xvar = np.arange(self.xdeg + 1)
        yvar = np.arange(self.ydeg + 1)
        for j in yvar:
            for i in xvar:
                c.append((i, j))
        return np.array(c[::-1])

    def invlex_coeff(self):
        coeff = []
        xvar = np.arange(self.xdeg + 1)
        yvar = np.arange(self.ydeg + 1)
        for j in yvar:
            for i in xvar:
                name = 'c' + str(i) + '_' + str(j)
                coeff.append(getattr(self, name))
        return np.array(coeff[::-1])

    def _alpha(self):
        invlexdeg = self._invlex()
        invlexdeg[:, 1] = invlexdeg[:, 1] + self.xdeg + 1
        nx = self.xdeg + 1
        ny = self.ydeg + 1
        alpha = np.zeros((ny * nx + 3, ny + nx))
        for n in range(len(invlexdeg)):
            alpha[n][invlexdeg[n]] = [1, 1]
            alpha[-2, 0] = 1
            alpha[-3, nx] = 1
        return alpha

    def imhorner(self, x, y, coeff):
        _coeff = list(coeff[:])
        _coeff.extend([0, 0, 0])
        alpha = self._alpha()
        r0 = _coeff[0]
        nalpha = len(alpha)

        karr = np.diff(alpha, axis=0)
        kfunc = self._fcache(x, y)
        xterms = self.xdeg + 1
        yterms = self.ydeg + 1
        nterms = xterms + yterms
        for n in range(1, nterms + 1 + 3):
            setattr(self, 'r' + str(n), 0.)

        for n in range(1, nalpha):
            k = karr[n - 1].nonzero()[0].max() + 1
            rsum = 0
            for i in range(1, k + 1):
                rsum = rsum + getattr(self, 'r' + str(i))
            val = kfunc[k - 1] * (r0 + rsum)
            setattr(self, 'r' + str(k), val)
            r0 = _coeff[n]
            for i in range(1, k):
                setattr(self, 'r' + str(i), 0.)
        result = r0
        for i in range(1, nterms + 1 + 3):
            result = result + getattr(self, 'r' + str(i))
        return result

    def _fcache(self, x, y):
        """
        To be implemented by subclasses
        """
        raise NotImplementedError("Subclasses should implement this")

    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : scalar, list or array
        y : scalar, lis or array

        """
        x, _ = _convert_input(x, self.param_dim)
        y, fmt = _convert_input(y, self.param_dim)
        assert x.shape == y.shape, \
            "Expected input arrays to have the same shape"
        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)
        invcoeff = self.invlex_coeff()

        result = self.imhorner(x, y, invcoeff)
        return _convert_output(result, fmt)

class Chebyshev1DModel(PolynomialModel):

    """

    1D Chebyshev polynomial of the 1st kind.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : list or None
    window : list or None
        If None, it is set to [-1,1]
        Fitters will remap the domain to this window
    param_dim : int
        number of parameter sets
    **pars : dict
        keyword : value pairs, representing parameter_name: value

    """
    def __init__(self, degree, domain=None, window=[-1, 1], param_dim=1, **pars):
        self.domain = domain
        self.window = window
        super(Chebyshev1DModel, self).__init__(degree, n_inputs=1, n_outputs=1,
                                               param_dim=param_dim, **pars)

    def clenshaw(self, x, coeff):
        """
        Evaluates the polynomial using Clenshaw's algorithm.
        """
        if isinstance(x, tuple) or isinstance(x, list):
            x = np.asarray(x)
        if len(coeff) == 1:
            c0 = coeff[0]
            c1 = 0
        elif len(coeff) == 2:
            c0 = coeff[0]
            c1 = coeff[1]
        else:
            x2 = 2 * x
            c0 = coeff[-2]
            c1 = coeff[-1]
            for i in range(3, len(coeff) + 1):
                tmp = c0
                c0 = coeff[-i] - c1
                c1 = tmp + c1 * x2
        return c0 + c1 * x

    def deriv(self, pars=None, x=None, y=None):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=x.dtype)
        v[0] = x * 0 + 1
        x2 = 2 * x
        v[1] = x
        for i in range(2, self.deg + 1):
            v[i] = v[i - 1] * x2 - v[i - 2]
        return np.rollaxis(v, 0, v.ndim)

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : scalar, list or array
            input

        """
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)#self.set_domain(x)
        x, fmt = _convert_input(x, self.param_dim)
        result = self.clenshaw(x, self.param_sets)
        return _convert_output(result, fmt)


class Legendre1DModel(PolynomialModel):

    """

    1D Legendre polynomial.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : list or None
    window : list or None
        If None, it is set to [-1,1]
        Fitters will remap the domain to this window
    param_dim : int
        number of parameter sets
    **pars : dict
        keyword: value pairs, representing parameter_name: value

    """
    def __init__(self, degree, domain=None, window=[-1, 1], param_dim=1, **pars):
        self.domain = domain
        self.window = window
        super(Legendre1DModel, self).__init__(degree, n_inputs=1, n_outputs=1,
                                              param_dim=param_dim, **pars)

    def clenshaw(self, x, coeff):
        if isinstance(x, tuple) or isinstance(x, list):
            x = np.asarray(x)
        if len(coeff) == 1:
            c0 = coeff[0]
            c1 = 0
        elif len(coeff) == 2:
            c0 = coeff[0]
            c1 = coeff[1]
        else:
            nd = len(coeff)
            c0 = coeff[-2]
            c1 = coeff[-1]
            for i in range(3, len(coeff) + 1):
                tmp = c0
                nd = nd - 1
                c0 = coeff[-i] - (c1 * (nd - 1)) / nd
                c1 = tmp + (c1 * x * (2 * nd - 1)) / nd
        return c0 + c1 * x

    def deriv(self, pars=None, x=None, y=None):
        """
        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix

        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=x.dtype)
        v[0] = x * 0 + 1
        v[1] = x
        for i in range(2, self.deg + 1):
            v[i] = (v[i - 1] * x * (2 * i - 1) - v[i - 2] * (i - 1)) / i
        return np.rollaxis(v, 0, v.ndim)

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : scalar, list or array
            input

        """
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)#self.set_domain(x)
        x, fmt = _convert_input(x, self.param_dim)
        result = self.clenshaw(x, self.param_sets)
        return _convert_output(result, fmt)


class Poly1DModel(PolynomialModel):

    """

    1D Polynomial model.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : list or None
    window : list or None
        If None, it is set to [-1,1]
        Fitters will remap the domain to this window
    param_dim : int
        number of parameter sets
    **pars : dict
        keyword: value pairs, representing parameter_name: value
    """
    def __init__(self, degree,
                 domain=[-1, 1], window=[-1, 1],
                 param_dim=1, **pars):
        self.domain = domain
        self.window = window
        super(Poly1DModel, self).__init__(degree, n_inputs=1, n_outputs=1,
                                          param_dim=param_dim, **pars)

    def deriv(self, pars=None, x=None, y=None):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix

        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.deg + 1,) + x.shape, dtype=np.float)
        v[0] = x * 0 + 1
        v[1] = x
        for i in range(2, self.deg + 1):
            v[i] = v[i - 1] * x
        return np.rollaxis(v, 0, v.ndim)

    def horner(self, x, coef):
        c0 = coef[-1] + x * 0
        for i in range(2, len(coef) + 1):
            c0 = coef[-i] + c0 * x
        return c0

    def __call__(self, x):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : scalar, list or array
            input

        """
        x, fmt = _convert_input(x, self.param_dim)
        result = self.horner(x, self.param_sets)
        return _convert_output(result, fmt)


class Poly2DModel(PolynomialModel):

    """
    2D Polynomial  model.

    Represents a general polynomial of degree n:

    .. math::
        P(x,y) = c_{0_0} + c_{1_0}*x + ...+ c_{n_0}*x^n + c_{0_1}*y + ...+ c_{0_n}*y^n
        + c_{1_1}*x*y + c_{1_2}*x*y^2 + ... + c_{1_(n-1)}*x*y^{n-1}+ ... + c_{(n-1)_1}*x^{n-1}*y

    Parameters
    ----------
    degree : int
        highest power of the polynomial,
        the number of terms is degree+1
    x_domain : list or None
        domain of the x independent variable
    y_domain : list or None
        domain of the y independent variable
    x_window : list or None
        range of the x independent variable
    y_window : list or None
        range of the y independent variable
    param_dim : int
        number of parameter sets
    pars : dict
        keyword: value pairs, representing parameter_name: value

    """
    def __init__(self, degree, x_domain=[-1, 1], y_domain=[-1, 1],
                 x_window=[-1, 1], y_window=[-1, 1],
                 param_dim=1, **pars):
        super(Poly2DModel, self).__init__(degree, n_inputs=2, n_outputs=1,
                                          param_dim=param_dim, **pars)
        self.x_domain = x_domain
        self.y_domain = y_domain
        self.x_window = x_window
        self.y_window = y_window

    def mhorner(self, x, y, coeff):
        """
        Multivariate Horner's scheme

        Parameters
        --------------
        x, y : array
        coeff : array of coefficients in inverse lexical order
        """
        alpha = np.array(self._invlex())
        r0 = coeff[0]
        r1 = r0 * 0.0
        r2 = r0 * 0.0
        karr = np.diff(alpha, axis=0)
        for n in range(len(karr)):
            if karr[n, 1] != 0:
                r2 = y * (r0 + r1 + r2)
                r1 = coeff[0] * 0.
            else:
                r1 = x * (r0 + r1)
            r0 = coeff[n + 1]
        return r0 + r1 + r2

    def deriv(self, pars=None, x=None, y=None, z=None):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : ndarray
            input
        z : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.ndim == 2:
            x = x.flatten()
        if y.ndim == 2:
            y = y.flatten()
        if x.size != y.size:
            raise ValueError('Expected x and y to be of equal size')

        designx = x[:, None] ** np.arange(self.deg + 1)
        designy = y[:, None] ** np.arange(1, self.deg + 1)

        designmixed = []
        for i in range(1, self.deg):
            for j in range(1, self.deg):
                if i + j <= self.deg:
                    designmixed.append((x ** i) * (y ** j))
        designmixed = np.array(designmixed).T
        if designmixed.any():
            v = np.hstack([designx, designy, designmixed])
        else:
            v = np.hstack([designx, designy])
        return v

    def invlex_coeff(self):
        coeff = []
        lencoeff = range(self.deg + 1)
        for i in lencoeff:
            for j in lencoeff:
                if i + j <= self.deg:
                    name = 'c' + str(j) + '_' + str(i)
                    coeff.append(getattr(self, name))
        return np.array(coeff[::-1])

    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        --------------
        x : scalar, list or array
            input
        y : scalar, list or array
            input

        """
        invcoeff = self.invlex_coeff()
        x, _ = _convert_input(x, self.param_dim)
        y, fmt = _convert_input(y, self.param_dim)
        assert x.shape == y.shape, \
            "Expected input arrays to have the same shape"

        result = self.mhorner(x, y, invcoeff)
        return _convert_output(result, fmt)


class Chebyshev2DModel(OrthogPolyBase):

    """
    2D Chebyshev polynomial of the 1st kind.

    It is defined as

    .. math:: P_{n_m}(x,y) = \sum C_{n_m}  T_n(x) T_m(y)

    Parameters
    ----------

    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : list or None
        domain of the x independent variable
    y_domain : list or None
        domain of the y independent variable
    x_window : list or None
        range of the x independent variable
    y_window : list or None
        range of the y independent variable
    param_dim : int
        number of parameter sets
    pars : dict
        keyword: value pairs, representing parameter_name: value

    """
    def __init__(self, x_degree, y_degree, x_domain=None, x_window=[-1, 1],
                 y_domain=None, y_window=[-1, 1], param_dim=1, **pars):
        super(Chebyshev2DModel, self).__init__(x_degree, y_degree,
                                               x_domain=x_domain, y_domain=y_domain,
                                               x_window=x_window, y_window=y_window,
                                               param_dim=param_dim, **pars)

    def _fcache(self, x, y):
        """
        Calculate the individual Chebyshev functions once
        and store them in a dictionary to be reused.
        """
        xterms = self.xdeg + 1
        yterms = self.ydeg + 1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[xterms] = np.ones(y.shape)
        kfunc[xterms + 1] = y.copy()
        for n in range(2, xterms):
            kfunc[n] = 2 * x * kfunc[n - 1] - kfunc[n - 2]
        for n in range(xterms + 2, xterms + yterms):
            kfunc[n] = 2 * y * kfunc[n - 1] - kfunc[n - 2]
        return kfunc

    def deriv(self, pars=None, x=None, y=None, z=None):
        """
        Derivatives with respect to the coefficients.
        This is an array with Chebyshev polynomials:

        Tx0Ty0  Tx1Ty0...TxnTy0...TxnTym

        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : ndarray
            input
        z : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix

        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.flatten()
        y = y.flatten()
        xderiv = self._chebderiv1d(x, self.xdeg + 1).T
        yderiv = self._chebderiv1d(y, self.ydeg + 1).T

        ij = []
        for i in range(self.ydeg + 1):
            for j in range(self.xdeg + 1):
                ij.append(xderiv[j] * yderiv[i])
        v = np.array(ij)
        return v.T

    def _chebderiv1d(self, x, deg):
        """
        Derivative of 1D Chebyshev series
        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        d = np.empty((deg + 1, len(x)), dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0:
            x2 = 2 * x
            d[1] = x
            for i in range(2, deg + 1):
                d[i] = d[i - 1] * x2 - d[i - 2]
        return np.rollaxis(d, 0, d.ndim)




class Legendre2DModel(OrthogPolyBase):

    """
    Legendre 2D polynomial.

    Defined as:

    .. math:: P_{nm}(x,y) = C_{n_m}  L_n(x ) L_m(y)


    Parameters
    ----------

    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : list or None
        domain of the x independent variable
    y_domain : list or None
        domain of the y independent variable
    x_window : list or None
        range of the x independent variable
    y_window : list or None
        range of the y independent variable
    param_dim : int
        number of parameter sets
    pars : dict
        keyword: value pairs, representing parameter_name: value

    """
    def __init__(self, x_degree, y_degree, x_domain=None, x_window=[-1, 1],
                 y_domain=None, y_window=[-1, 1], param_dim=1, **pars):
        super(Legendre2DModel, self).__init__(x_degree, y_degree,
                                              x_domain=x_domain, y_domain=y_domain,
                                              x_window=x_window, y_window=y_window,
                                              param_dim=param_dim, **pars)

    def _fcache(self, x, y):
        """
        Calculate the individual Legendre functions once
        and store them in a dictionary to be reused.
        """
        xterms = self.xdeg + 1
        yterms = self.ydeg + 1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[xterms] = np.ones(y.shape)
        kfunc[xterms + 1] = y.copy()
        for n in range(2, xterms):
            kfunc[n] = ((2 * (n-1) + 1) * x * kfunc[n - 1] - (n-1) * kfunc[n - 2]) / n
        for n in range(2, yterms):
            kfunc[n + xterms] = ((2 * (n-1) + 1) * y * kfunc[n + xterms - 1] -
                                 (n-1) * kfunc[n + xterms - 2]) / (n)
        return kfunc

    def deriv(self, pars=None, x=None, y=None, z=None):
        """
        Derivatives with repect to the coefficients.
        This is an array with Legendre polynomials:

        Lx0Ly0  Lx1Ly0...LxnLy0...LxnLym

        Parameters
        ----------
        pars : throw away parameter
            parameter list returned by non-linear fitters
        x : ndarray
            input
        y : ndarray
            input
        z : throw away parameter
            Present here so that the non-linear fitting algorithms can work

        Returns
        -------
        result : ndarray
            The Vandermonde matrix

        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.flatten()
        y = y.flatten()
        xderiv = self._legendderiv1d(x, self.xdeg + 1).T
        yderiv = self._legendderiv1d(y, self.ydeg + 1).T

        ij = []
        for i in range(self.ydeg + 1):
            for j in range(self.xdeg + 1):
                ij.append(xderiv[j] * yderiv[i])

        v = np.array(ij)
        return v.T

    def _legendderiv1d(self, x, deg):
        """
        Derivative of 1D Legendre polynomial
        """
        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        d = np.empty((deg + 1,) + x.shape, dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0:
            d[1] = x
            for i in range(2, deg + 1):
                d[i] = (d[i - 1] * x * (2 * i - 1) - d[i - 2] * (i - 1)) / i
        return np.rollaxis(d, 0, d.ndim)


class _SIP1D(Model):

    """
    This implements the Simple Imaging Protocol Model (SIP) in 1D.

    It's unlikely it will be used in 1D so this class is private
    and SIPModel should be used instead.

    """
    def __init__(self, order, coeff_prefix, param_dim=1, **pars):
        self.order = order
        self.coeff_prefix = coeff_prefix
        self.param_names = self._generate_coeff_names(coeff_prefix)

        if not pars:
            self.set_coeff(pardim=param_dim)
        else:
            p = pars.get('{0}0_2'.format(coeff_prefix, None))
            if isinstance(p, collections.Sequence):
                lenpars = len(p)
            else:
                lenpars = 1
            if param_dim != lenpars:
                if param_dim == 1:
                    log.info("Inferred {0} dimensions when creating a {1} model. "
                             "Resetting param_dim to {2}".format(lenpars,
                                                                self.__class__.__name__,
                                                                lenpars))
                    param_dim = lenpars
                else:
                    raise ValueError("Number of coefficient sets ({0}) does not match number "
                                     "of parameter sets ({1}).".format(lenpars, param_dim))
            self._validate_pars(ndim=2, **pars)
            self.set_coeff(pardim=param_dim, **pars)

        super(_SIP1D, self).__init__(self.param_names, n_inputs=2, n_outputs=1,
                                     param_dim=param_dim)

    def __repr__(self):
        fmt = """
        Model: {0}
        order: {1}
        coeff_prefix: {2}
        param_dim: {3}
        pars: {4}
        """.format(
              self.__class__.__name__,
              self.order,
              self.coeff_prefix,
              self._param_dim,
              "\n                   ".join(i + ':  ' + str(getattr(self, i)) for
                                           i in self.param_names)
        )
        return fmt

    def __str__(self):
        fmt = """
        Model: {0}
        order: {1}
        coeff_prefix: {2}
        Parameter sets: {3}
        Parameters:
                   {3}
        """.format(
              self.__class__.__name__,
              self.order,
              self.coeff_prefix,
              self.param_dim,
              "\n                   ".join(i + ':  ' + str(getattr(self, i)) for
                                           i in self.param_names)
        )
        return fmt

    def get_numcoeff(self, ndim):
        """
        Return the number of coefficients in one parset
        """
        if self.order < 2 or self.order > 9:
            raise ValueError("Degree of polynomial must be 2< deg < 9")
        nmixed = comb(self.order - 1, ndim)
        numc = self.order * ndim + nmixed + 1
        return numc

    def _generate_coeff_names(self, coeff_prefix):
        names = []
        for i in range(2, self.order + 1):
            names.append('{0}{1}_{2}'.format(coeff_prefix, i, 0))
        for i in range(2, self.order + 1):
            names.append('{0}{1}_{2}'.format(coeff_prefix, 0, i))
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    names.append('{0}{1}_{2}'.format(coeff_prefix, i, j))
        return names

    def set_coeff(self, pardim=1, **pars):
        if not pars:
            # default values
            for name in self.param_names:
                if pardim == 1:
                    self.__setattr__('_' + name,
                                     parameters.Parameter(name, 0, self, 1))
                else:
                    self.__setattr__('_' + name,
                                     parameters.Parameter(name, [0] * pardim,
                                                          self, pardim))
        else:
            for name in self.param_names:
                self.__setattr__('_' + name,
                                 parameters.Parameter(name, pars[name],
                                                      self, pardim))

    def _validate_pars(self, ndim, **pars):
        numcoeff = self.get_numcoeff(ndim)
        assert(len(pars) == numcoeff)

    def _coef_matrix(self, coeff_prefix):
        mat = np.zeros((self.order + 1, self.order + 1))
        for i in range(2, self.order + 1):
            mat[i, 0] = getattr(self, '{0}{1}_{2}'.format(coeff_prefix, i, 0))[0]
        for i in range(2, self.order + 1):
            mat[0, i] = getattr(self, '{0}{1}_{2}'.format(coeff_prefix, 0, i))[0]
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    mat[i, j] = getattr(self, '{0}{1}_{2}'.format(coeff_prefix, i, j))[0]
        return mat

    def _eval_sip(self, x, y, coef):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if self.coeff_prefix == 'a':
            result = np.zeros(x.shape)
        else:
            result = np.zeros(y.shape)

        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]):
                if i + j > 1 and i + j < self.order + 1:
                    result = result + coef[i, j] * x ** i * y ** j
        return result

    def __call__(self, x, y):
        mcoef = self._coef_matrix(self.coeff_prefix)
        return self._eval_sip(x, y, mcoef)


class SIPModel(SCompositeModel):

    """

    Simple Imaging Polynomial (SIP) model.

    The SIP convention is used to represent distortions in FITS image headers.
    See [1]_ for a description of the SIP convention.

    Parameters
    ----------
    crpix : list or ndarray of length(2)
        CRPIX values
    a_order : int
        SIP polynomial order for first axis
    a_coeff : dict
        SIP coefficients for first axis
    coeff_prefix : string: 'a', 'b', 'A' or 'B'
        SIP coefficient preffix
    b_order : int
        SIP order for second axis
    b_coeff : dict
        SIP coefficients for the second axis
    a_inv_order : int
        order for the inverse transformation (AP coefficients)
    a_inv_coeff : dict
        coefficients for the inverse transform
    b_inv_order : int
        order for the inverse transformation (BP coefficients)
    b_inv_coeff : dict
        coefficients for the inverse transform
    param_dim : int
        number of parameter sets

    References
    ----------
    .. [1] `David Shupe, et al, ADASS, ASP Conference Series, Vol. 347, 2005 <http://adsabs.harvard.edu/abs/2005ASPC..347..491S>`_

    """
    def __init__(self, crpix, a_order, a_coeff, b_order, b_coeff,
        a_inv_order=None, a_inv_coeff=None,b_inv_order=None, b_inv_coeff=None,
        param_dim=1):
        self.shifta = ShiftModel(-crpix[0])
        self.shiftb = ShiftModel(-crpix[1])
        self.sip1da = _SIP1D(a_order, coeff_prefix='A',
                             param_dim=param_dim, **a_coeff)
        self.sip1db = _SIP1D(b_order, coeff_prefix='B',
                             param_dim=param_dim, **b_coeff)
        if a_inv_order is not None and a_inv_coeff is not None and \
            b_inv_order is not None and b_inv_coeff is not None:
            self.inversea = _SIP1D(a_inv_order, coeff_prefix='A', **a_inv_coeff)
            self.inverseb = _SIP1D(b_inv_order, coeff_prefix='BP', **b_inv_coeff)
            self.inverse = True
        else:
            self.inverse = None
        super(SIPModel, self).__init__([self.shifta, self.shiftb, self.sip1da, self.sip1db],
                                       inmap=[['x'], ['y'], ['x', 'y'], ['x', 'y']],
                                       outmap=[['x'], ['y'], ['x1'], ['y1']], n_inputs=2,
                                       n_outputs=2)

    def __repr__(self):
        models = [self.shifta, self.shiftb, self.sip1da, self.sip1db]
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__, self.sip1d.coeff_prefix)
        fmt1 = " %s  " * len(models) % tuple([repr(model) for model in models])
        fmt = fmt + fmt1
        return fmt

    def __str__(self):
        models = [self.shifta, self.shiftb, self.sip1da, sip1db]
        fmt = """
            Model:  {0}
            """.format(self.__class__.__name__, self.sip1d.coeff_prefix)
        fmt1 = " %s  " * len(models) % tuple([str(model) for model in models])
        fmt = fmt + fmt1
        return fmt

    def __call__(self, x, y):
        """
        Transforms data using this model.

        Parameters
        ----------
        x : scalar, list ot array
            input
        y : scalar, list or array
            input

        """
        labeled_input = LabeledInput([x, y], ['x', 'y'])
        result = SCompositeModel.__call__(self, labeled_input)
        return result.x1, result.y1

