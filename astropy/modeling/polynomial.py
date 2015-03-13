# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains predefined polynomial models.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from .core import FittableModel, Model
from .functional_models import Shift
from .parameters import Parameter
from .utils import poly_map_domain, comb
from ..utils import lazyproperty, indent


__all__ = [
    'Chebyshev1D', 'Chebyshev2D', 'InverseSIP',
    'Legendre1D', 'Legendre2D', 'Polynomial1D',
    'Polynomial2D', 'SIP', 'OrthoPolynomialBase',
    'PolynomialModel'
]


class PolynomialBase(FittableModel):
    """
    Base class for all polynomial-like models with an arbitrary number of
    parameters in the form of coefficients.

    In this case Parameter instances are returned through the class's
    ``__getattr__`` rather than through class descriptors.
    """

    # Default _param_names list; this will be filled in by the implementation's
    # __init__
    _param_names = ()

    linear = True
    col_fit_deriv = False

    @lazyproperty
    def param_names(self):
        """Coefficient names generated based on the model's polynomial degree
        and number of dimensions.

        Subclasses should implement this to return parameter names in the
        desired format.

        On most `Model` classes this is a class attribute, but for polynomial
        models it is an instance attribute since each polynomial model instance
        can have different parameters depending on the degree of the polynomial
        and the number of dimensions, for example.
        """

        return self._param_names

    def __getattr__(self, attr):
        if self._param_names and attr in self._param_names:
            return Parameter(attr, default=0.0, model=self)

        raise AttributeError(attr)

    def __setattr__(self, attr, value):
        # TODO: Support a means of specifying default values for coefficients
        # Check for self._ndim first--if it hasn't been defined then the
        # instance hasn't been initialized yet and self.param_names probably
        # won't work.
        # This has to vaguely duplicate the functionality of
        # Parameter.__set__.
        # TODO: I wonder if there might be a way around that though...
        if attr[0] != '_' and self._param_names and attr in self._param_names:
            param = Parameter(attr, default=0.0, model=self)
            # This is a little hackish, but we can actually reuse the
            # Parameter.__set__ method here
            param.__set__(self, value)
        else:
            super(PolynomialBase, self).__setattr__(attr, value)


class PolynomialModel(PolynomialBase):
    """
    Base class for polynomial models.

    Its main purpose is to determine how many coefficients are needed
    based on the polynomial order and dimension and to provide their
    default values, names and ordering.
    """

    def __init__(self, degree, n_models=None, model_set_axis=None,
                 name=None, meta=None, **params):
        self._degree = degree
        self._order = self.get_num_coeff(self.n_inputs)
        self._param_names = self._generate_coeff_names(self.n_inputs)

        super(PolynomialModel, self).__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name,
            meta=meta, **params)

    def __repr__(self):
        return self._format_repr([self.degree])

    def __str__(self):
        return self._format_str([('Degree', self.degree)])

    @property
    def degree(self):
        """Degree of polynomial."""

        return self._degree

    def get_num_coeff(self, ndim):
        """
        Return the number of coefficients in one parameter set
        """

        if self.degree < 0:
            raise ValueError("Degree of polynomial must be positive or null")
        # deg+1 is used to account for the difference between iraf using
        # degree and numpy using exact degree
        if ndim != 1:
            nmixed = comb(self.degree, ndim)
        else:
            nmixed = 0
        numc = self.degree * ndim + nmixed + 1
        return numc

    def _invlex(self):
        c = []
        lencoeff = self.degree + 1
        for i in range(lencoeff):
            for j in range(lencoeff):
                if i + j <= self.degree:
                    c.append((j, i))
        return c[::-1]

    def _generate_coeff_names(self, ndim):
        names = []
        if ndim == 1:
            for n in range(self._order):
                names.append('c{0}'.format(n))
        else:
            for i in range(self.degree + 1):
                names.append('c{0}_{1}'.format(i, 0))
            for i in range(1, self.degree + 1):
                names.append('c{0}_{1}'.format(0, i))
            for i in range(1, self.degree):
                for j in range(1, self.degree):
                    if i + j < self.degree + 1:
                        names.append('c{0}_{1}'.format(i, j))
        return tuple(names)


class OrthoPolynomialBase(PolynomialBase):
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
    x_window : list or None
        range of the x independent variable
    y_domain : list or None
        domain of the y independent variable
    y_window : list or None
        range of the y independent variable
    param_dim : int
        number of parameter sets
    **params : dict
        {keyword: value} pairs, representing {parameter_name: value}
    """

    inputs = ('x', 'y')
    outputs = ('z',)

    def __init__(self, x_degree, y_degree, x_domain=None, x_window=None,
                 y_domain=None, y_window=None, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        # TODO: Perhaps some of these other parameters should be properties?
        # TODO: An awful lot of the functionality in this method is still
        # shared by PolynomialModel; perhaps some of it can be generalized in
        # PolynomialBase
        self.x_degree = x_degree
        self.y_degree = y_degree
        self._order = self.get_num_coeff()
        self.x_domain = x_domain
        self.y_domain = y_domain
        self.x_window = x_window
        self.y_window = y_window
        self._param_names = self._generate_coeff_names()

        super(OrthoPolynomialBase, self).__init__(
            n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def __repr__(self):
        return self._format_repr([self.x_degree, self.y_degree])

    def __str__(self):
        return self._format_str(
            [('X-Degree', self.x_degree),
             ('Y-Degree', self.y_degree)])

    def get_num_coeff(self):
        """
        Determine how many coefficients are needed

        Returns
        -------
        numc : int
            number of coefficients
        """

        return (self.x_degree + 1) * (self.y_degree + 1)

    def _invlex(self):
        # TODO: This is a very slow way to do this; fix it and related methods
        # like _alpha
        c = []
        xvar = np.arange(self.x_degree + 1)
        yvar = np.arange(self.y_degree + 1)
        for j in yvar:
            for i in xvar:
                c.append((i, j))
        return np.array(c[::-1])

    def invlex_coeff(self):
        coeff = []
        xvar = np.arange(self.x_degree + 1)
        yvar = np.arange(self.y_degree + 1)
        for j in yvar:
            for i in xvar:
                name = 'c{0}_{1}'.format(i, j)
                coeff.append(getattr(self, name))
        return np.array(coeff[::-1])

    def _alpha(self):
        invlexdeg = self._invlex()
        invlexdeg[:, 1] = invlexdeg[:, 1] + self.x_degree + 1
        nx = self.x_degree + 1
        ny = self.y_degree + 1
        alpha = np.zeros((ny * nx + 3, ny + nx))
        for n in range(len(invlexdeg)):
            alpha[n][invlexdeg[n]] = [1, 1]
            alpha[-2, 0] = 1
            alpha[-3, nx] = 1
        return alpha

    def imhorner(self, x, y, coeff):
        _coeff = list(coeff)
        _coeff.extend([0, 0, 0])
        alpha = self._alpha()
        r0 = _coeff[0]
        nalpha = len(alpha)

        karr = np.diff(alpha, axis=0)
        kfunc = self._fcache(x, y)
        x_terms = self.x_degree + 1
        y_terms = self.y_degree + 1
        nterms = x_terms + y_terms
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

    def _generate_coeff_names(self):
        names = []
        for j in range(self.y_degree + 1):
            for i in range(self.x_degree + 1):
                names.append('c{0}_{1}'.format(i, j))
        return tuple(names)

    def _fcache(self, x, y):
        # TODO: Write a docstring explaining the actual purpose of this method
        """To be implemented by subclasses"""

        raise NotImplementedError("Subclasses should implement this")

    def evaluate(self, x, y, *coeffs):
        invcoeff = self.invlex_coeff()
        return self.imhorner(x, y, invcoeff)

    def prepare_inputs(self, x, y, **kwargs):
        inputs, format_info = \
                super(OrthoPolynomialBase, self).prepare_inputs(x, y, **kwargs)

        x, y = inputs

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")
        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)

        return (x, y), format_info


class Chebyshev1D(PolynomialModel):
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
    **params : dict
        keyword : value pairs, representing parameter_name: value
    """

    inputs = ('x',)
    outputs = ('y',)

    def __init__(self, degree, domain=None, window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.domain = domain
        self.window = window
        super(Chebyshev1D, self).__init__(
            degree, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """

        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.degree + 1,) + x.shape, dtype=x.dtype)
        v[0] = 1
        if self.degree > 0:
            x2 = 2 * x
            v[1] = x
            for i in range(2, self.degree + 1):
                v[i] = v[i - 1] * x2 - v[i - 2]
        return np.rollaxis(v, 0, v.ndim)

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = \
                super(PolynomialModel, self).prepare_inputs(x, **kwargs)

        x = inputs[0]

        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)

        return (x,), format_info

    @classmethod
    def evaluate(cls, x, *coeffs):
        return cls.clenshaw(x, coeffs)

    @staticmethod
    def clenshaw(x, coeffs):
        """Evaluates the polynomial using Clenshaw's algorithm."""

        if len(coeffs) == 1:
            c0 = coeffs[0]
            c1 = 0
        elif len(coeffs) == 2:
            c0 = coeffs[0]
            c1 = coeffs[1]
        else:
            x2 = 2 * x
            c0 = coeffs[-2]
            c1 = coeffs[-1]
            for i in range(3, len(coeffs) + 1):
                tmp = c0
                c0 = coeffs[-i] - c1
                c1 = tmp + c1 * x2
        return c0 + c1 * x


class Legendre1D(PolynomialModel):
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
    **params : dict
        keyword: value pairs, representing parameter_name: value
    """

    inputs = ('x',)
    outputs = ('y',)

    def __init__(self, degree, domain=None, window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.domain = domain
        self.window = window
        super(Legendre1D, self).__init__(
            degree, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = \
                super(PolynomialModel, self).prepare_inputs(x, **kwargs)

        x = inputs[0]

        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)

        return (x,), format_info

    @classmethod
    def evaluate(cls, x, *coeffs):
        return cls.clenshaw(x, coeffs)

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """

        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        v = np.empty((self.degree + 1,) + x.shape, dtype=x.dtype)
        v[0] = 1
        if self.degree > 0:
            v[1] = x
            for i in range(2, self.degree + 1):
                v[i] = (v[i - 1] * x * (2 * i - 1) - v[i - 2] * (i - 1)) / i
        return np.rollaxis(v, 0, v.ndim)

    @staticmethod
    def clenshaw(x, coeffs):
        if len(coeffs) == 1:
            c0 = coeffs[0]
            c1 = 0
        elif len(coeffs) == 2:
            c0 = coeffs[0]
            c1 = coeffs[1]
        else:
            nd = len(coeffs)
            c0 = coeffs[-2]
            c1 = coeffs[-1]
            for i in range(3, len(coeffs) + 1):
                tmp = c0
                nd = nd - 1
                c0 = coeffs[-i] - (c1 * (nd - 1)) / nd
                c1 = tmp + (c1 * x * (2 * nd - 1)) / nd
        return c0 + c1 * x


class Polynomial1D(PolynomialModel):
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
    **params : dict
        keyword: value pairs, representing parameter_name: value
    """

    inputs = ('x',)
    outputs = ('y',)

    def __init__(self, degree, domain=[-1, 1], window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.domain = domain
        self.window = window
        super(Polynomial1D, self).__init__(
            degree, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)

    def prepare_inputs(self, x, **kwargs):
        inputs, format_info = \
                super(Polynomial1D, self).prepare_inputs(x, **kwargs)

        x = inputs[0]

        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)

        return (x,), format_info

    @classmethod
    def evaluate(cls, x, *coeffs):
        return cls.horner(x, coeffs)

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """

        v = np.empty((self.degree + 1,) + x.shape, dtype=np.float)
        v[0] = 1
        if self.degree > 0:
            v[1] = x
            for i in range(2, self.degree + 1):
                v[i] = v[i - 1] * x
        return np.rollaxis(v, 0, v.ndim)

    @staticmethod
    def horner(x, coeffs):
        c0 = coeffs[-1] + x * 0
        for i in range(2, len(coeffs) + 1):
            c0 = coeffs[-i] + c0 * x
        return c0


class Polynomial2D(PolynomialModel):
    """
    2D Polynomial  model.

    Represents a general polynomial of degree n:

    .. math::

        P(x,y) = c_{00} + c_{10}x + ...+ c_{n0}x^n + c_{01}y + ...+ c_{0n}y^n
        + c_{11}xy + c_{12}xy^2 + ... + c_{1(n-1)}xy^{n-1}+ ... + c_{(n-1)1}x^{n-1}y

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
    **params : dict
        keyword: value pairs, representing parameter_name: value
    """

    inputs = ('x', 'y')
    outputs = ('z',)

    def __init__(self, degree, x_domain=[-1, 1], y_domain=[-1, 1],
                 x_window=[-1, 1], y_window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Polynomial2D, self).__init__(
            degree, n_models=n_models, model_set_axis=model_set_axis,
            name=name, meta=meta, **params)
        self.x_domain = x_domain
        self.y_domain = y_domain
        self.x_window = x_window
        self.y_window = y_window

    def prepare_inputs(self, x, y, **kwargs):
        inputs, format_info = \
                super(Polynomial2D, self).prepare_inputs(x, y, **kwargs)

        x, y = inputs

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)

        return (x, y), format_info

    def evaluate(self, x, y, *coeffs):
        invcoeff = self.invlex_coeff(coeffs)
        return self.multivariate_horner(x, y, invcoeff)

    def fit_deriv(self, x, y, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

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

        designx = x[:, None] ** np.arange(self.degree + 1)
        designy = y[:, None] ** np.arange(1, self.degree + 1)

        designmixed = []
        for i in range(1, self.degree):
            for j in range(1, self.degree):
                if i + j <= self.degree:
                    designmixed.append((x ** i) * (y ** j))
        designmixed = np.array(designmixed).T
        if designmixed.any():
            v = np.hstack([designx, designy, designmixed])
        else:
            v = np.hstack([designx, designy])
        return v

    def invlex_coeff(self, coeffs):
        invlex_coeffs = []
        lencoeff = range(self.degree + 1)
        for i in lencoeff:
            for j in lencoeff:
                if i + j <= self.degree:
                    name = 'c{0}_{1}'.format(j, i)
                    coeff = coeffs[self.param_names.index(name)]
                    invlex_coeffs.append(coeff)
        return np.array(invlex_coeffs[::-1])

    def multivariate_horner(self, x, y, coeffs):
        """
        Multivariate Horner's scheme

        Parameters
        --------------
        x, y : array
        coeff : array of coefficients in inverse lexical order
        """

        alpha = np.array(self._invlex())
        r0 = coeffs[0]
        r1 = r0 * 0.0
        r2 = r0 * 0.0
        karr = np.diff(alpha, axis=0)
        for n in range(len(karr)):
            if karr[n, 1] != 0:
                r2 = y * (r0 + r1 + r2)
                r1 = coeffs[0] * 0.
            else:
                r1 = x * (r0 + r1)
            r0 = coeffs[n + 1]
        return r0 + r1 + r2


class Chebyshev2D(OrthoPolynomialBase):
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
    **params : dict
        keyword: value pairs, representing parameter_name: value

    """

    def __init__(self, x_degree, y_degree, x_domain=None, x_window=[-1, 1],
                 y_domain=None, y_window=[-1,1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Chebyshev2D, self).__init__(
            x_degree, y_degree, x_domain=x_domain, y_domain=y_domain,
            x_window=x_window, y_window=y_window, n_models=n_models,
            model_set_axis=model_set_axis, name=name, meta=meta, **params)

    def _fcache(self, x, y):
        """
        Calculate the individual Chebyshev functions once and store them in a
        dictionary to be reused.
        """

        x_terms = self.x_degree + 1
        y_terms = self.y_degree + 1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[x_terms] = np.ones(y.shape)
        kfunc[x_terms + 1] = y.copy()
        for n in range(2, x_terms):
            kfunc[n] = 2 * x * kfunc[n - 1] - kfunc[n - 2]
        for n in range(x_terms + 2, x_terms + y_terms):
            kfunc[n] = 2 * y * kfunc[n - 1] - kfunc[n - 2]
        return kfunc

    def fit_deriv(self, x, y, *params):
        """
        Derivatives with respect to the coefficients.

        This is an array with Chebyshev polynomials:

        .. math::

            T_{x_0}T_{y_0}, T_{x_1}T_{y_0}...T_{x_n}T_{y_0}...T_{x_n}T_{y_m}

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """

        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")

        x = x.flatten()
        y = y.flatten()
        x_deriv = self._chebderiv1d(x, self.x_degree + 1).T
        y_deriv = self._chebderiv1d(y, self.y_degree + 1).T

        ij = []
        for i in range(self.y_degree + 1):
            for j in range(self.x_degree + 1):
                ij.append(x_deriv[j] * y_deriv[i])

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


class Legendre2D(OrthoPolynomialBase):
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
    **params : dict
        keyword: value pairs, representing parameter_name: value
    """

    def __init__(self, x_degree, y_degree, x_domain=None, x_window=[-1, 1],
                 y_domain=None, y_window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Legendre2D, self).__init__(
            x_degree, y_degree, x_domain=x_domain, y_domain=y_domain,
            x_window=x_window, y_window=y_window, n_models=n_models,
            model_set_axis=model_set_axis, name=name, meta=meta, **params)

    def _fcache(self, x, y):
        """
        Calculate the individual Legendre functions once and store them in a
        dictionary to be reused.
        """

        x_terms = self.x_degree + 1
        y_terms = self.y_degree + 1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = x.copy()
        kfunc[x_terms] = np.ones(y.shape)
        kfunc[x_terms + 1] = y.copy()
        for n in range(2, x_terms):
            kfunc[n] = (((2 * (n - 1) + 1) * x * kfunc[n - 1] -
                        (n - 1) * kfunc[n - 2]) / n)
        for n in range(2, y_terms):
            kfunc[n + x_terms] = ((2 * (n - 1) + 1) * y * kfunc[n + x_terms - 1] -
                                  (n - 1) * kfunc[n + x_terms - 2]) / (n)
        return kfunc

    def fit_deriv(self, x, y, *params):
        """
        Derivatives with respect to the coefficients.
        This is an array with Legendre polynomials:

        Lx0Ly0  Lx1Ly0...LxnLy0...LxnLym

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        params : throw away parameter
            parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.flatten()
        y = y.flatten()
        x_deriv = self._legendderiv1d(x, self.x_degree + 1).T
        y_deriv = self._legendderiv1d(y, self.y_degree + 1).T

        ij = []
        for i in range(self.y_degree + 1):
            for j in range(self.x_degree + 1):
                ij.append(x_deriv[j] * y_deriv[i])

        v = np.array(ij)
        return v.T

    def _legendderiv1d(self, x, deg):
        """Derivative of 1D Legendre polynomial"""

        x = np.array(x, dtype=np.float, copy=False, ndmin=1)
        d = np.empty((deg + 1,) + x.shape, dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0:
            d[1] = x
            for i in range(2, deg + 1):
                d[i] = (d[i - 1] * x * (2 * i - 1) - d[i - 2] * (i - 1)) / i
        return np.rollaxis(d, 0, d.ndim)


class _SIP1D(PolynomialBase):
    """
    This implements the Simple Imaging Polynomial Model (SIP) in 1D.

    It's unlikely it will be used in 1D so this class is private
    and SIP should be used instead.
    """

    inputs = ('u', 'v')
    outputs = ('w',)

    def __init__(self, order, coeff_prefix, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.order = order
        self.coeff_prefix = coeff_prefix
        self._param_names = self._generate_coeff_names(coeff_prefix)

        super(_SIP1D, self).__init__(n_models=n_models,
                                     model_set_axis=model_set_axis,
                                     name=name, meta=meta, **params)

    def __repr__(self):
        return self._format_repr(args=[self.order, self.coeff_prefix])

    def __str__(self):
        return self._format_str(
            [('Order', self.order),
             ('Coeff. Prefix', self.coeff_prefix)])

    def evaluate(self, x, y, *coeffs):
        # TODO: Rewrite this so that it uses a simpler method of determining
        # the matrix based on the number of given coefficients.
        mcoef = self._coeff_matrix(self.coeff_prefix, coeffs)
        return self._eval_sip(x, y, mcoef)

    def get_num_coeff(self, ndim):
        """
        Return the number of coefficients in one param set
        """

        if self.order < 2 or self.order > 9:
            raise ValueError("Degree of polynomial must be 2< deg < 9")

        nmixed = comb(self.order, ndim)
        # remove 3 terms because SIP deg >= 2
        numc = self.order * ndim + nmixed - 2
        return numc

    def _generate_coeff_names(self, coeff_prefix):
        names = []
        for i in range(2, self.order + 1):
            names.append('{0}_{1}_{2}'.format(coeff_prefix, i, 0))
        for i in range(2, self.order + 1):
            names.append('{0}_{1}_{2}'.format(coeff_prefix, 0, i))
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    names.append('{0}_{1}_{2}'.format(coeff_prefix, i, j))
        return names

    def _coeff_matrix(self, coeff_prefix, coeffs):
        mat = np.zeros((self.order + 1, self.order + 1))
        for i in range(2, self.order + 1):
            attr = '{0}_{1}_{2}'.format(coeff_prefix, i, 0)
            mat[i, 0] = coeffs[self.param_names.index(attr)]
        for i in range(2, self.order + 1):
            attr = '{0}_{1}_{2}'.format(coeff_prefix, 0, i)
            mat[0, i] = coeffs[self.param_names.index(attr)]
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    attr = '{0}_{1}_{2}'.format(coeff_prefix, i, j)
                    mat[i, j] = coeffs[self.param_names.index(attr)]
        return mat

    def _eval_sip(self, x, y, coef):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if self.coeff_prefix == 'A':
            result = np.zeros(x.shape)
        else:
            result = np.zeros(y.shape)

        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]):
                if i + j > 1 and i + j < self.order + 1:
                    result = result + coef[i, j] * x ** i * y ** j
        return result


class SIP(Model):
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
    b_order : int
        SIP order for second axis
    a_coeff : dict
        SIP coefficients for first axis
    b_coeff : dict
        SIP coefficients for the second axis
    ap_order : int
        order for the inverse transformation (AP coefficients)
    bp_order : int
        order for the inverse transformation (BP coefficients)
    ap_coeff : dict
        coefficients for the inverse transform
    bp_coeff : dict
        coefficients for the inverse transform
    param_dim : int
        number of parameter sets

    References
    ----------
    .. [1] `David Shupe, et al, ADASS, ASP Conference Series, Vol. 347, 2005 <http://adsabs.harvard.edu/abs/2005ASPC..347..491S>`_
    """


    inputs = ('u', 'v')
    outputs = ('x', 'y')

    def __init__(self, crpix, a_order, b_order, a_coeff={}, b_coeff={},
                 ap_order=None, bp_order=None, ap_coeff={}, bp_coeff={},
                 n_models=None, model_set_axis=None, name=None, meta=None):
        self._crpix = crpix
        self._a_order = a_order
        self._b_order = b_order
        self._a_coeff = a_coeff
        self._b_coeff = b_coeff
        self._ap_order = ap_order
        self._bp_order = bp_order
        self._ap_coeff = ap_coeff
        self._bp_coeff = bp_coeff
        self.shift_a = Shift(-crpix[0])
        self.shift_b = Shift(-crpix[1])
        self.sip1d_a = _SIP1D(a_order, coeff_prefix='A', n_models=n_models,
                              model_set_axis=model_set_axis, **a_coeff)
        self.sip1d_b = _SIP1D(b_order, coeff_prefix='B', n_models=n_models,
                              model_set_axis=model_set_axis, **b_coeff)
        super(SIP, self).__init__(n_models=n_models,
                                  model_set_axis=model_set_axis, name=name,
                                  meta=meta)

    def __repr__(self):
        return '<{0}({1!r})>'.format(self.__class__.__name__,
            [self.shift_a, self.shift_b, self.sip1d_a, self.sip1d_b])

    def __str__(self):
        parts = ['Model: {0}'.format(self.__class__.__name__)]
        for model in [self.shift_a, self.shift_b, self.sip1d_a, self.sip1d_b]:
            parts.append(indent(str(model), width=4))
            parts.append('')

        return '\n'.join(parts)

    @property
    def inverse(self):
        if (self._ap_order is not None and self._bp_order is not None):
            return InverseSIP(self._ap_order, self._bp_order,
                              self._ap_coeff, self._bp_coeff)
        else:
            raise NotImplementedError("SIP inverse coefficients are not available.")

    def evaluate(self, x, y):
        u = self.shift_a.evaluate(x, *self.shift_a.param_sets)
        v = self.shift_b.evaluate(y, *self.shift_b.param_sets)
        f = self.sip1d_a.evaluate(u, v, *self.sip1d_a.param_sets)
        g = self.sip1d_b.evaluate(u, v, *self.sip1d_b.param_sets)
        return f, g


class InverseSIP(Model):
    """
    Inverse Simple Imaging Polynomial

    Parameters
    ----------
    ap_order : int
        order for the inverse transformation (AP coefficients)
    bp_order : int
        order for the inverse transformation (BP coefficients)
    ap_coeff : dict
        coefficients for the inverse transform
    bp_coeff : dict
        coefficients for the inverse transform
    param_dim : int
        number of parameter sets

    """

    inputs = ('x', 'y')
    outputs = ('u', 'v')

    def __init__(self, ap_order, bp_order, ap_coeff={}, bp_coeff={},
                 n_models=None, model_set_axis=None, name=None, meta=None):
        self._ap_order = ap_order
        self._bp_order = bp_order
        self._ap_coeff = ap_coeff
        self._bp_coeff = bp_coeff

        # define the 0th term in order to use Polynomial2D
        ap_coeff.setdefault('AP_0_0', 0)
        bp_coeff.setdefault('BP_0_0', 0)

        ap_coeff_params = dict((k.replace('AP_', 'c'), v)
                               for k, v in ap_coeff.items())
        bp_coeff_params = dict((k.replace('BP_', 'c'), v)
                               for k, v in bp_coeff.items())

        self.sip1d_ap = Polynomial2D(degree=ap_order,
                                     model_set_axis=model_set_axis,
                                     **ap_coeff_params)
        self.sip1d_bp = Polynomial2D(degree=bp_order,
                                     model_set_axis=model_set_axis,
                                     **bp_coeff_params)
        super(InverseSIP, self).__init__(n_models=n_models,
                                         model_set_axis=model_set_axis,
                                         name=name, meta=meta)

    def __repr__(self):
        return '<{0}({1!r})>'.format(self.__class__.__name__,
            [self.sip1d_ap, self.sip1d_bp])

    def __str__(self):
        parts = ['Model: {0}'.format(self.__class__.__name__)]
        for model in [self.sip1d_ap, self.sip1d_bp]:
            parts.append(indent(str(model), width=4))
            parts.append('')

        return '\n'.join(parts)

    def evaluate(self, x, y):
        x1 = self.sip1d_ap.evaluate(x, y, *self.sip1d_ap.param_sets)
        y1 = self.sip1d_bp.evaluate(x, y, *self.sip1d_bp.param_sets)
        return x1, y1
