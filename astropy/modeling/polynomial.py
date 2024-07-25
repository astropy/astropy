# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains models representing polynomials and polynomial series.

"""

# pylint: disable=invalid-name
from math import comb
from textwrap import indent

import numpy as np

from astropy.utils import check_broadcast
from astropy.utils.compat import COPY_IF_NEEDED

from .core import FittableModel, Model
from .functional_models import Shift
from .parameters import Parameter
from .utils import _validate_domain_window, poly_map_domain

__all__ = [
    "Chebyshev1D",
    "Chebyshev2D",
    "Hermite1D",
    "Hermite2D",
    "InverseSIP",
    "Legendre1D",
    "Legendre2D",
    "Polynomial1D",
    "Polynomial2D",
    "SIP",
    "OrthoPolynomialBase",
    "PolynomialModel",
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

    @property
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


class PolynomialModel(PolynomialBase):
    """
    Base class for polynomial models.

    Its main purpose is to determine how many coefficients are needed
    based on the polynomial order and dimension and to provide their
    default values, names and ordering.
    """

    def __init__(
        self, degree, n_models=None, model_set_axis=None, name=None, meta=None, **params
    ):
        self._degree = degree
        self._order = self.get_num_coeff(self.n_inputs)
        self._param_names = self._generate_coeff_names(self.n_inputs)
        if n_models:
            if model_set_axis is None:
                model_set_axis = 0
            minshape = (1,) * model_set_axis + (n_models,)
        else:
            minshape = ()
        for param_name in self._param_names:
            self._parameters_[param_name] = Parameter(
                param_name, default=np.zeros(minshape)
            )

        super().__init__(
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    @property
    def degree(self):
        """Degree of polynomial."""
        return self._degree

    def get_num_coeff(self, ndim):
        """
        Return the number of coefficients in one parameter set.
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
                names.append(f"c{n}")
        else:
            for i in range(self.degree + 1):
                names.append(f"c{i}_{0}")
            for i in range(1, self.degree + 1):
                names.append(f"c{0}_{i}")
            for i in range(1, self.degree):
                for j in range(1, self.degree):
                    if i + j < self.degree + 1:
                        names.append(f"c{i}_{j}")
        return tuple(names)


class _PolyDomainWindow1D(PolynomialModel):
    """
    This class sets ``domain`` and ``window`` of 1D polynomials.
    """

    def __init__(
        self,
        degree,
        domain=None,
        window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree, n_models, model_set_axis, name=name, meta=meta, **params
        )

        self._set_default_domain_window(domain, window)

    @property
    def window(self):
        return self._window

    @window.setter
    def window(self, val):
        self._window = _validate_domain_window(val)

    @property
    def domain(self):
        return self._domain

    @domain.setter
    def domain(self, val):
        self._domain = _validate_domain_window(val)

    def _set_default_domain_window(self, domain, window):
        """
        This method sets the ``domain`` and ``window`` attributes on 1D subclasses.

        """
        self._default_domain_window = {"domain": None, "window": (-1, 1)}
        self.window = window or (-1, 1)
        self.domain = domain

    def __repr__(self):
        return self._format_repr(
            [self.degree],
            kwargs={"domain": self.domain, "window": self.window},
            defaults=self._default_domain_window,
        )

    def __str__(self):
        return self._format_str(
            [("Degree", self.degree), ("Domain", self.domain), ("Window", self.window)],
            self._default_domain_window,
        )


class OrthoPolynomialBase(PolynomialBase):
    """
    This is a base class for the 2D Chebyshev and Legendre models.

    The polynomials implemented here require a maximum degree in x and y.

    For explanation of ``x_domain``, ``y_domain``, ```x_window`` and ```y_window``
    see :ref:`Notes regarding usage of domain and window <astropy:domain-window-note>`.


    Parameters
    ----------
    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : tuple or None, optional
        domain of the x independent variable
    x_window : tuple or None, optional
        range of the x independent variable
    y_domain : tuple or None, optional
        domain of the y independent variable
    y_window : tuple or None, optional
        range of the y independent variable
    **params : dict
        {keyword: value} pairs, representing {parameter_name: value}
    """

    n_inputs = 2
    n_outputs = 1

    def __init__(
        self,
        x_degree,
        y_degree,
        x_domain=None,
        x_window=None,
        y_domain=None,
        y_window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        self.x_degree = x_degree
        self.y_degree = y_degree
        self._order = self.get_num_coeff()
        # Set the ``x/y_domain`` and ``x/y_wndow`` attributes in subclasses.
        self._default_domain_window = {
            "x_window": (-1, 1),
            "y_window": (-1, 1),
            "x_domain": None,
            "y_domain": None,
        }

        self.x_window = x_window or self._default_domain_window["x_window"]
        self.y_window = y_window or self._default_domain_window["y_window"]
        self.x_domain = x_domain
        self.y_domain = y_domain

        self._param_names = self._generate_coeff_names()
        if n_models:
            if model_set_axis is None:
                model_set_axis = 0
            minshape = (1,) * model_set_axis + (n_models,)
        else:
            minshape = ()

        for param_name in self._param_names:
            self._parameters_[param_name] = Parameter(
                param_name, default=np.zeros(minshape)
            )
        super().__init__(
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    @property
    def x_domain(self):
        return self._x_domain

    @x_domain.setter
    def x_domain(self, val):
        self._x_domain = _validate_domain_window(val)

    @property
    def y_domain(self):
        return self._y_domain

    @y_domain.setter
    def y_domain(self, val):
        self._y_domain = _validate_domain_window(val)

    @property
    def x_window(self):
        return self._x_window

    @x_window.setter
    def x_window(self, val):
        self._x_window = _validate_domain_window(val)

    @property
    def y_window(self):
        return self._y_window

    @y_window.setter
    def y_window(self, val):
        self._y_window = _validate_domain_window(val)

    def __repr__(self):
        return self._format_repr(
            [self.x_degree, self.y_degree],
            kwargs={
                "x_domain": self.x_domain,
                "y_domain": self.y_domain,
                "x_window": self.x_window,
                "y_window": self.y_window,
            },
            defaults=self._default_domain_window,
        )

    def __str__(self):
        return self._format_str(
            [
                ("X_Degree", self.x_degree),
                ("Y_Degree", self.y_degree),
                ("X_Domain", self.x_domain),
                ("Y_Domain", self.y_domain),
                ("X_Window", self.x_window),
                ("Y_Window", self.y_window),
            ],
            self._default_domain_window,
        )

    def get_num_coeff(self):
        """
        Determine how many coefficients are needed.

        Returns
        -------
        numc : int
            number of coefficients
        """
        if self.x_degree < 0 or self.y_degree < 0:
            raise ValueError("Degree of polynomial must be positive or null")

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

    def invlex_coeff(self, coeffs):
        invlex_coeffs = []
        xvar = np.arange(self.x_degree + 1)
        yvar = np.arange(self.y_degree + 1)
        for j in yvar:
            for i in xvar:
                name = f"c{i}_{j}"
                coeff = coeffs[self.param_names.index(name)]
                invlex_coeffs.append(coeff)
        return np.array(invlex_coeffs[::-1])

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
            setattr(self, "r" + str(n), 0.0)

        for n in range(1, nalpha):
            k = karr[n - 1].nonzero()[0].max() + 1
            rsum = 0
            for i in range(1, k + 1):
                rsum = rsum + getattr(self, "r" + str(i))
            val = kfunc[k - 1] * (r0 + rsum)
            setattr(self, "r" + str(k), val)
            r0 = _coeff[n]
            for i in range(1, k):
                setattr(self, "r" + str(i), 0.0)
        result = r0
        for i in range(1, nterms + 1 + 3):
            result = result + getattr(self, "r" + str(i))
        return result

    def _generate_coeff_names(self):
        names = []
        for j in range(self.y_degree + 1):
            for i in range(self.x_degree + 1):
                names.append(f"c{i}_{j}")
        return tuple(names)

    def _fcache(self, x, y):
        """
        Computation and store the individual functions.

        To be implemented by subclasses"
        """
        raise NotImplementedError("Subclasses should implement this")

    def evaluate(self, x, y, *coeffs):
        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)
        invcoeff = self.invlex_coeff(coeffs)
        return self.imhorner(x, y, invcoeff)

    def prepare_inputs(self, x, y, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, y, **kwargs)

        x, y = inputs

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        return (x, y), broadcasted_shapes


class Chebyshev1D(_PolyDomainWindow1D):
    r"""
    Univariate Chebyshev series.

    It is defined as:

    .. math::

        P(x) = \sum_{i=0}^{i=n}C_{i} * T_{i}(x)

    where ``T_i(x)`` is the corresponding Chebyshev polynomial of the 1st kind.

    For explanation of ```domain``, and ``window`` see
    :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : tuple or None, optional
    window : tuple or None, optional
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window.
    **params : dict
        keyword : value pairs, representing parameter_name: value

    Notes
    -----
    This model does not support the use of units/quantities, because each term
    in the sum of Chebyshev polynomials is a polynomial in x - since the
    coefficients within each Chebyshev polynomial are fixed, we can't use
    quantities for x since the units would not be compatible. For example, the
    third Chebyshev polynomial (T2) is 2x^2-1, but if x was specified with
    units, 2x^2 and -1 would have incompatible units.
    """

    n_inputs = 1
    n_outputs = 1

    _separable = True

    def __init__(
        self,
        degree,
        domain=None,
        window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree,
            domain=domain,
            window=window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
        v = np.empty((self.degree + 1,) + x.shape, dtype=x.dtype)
        v[0] = 1
        if self.degree > 0:
            x2 = 2 * x
            v[1] = x
            for i in range(2, self.degree + 1):
                v[i] = v[i - 1] * x2 - v[i - 2]
        return np.rollaxis(v, 0, v.ndim)

    def prepare_inputs(self, x, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, **kwargs)

        x = inputs[0]

        return (x,), broadcasted_shapes

    def evaluate(self, x, *coeffs):
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)
        return self.clenshaw(x, coeffs)

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


class Hermite1D(_PolyDomainWindow1D):
    r"""
    Univariate Hermite series.

    It is defined as:

    .. math::

        P(x) = \sum_{i=0}^{i=n}C_{i} * H_{i}(x)

    where ``H_i(x)`` is the corresponding Hermite polynomial ("Physicist's kind").

    For explanation of ``domain``, and ``window`` see
    :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : tuple or None, optional
    window : tuple or None, optional
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    **params : dict
        keyword : value pairs, representing parameter_name: value

    Notes
    -----
    This model does not support the use of units/quantities, because each term
    in the sum of Hermite polynomials is a polynomial in x - since the
    coefficients within each Hermite polynomial are fixed, we can't use
    quantities for x since the units would not be compatible. For example, the
    third Hermite polynomial (H2) is 4x^2-2, but if x was specified with units,
    4x^2 and -2 would have incompatible units.
    """

    n_inputs = 1
    n_outputs = 1

    _separable = True

    def __init__(
        self,
        degree,
        domain=None,
        window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree,
            domain,
            window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
        v = np.empty((self.degree + 1,) + x.shape, dtype=x.dtype)
        v[0] = 1
        if self.degree > 0:
            x2 = 2 * x
            v[1] = 2 * x
            for i in range(2, self.degree + 1):
                v[i] = x2 * v[i - 1] - 2 * (i - 1) * v[i - 2]
        return np.rollaxis(v, 0, v.ndim)

    def prepare_inputs(self, x, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, **kwargs)

        x = inputs[0]

        return (x,), broadcasted_shapes

    def evaluate(self, x, *coeffs):
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)
        return self.clenshaw(x, coeffs)

    @staticmethod
    def clenshaw(x, coeffs):
        x2 = x * 2
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
                temp = c0
                nd = nd - 1
                c0 = coeffs[-i] - c1 * (2 * (nd - 1))
                c1 = temp + c1 * x2
        return c0 + c1 * x2


class Hermite2D(OrthoPolynomialBase):
    r"""
    Bivariate Hermite series.

    It is defined as

    .. math:: P_{nm}(x,y) = \sum_{n,m=0}^{n=d,m=d}C_{nm} H_n(x) H_m(y)

    where ``H_n(x)`` and ``H_m(y)`` are Hermite polynomials.

    For explanation of ``x_domain``, ``y_domain``, ``x_window`` and ``y_window``
    see :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : tuple or None, optional
        domain of the x independent variable
    y_domain : tuple or None, optional
        domain of the y independent variable
    x_window : tuple or None, optional
        range of the x independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    y_window : tuple or None, optional
        range of the y independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    **params : dict
        keyword: value pairs, representing parameter_name: value

    Notes
    -----
    This model does not support the use of units/quantities, because each term
    in the sum of Hermite polynomials is a polynomial in x and/or y - since the
    coefficients within each Hermite polynomial are fixed, we can't use
    quantities for x and/or y since the units would not be compatible. For
    example, the third Hermite polynomial (H2) is 4x^2-2, but if x was
    specified with units, 4x^2 and -2 would have incompatible units.
    """

    _separable = False

    def __init__(
        self,
        x_degree,
        y_degree,
        x_domain=None,
        x_window=None,
        y_domain=None,
        y_window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            x_degree,
            y_degree,
            x_domain=x_domain,
            y_domain=y_domain,
            x_window=x_window,
            y_window=y_window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    def _fcache(self, x, y):
        """
        Calculate the individual Hermite functions once and store them in a
        dictionary to be reused.
        """
        x_terms = self.x_degree + 1
        y_terms = self.y_degree + 1
        kfunc = {}
        kfunc[0] = np.ones(x.shape)
        kfunc[1] = 2 * x.copy()
        kfunc[x_terms] = np.ones(y.shape)
        kfunc[x_terms + 1] = 2 * y.copy()
        for n in range(2, x_terms):
            kfunc[n] = 2 * x * kfunc[n - 1] - 2 * (n - 1) * kfunc[n - 2]
        for n in range(x_terms + 2, x_terms + y_terms):
            kfunc[n] = 2 * y * kfunc[n - 1] - 2 * (n - 1) * kfunc[n - 2]
        return kfunc

    def fit_deriv(self, x, y, *params):
        """
        Derivatives with respect to the coefficients.

        This is an array with Hermite polynomials:

        .. math::

            H_{x_0}H_{y_0}, H_{x_1}H_{y_0}...H_{x_n}H_{y_0}...H_{x_n}H_{y_m}

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")

        x = x.ravel()
        y = y.ravel()
        x_deriv = self._hermderiv1d(x, self.x_degree + 1).T
        y_deriv = self._hermderiv1d(y, self.y_degree + 1).T

        ij = []
        for i in range(self.y_degree + 1):
            for j in range(self.x_degree + 1):
                ij.append(x_deriv[j] * y_deriv[i])

        v = np.array(ij)
        return v.T

    def _hermderiv1d(self, x, deg):
        """
        Derivative of 1D Hermite series.
        """
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
        d = np.empty((deg + 1, len(x)), dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0:
            x2 = 2 * x
            d[1] = x2
            for i in range(2, deg + 1):
                d[i] = x2 * d[i - 1] - 2 * (i - 1) * d[i - 2]
        return np.rollaxis(d, 0, d.ndim)


class Legendre1D(_PolyDomainWindow1D):
    r"""
    Univariate Legendre series.

    It is defined as:

    .. math::

        P(x) = \sum_{i=0}^{i=n}C_{i} * L_{i}(x)

    where ``L_i(x)`` is the corresponding Legendre polynomial.

    For explanation of ``domain``, and ``window`` see
    :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : tuple or None, optional
    window : tuple or None, optional
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    **params : dict
        keyword: value pairs, representing parameter_name: value


    Notes
    -----
    This model does not support the use of units/quantities, because each term
    in the sum of Legendre polynomials is a polynomial in x - since the
    coefficients within each Legendre polynomial are fixed, we can't use
    quantities for x since the units would not be compatible. For example, the
    third Legendre polynomial (P2) is 1.5x^2-0.5, but if x was specified with
    units, 1.5x^2 and -0.5 would have incompatible units.
    """

    n_inputs = 1
    n_outputs = 1

    _separable = True

    def __init__(
        self,
        degree,
        domain=None,
        window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree,
            domain,
            window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    def prepare_inputs(self, x, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, **kwargs)

        x = inputs[0]

        return (x,), broadcasted_shapes

    def evaluate(self, x, *coeffs):
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)
        return self.clenshaw(x, coeffs)

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
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


class Polynomial1D(_PolyDomainWindow1D):
    r"""
    1D Polynomial model.

    It is defined as:

    .. math::

        P = \sum_{i=0}^{i=n}C_{i} * x^{i}

    For explanation of ``domain``, and ``window`` see
    :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    degree : int
        degree of the series
    domain : tuple or None, optional
        If None, it is set to (-1, 1)
    window : tuple or None, optional
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    **params : dict
        keyword: value pairs, representing parameter_name: value

    """

    n_inputs = 1
    n_outputs = 1

    _separable = True

    def __init__(
        self,
        degree,
        domain=None,
        window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree,
            domain,
            window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

        # Set domain separately because it's different from
        # the orthogonal polynomials.
        self._default_domain_window = {
            "domain": (-1, 1),
            "window": (-1, 1),
        }
        self.domain = domain or self._default_domain_window["domain"]
        self.window = window or self._default_domain_window["window"]

    def prepare_inputs(self, x, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, **kwargs)

        x = inputs[0]
        return (x,), broadcasted_shapes

    def evaluate(self, x, *coeffs):
        if self.domain is not None:
            x = poly_map_domain(x, self.domain, self.window)
        return self.horner(x, coeffs)

    def fit_deriv(self, x, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        v = np.empty((self.degree + 1,) + x.shape, dtype=float)
        v[0] = 1
        if self.degree > 0:
            v[1] = x
            for i in range(2, self.degree + 1):
                v[i] = v[i - 1] * x
        return np.rollaxis(v, 0, v.ndim)

    @staticmethod
    def horner(x, coeffs):
        if len(coeffs) == 1:
            c0 = coeffs[-1] * np.ones_like(x, subok=False)
        else:
            c0 = coeffs[-1]
            for i in range(2, len(coeffs) + 1):
                c0 = coeffs[-i] + c0 * x
        return c0

    @property
    def input_units(self):
        if self.degree == 0 or self.c1.input_unit is None:
            return None
        else:
            return {self.inputs[0]: self.c0.input_unit / self.c1.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        mapping = {}
        for i in range(self.degree + 1):
            par = getattr(self, f"c{i}")
            mapping[par.name] = (
                outputs_unit[self.outputs[0]] / inputs_unit[self.inputs[0]] ** i
            )
        return mapping


class Polynomial2D(PolynomialModel):
    r"""
    2D Polynomial  model.

    Represents a general polynomial of degree n:

    .. math::

        P(x,y) = c_{00} + c_{10}x + ...+ c_{n0}x^n + c_{01}y + ...+ c_{0n}y^n
        + c_{11}xy + c_{12}xy^2 + ... + c_{1(n-1)}xy^{n-1}+ ... + c_{(n-1)1}x^{n-1}y

    For explanation of ``x_domain``, ``y_domain``, ``x_window`` and ``y_window``
    see :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    degree : int
        Polynomial degree: largest sum of exponents (:math:`i + j`) of
        variables in each monomial term of the form :math:`x^i y^j`. The
        number of terms in a 2D polynomial of degree ``n`` is given by binomial
        coefficient :math:`C(n + 2, 2) = (n + 2)! / (2!\,n!) = (n + 1)(n + 2) / 2`.
    x_domain : tuple or None, optional
        domain of the x independent variable
        If None, it is set to (-1, 1)
    y_domain : tuple or None, optional
        domain of the y independent variable
        If None, it is set to (-1, 1)
    x_window : tuple or None, optional
        range of the x independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the x_domain to x_window
    y_window : tuple or None, optional
        range of the y independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the y_domain to y_window
    **params : dict
        keyword: value pairs, representing parameter_name: value
    """

    n_inputs = 2
    n_outputs = 1

    _separable = False

    def __init__(
        self,
        degree,
        x_domain=None,
        y_domain=None,
        x_window=None,
        y_window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            degree,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

        self._default_domain_window = {
            "x_domain": (-1, 1),
            "y_domain": (-1, 1),
            "x_window": (-1, 1),
            "y_window": (-1, 1),
        }

        self.x_domain = x_domain or self._default_domain_window["x_domain"]
        self.y_domain = y_domain or self._default_domain_window["y_domain"]
        self.x_window = x_window or self._default_domain_window["x_window"]
        self.y_window = y_window or self._default_domain_window["y_window"]

    def prepare_inputs(self, x, y, **kwargs):
        inputs, broadcasted_shapes = super().prepare_inputs(x, y, **kwargs)

        x, y = inputs
        return (x, y), broadcasted_shapes

    def evaluate(self, x, y, *coeffs):
        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)
        invcoeff = self.invlex_coeff(coeffs)
        result = self.multivariate_horner(x, y, invcoeff)

        # Special case for degree==0 to ensure that the shape of the output is
        # still as expected by the broadcasting rules, even though the x and y
        # inputs are not used in the evaluation
        if self.degree == 0:
            output_shape = check_broadcast(np.shape(coeffs[0]), x.shape)
            if output_shape:
                new_result = np.empty(output_shape)
                new_result[:] = result
                result = new_result

        return result

    def __repr__(self):
        return self._format_repr(
            [self.degree],
            kwargs={
                "x_domain": self.x_domain,
                "y_domain": self.y_domain,
                "x_window": self.x_window,
                "y_window": self.y_window,
            },
            defaults=self._default_domain_window,
        )

    def __str__(self):
        return self._format_str(
            [
                ("Degree", self.degree),
                ("X_Domain", self.x_domain),
                ("Y_Domain", self.y_domain),
                ("X_Window", self.x_window),
                ("Y_Window", self.y_window),
            ],
            self._default_domain_window,
        )

    def fit_deriv(self, x, y, *params):
        """
        Computes the Vandermonde matrix.

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.ndim == 2:
            x = x.ravel()
        if y.ndim == 2:
            y = y.ravel()
        if x.size != y.size:
            raise ValueError("Expected x and y to be of equal size")

        designx = x[:, None] ** np.arange(self.degree + 1)
        designy = y[:, None] ** np.arange(1, self.degree + 1)

        designmixed = []
        for i in range(1, self.degree):
            for j in range(1, self.degree):
                if i + j <= self.degree:
                    designmixed.append((x**i) * (y**j))
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
                    name = f"c{j}_{i}"
                    coeff = coeffs[self.param_names.index(name)]
                    invlex_coeffs.append(coeff)
        return invlex_coeffs[::-1]

    def multivariate_horner(self, x, y, coeffs):
        """
        Multivariate Horner's scheme.

        Parameters
        ----------
        x, y : array
        coeffs : array
            Coefficients in inverse lexical order.
        """
        alpha = self._invlex()
        r0 = coeffs[0]
        r1 = r0 * 0.0
        r2 = r0 * 0.0
        karr = np.diff(alpha, axis=0)

        for n in range(len(karr)):
            if karr[n, 1] != 0:
                r2 = y * (r0 + r1 + r2)
                r1 = np.zeros_like(coeffs[0], subok=False)
            else:
                r1 = x * (r0 + r1)
            r0 = coeffs[n + 1]
        return r0 + r1 + r2

    @property
    def input_units(self):
        if self.degree == 0 or (
            self.c1_0.input_unit is None and self.c0_1.input_unit is None
        ):
            return None
        return {
            self.inputs[0]: self.c0_0.input_unit / self.c1_0.input_unit,
            self.inputs[1]: self.c0_0.input_unit / self.c0_1.input_unit,
        }

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        mapping = {}
        for i in range(self.degree + 1):
            for j in range(self.degree + 1):
                if i + j > 2:
                    continue
                par = getattr(self, f"c{i}_{j}")
                mapping[par.name] = (
                    outputs_unit[self.outputs[0]]
                    / inputs_unit[self.inputs[0]] ** i
                    / inputs_unit[self.inputs[1]] ** j
                )
        return mapping

    @property
    def x_domain(self):
        return self._x_domain

    @x_domain.setter
    def x_domain(self, val):
        self._x_domain = _validate_domain_window(val)

    @property
    def y_domain(self):
        return self._y_domain

    @y_domain.setter
    def y_domain(self, val):
        self._y_domain = _validate_domain_window(val)

    @property
    def x_window(self):
        return self._x_window

    @x_window.setter
    def x_window(self, val):
        self._x_window = _validate_domain_window(val)

    @property
    def y_window(self):
        return self._y_window

    @y_window.setter
    def y_window(self, val):
        self._y_window = _validate_domain_window(val)


class Chebyshev2D(OrthoPolynomialBase):
    r"""
    Bivariate Chebyshev series..

    It is defined as

    .. math:: P_{nm}(x,y) = \sum_{n,m=0}^{n=d,m=d}C_{nm}  T_n(x ) T_m(y)

    where ``T_n(x)`` and ``T_m(y)`` are Chebyshev polynomials of the first kind.

    For explanation of ``x_domain``, ``y_domain``, ``x_window`` and ``y_window``
    see :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : tuple or None, optional
        domain of the x independent variable
    y_domain : tuple or None, optional
        domain of the y independent variable
    x_window : tuple or None, optional
        range of the x independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    y_window : tuple or None, optional
        range of the y independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window

    **params : dict
        keyword: value pairs, representing parameter_name: value

    Notes
    -----
    This model does not support the use of units/quantities, because each term
    in the sum of Chebyshev polynomials is a polynomial in x and/or y - since
    the coefficients within each Chebyshev polynomial are fixed, we can't use
    quantities for x and/or y since the units would not be compatible. For
    example, the third Chebyshev polynomial (T2) is 2x^2-1, but if x was
    specified with units, 2x^2 and -1 would have incompatible units.
    """

    _separable = False

    def __init__(
        self,
        x_degree,
        y_degree,
        x_domain=None,
        x_window=None,
        y_domain=None,
        y_window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            x_degree,
            y_degree,
            x_domain=x_domain,
            y_domain=y_domain,
            x_window=x_window,
            y_window=y_window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

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
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")

        x = x.ravel()
        y = y.ravel()
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
        Derivative of 1D Chebyshev series.
        """
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
        d = np.empty((deg + 1, len(x)), dtype=x.dtype)
        d[0] = x * 0 + 1
        if deg > 0:
            x2 = 2 * x
            d[1] = x
            for i in range(2, deg + 1):
                d[i] = d[i - 1] * x2 - d[i - 2]
        return np.rollaxis(d, 0, d.ndim)


class Legendre2D(OrthoPolynomialBase):
    r"""
    Bivariate Legendre series.

    Defined as:

    .. math:: P_{n_m}(x,y) = \sum_{n,m=0}^{n=d,m=d}C_{nm}  L_n(x ) L_m(y)

    where ``L_n(x)`` and ``L_m(y)`` are Legendre polynomials.

    For explanation of ``x_domain``, ``y_domain``, ``x_window`` and ``y_window``
    see :ref:`Notes regarding usage of domain and window <domain-window-note>`.

    Parameters
    ----------
    x_degree : int
        degree in x
    y_degree : int
        degree in y
    x_domain : tuple or None, optional
        domain of the x independent variable
    y_domain : tuple or None, optional
        domain of the y independent variable
    x_window : tuple or None, optional
        range of the x independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    y_window : tuple or None, optional
        range of the y independent variable
        If None, it is set to (-1, 1)
        Fitters will remap the domain to this window
    **params : dict
        keyword: value pairs, representing parameter_name: value

    Notes
    -----
    Model formula:

    .. math::

        P(x) = \sum_{i=0}^{i=n}C_{i} * L_{i}(x)

    where ``L_{i}`` is the corresponding Legendre polynomial.

    This model does not support the use of units/quantities, because each term
    in the sum of Legendre polynomials is a polynomial in x - since the
    coefficients within each Legendre polynomial are fixed, we can't use
    quantities for x since the units would not be compatible. For example, the
    third Legendre polynomial (P2) is 1.5x^2-0.5, but if x was specified with
    units, 1.5x^2 and -0.5 would have incompatible units.
    """

    _separable = False

    def __init__(
        self,
        x_degree,
        y_degree,
        x_domain=None,
        x_window=None,
        y_domain=None,
        y_window=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        super().__init__(
            x_degree,
            y_degree,
            x_domain=x_domain,
            y_domain=y_domain,
            x_window=x_window,
            y_window=y_window,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

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
            kfunc[n] = (
                (2 * (n - 1) + 1) * x * kfunc[n - 1] - (n - 1) * kfunc[n - 2]
            ) / n
        for n in range(2, y_terms):
            kfunc[n + x_terms] = (
                (2 * (n - 1) + 1) * y * kfunc[n + x_terms - 1]
                - (n - 1) * kfunc[n + x_terms - 2]
            ) / (n)
        return kfunc

    def fit_deriv(self, x, y, *params):
        """Derivatives with respect to the coefficients.

        This is an array with Legendre polynomials:

        Lx0Ly0  Lx1Ly0...LxnLy0...LxnLym

        Parameters
        ----------
        x : ndarray
            input
        y : ndarray
            input
        *params
            throw-away parameter list returned by non-linear fitters

        Returns
        -------
        result : ndarray
            The Vandermonde matrix
        """
        if x.shape != y.shape:
            raise ValueError("x and y must have the same shape")
        x = x.ravel()
        y = y.ravel()
        x_deriv = self._legendderiv1d(x, self.x_degree + 1).T
        y_deriv = self._legendderiv1d(y, self.y_degree + 1).T

        ij = []
        for i in range(self.y_degree + 1):
            for j in range(self.x_degree + 1):
                ij.append(x_deriv[j] * y_deriv[i])

        v = np.array(ij)
        return v.T

    def _legendderiv1d(self, x, deg):
        """Derivative of 1D Legendre polynomial."""
        x = np.array(x, dtype=float, copy=COPY_IF_NEEDED, ndmin=1)
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

    n_inputs = 2
    n_outputs = 1

    _separable = False

    def __init__(
        self,
        order,
        coeff_prefix,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
        **params,
    ):
        self.order = order
        self.coeff_prefix = coeff_prefix
        self._param_names = self._generate_coeff_names(coeff_prefix)

        if n_models:
            if model_set_axis is None:
                model_set_axis = 0
            minshape = (1,) * model_set_axis + (n_models,)
        else:
            minshape = ()
        for param_name in self._param_names:
            self._parameters_[param_name] = Parameter(
                param_name, default=np.zeros(minshape)
            )
        super().__init__(
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
            **params,
        )

    def __repr__(self):
        return self._format_repr(args=[self.order, self.coeff_prefix])

    def __str__(self):
        return self._format_str(
            [("Order", self.order), ("Coeff. Prefix", self.coeff_prefix)]
        )

    def evaluate(self, x, y, *coeffs):
        # TODO: Rewrite this so that it uses a simpler method of determining
        # the matrix based on the number of given coefficients.
        mcoef = self._coeff_matrix(self.coeff_prefix, coeffs)
        return self._eval_sip(x, y, mcoef)

    def get_num_coeff(self, ndim):
        """
        Return the number of coefficients in one param set.
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
            names.append(f"{coeff_prefix}_{i}_{0}")
        for i in range(2, self.order + 1):
            names.append(f"{coeff_prefix}_{0}_{i}")
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    names.append(f"{coeff_prefix}_{i}_{j}")
        return tuple(names)

    def _coeff_matrix(self, coeff_prefix, coeffs):
        mat = np.zeros((self.order + 1, self.order + 1))
        for i in range(2, self.order + 1):
            attr = f"{coeff_prefix}_{i}_{0}"
            mat[i, 0] = coeffs[self.param_names.index(attr)][0]
        for i in range(2, self.order + 1):
            attr = f"{coeff_prefix}_{0}_{i}"
            mat[0, i] = coeffs[self.param_names.index(attr)][0]
        for i in range(1, self.order):
            for j in range(1, self.order):
                if i + j < self.order + 1:
                    attr = f"{coeff_prefix}_{i}_{j}"
                    mat[i, j] = coeffs[self.param_names.index(attr)][0]
        return mat

    def _eval_sip(self, x, y, coef):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if self.coeff_prefix == "A":
            result = np.zeros(x.shape)
        else:
            result = np.zeros(y.shape)

        for i in range(coef.shape[0]):
            for j in range(coef.shape[1]):
                if 1 < i + j < self.order + 1:
                    result = result + coef[i, j] * x**i * y**j
        return result


class SIP(Model):
    """
    Simple Imaging Polynomial (SIP) model.

    The SIP convention is used to represent distortions in FITS image headers.
    See [1]_ for a description of the SIP convention.

    Parameters
    ----------
    crpix : list or (2,) ndarray
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

    References
    ----------
    .. [1] `David Shupe, et al, ADASS, ASP Conference Series, Vol. 347, 2005
        <https://ui.adsabs.harvard.edu/abs/2005ASPC..347..491S>`_
    """

    n_inputs = 2
    n_outputs = 2

    _separable = False

    def __init__(
        self,
        crpix,
        a_order,
        b_order,
        a_coeff={},
        b_coeff={},
        ap_order=None,
        bp_order=None,
        ap_coeff={},
        bp_coeff={},
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
    ):
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
        self.sip1d_a = _SIP1D(
            a_order,
            coeff_prefix="A",
            n_models=n_models,
            model_set_axis=model_set_axis,
            **a_coeff,
        )
        self.sip1d_b = _SIP1D(
            b_order,
            coeff_prefix="B",
            n_models=n_models,
            model_set_axis=model_set_axis,
            **b_coeff,
        )
        super().__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name, meta=meta
        )
        self._inputs = ("u", "v")
        self._outputs = ("x", "y")

    def __repr__(self):
        return (
            f"<{self.__class__.__name__}"
            f"({[self.shift_a, self.shift_b, self.sip1d_a, self.sip1d_b]!r})>"
        )

    def __str__(self):
        parts = [f"Model: {self.__class__.__name__}"]
        for model in [self.shift_a, self.shift_b, self.sip1d_a, self.sip1d_b]:
            parts.append(indent(str(model), 4 * " "))
            parts.append("")

        return "\n".join(parts)

    @property
    def inverse(self):
        if self._ap_order is not None and self._bp_order is not None:
            return InverseSIP(
                self._ap_order, self._bp_order, self._ap_coeff, self._bp_coeff
            )
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
    Inverse Simple Imaging Polynomial.

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

    """

    n_inputs = 2
    n_outputs = 2

    _separable = False

    def __init__(
        self,
        ap_order,
        bp_order,
        ap_coeff={},
        bp_coeff={},
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
    ):
        self._ap_order = ap_order
        self._bp_order = bp_order
        self._ap_coeff = ap_coeff
        self._bp_coeff = bp_coeff

        # define the 0th term in order to use Polynomial2D
        ap_coeff.setdefault("AP_0_0", 0)
        bp_coeff.setdefault("BP_0_0", 0)

        ap_coeff_params = {k.replace("AP_", "c"): v for k, v in ap_coeff.items()}
        bp_coeff_params = {k.replace("BP_", "c"): v for k, v in bp_coeff.items()}

        self.sip1d_ap = Polynomial2D(
            degree=ap_order, model_set_axis=model_set_axis, **ap_coeff_params
        )
        self.sip1d_bp = Polynomial2D(
            degree=bp_order, model_set_axis=model_set_axis, **bp_coeff_params
        )
        super().__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name, meta=meta
        )

    def __repr__(self):
        return f"<{self.__class__.__name__}({[self.sip1d_ap, self.sip1d_bp]!r})>"

    def __str__(self):
        parts = [f"Model: {self.__class__.__name__}"]
        for model in [self.sip1d_ap, self.sip1d_bp]:
            parts.append(indent(str(model), 4 * " "))
            parts.append("")

        return "\n".join(parts)

    def evaluate(self, x, y):
        x1 = self.sip1d_ap.evaluate(x, y, *self.sip1d_ap.param_sets)
        y1 = self.sip1d_bp.evaluate(x, y, *self.sip1d_bp.param_sets)
        return x1, y1
