# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains predefined polynomial models.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import abc
import itertools

import numpy as np

from .core import ModelDefinitionError, FittableModel, Model, _ModelMeta
from .functional_models import Shift
from .parameters import Parameter
from .utils import poly_map_domain, comb, check_broadcast, polynomial_exponents
from ..extern import six
from ..extern.six.moves import copyreg
from ..utils import lazyproperty, classproperty, indent, deprecated
from ..utils.compat import funcsigs
from ..utils.compat.functools import lru_cache


__all__ = [
    'Chebyshev1D', 'Chebyshev2D', 'InverseSIP',
    'Legendre1D', 'Legendre2D', 'Polynomial1D',
    'Polynomial2D', 'SIP', 'PolynomialBase', 'OrthoPolynomialBase'
]


class _PolynomialModelMeta(_ModelMeta):
    def __init__(cls, name, bases, members):
        # Skip the _get_subclass_for_degree setup if this is already a concrete
        # subclass
        if members.get('_is_dynamic'):
            return

        # This is so that each subclass gets its own lru_cache for per-degree
        # subclasses
        cache_wrapper = lru_cache(maxsize=cls._degree_subclass_cache_maxsize)

        # Use __func__--otherwise lru_cache will wrap the existing class
        # method, which makes re-wrapping it as a classmethod on the
        # metaclass instance tricky...
        if hasattr(cls._get_subclass_for_degree, '__wrapper__'):
            getter = cls._get_subclass_for_degree.__wrapper__.__func__
        else:
            getter = cls._get_subclass_for_degree.__func__

        # This goes directly into the class __dict__, so we have to make it
        # a classmethod object (otherwise lookup on this method won't go
        # through the `type` class's method-getter).
        cls._get_subclass_for_degree = classmethod(cache_wrapper(getter))

        # Add the degree attributes as classproperties (making them read-only)
        degree_attrs = cls._get_degree_attrs()

        for name, doc in cls._get_degree_attrs():

            def _degree_getter(cls, name=name):
                return getattr(cls, '_' + name)

            setattr(cls, '_' + name, None)
            setattr(cls, name, classproperty(_degree_getter, doc=doc))

    def _get_subclass_for_degree(cls, *degrees):
        # Creates a subclass of the given Polynomial class that has its degree
        # fixed (and has the appropriate coefficient Parameters to go with it)
        # If only one degree is given this is the maximum degree for any term
        # in the polynomial (regardless of number of variables).
        # If multiple degrees are given these correspond with each independent
        # variable individually (as in the orthogonal polynomial classes)

        # 'degree' attributes are implemented as classproperties so that we
        # can given them docstrings too

        # Internal '_degrees' attribute mostly just for book-keeping
        members = {
            '_degrees': degrees,
            '_is_dynamic': True  # See _ModelMeta.__reduce__
        }

        # See the degree attribute values (on their internal attributes)
        for attr, degree in zip(cls._get_degree_attrs(), degrees):
            members['_' + attr[0]] = degree

        for exponents in cls._term_powers(*degrees):
            coeff_name = cls._coefficient_name(*exponents)
            parameter = cls._make_coefficient_parameter(coeff_name,
                                                        *exponents)
            members[coeff_name] = parameter

        classname = 'Degree_{0}_{1}'.format('_'.join(str(x) for x in degrees),
                                            cls.__name__)

        return type(str(classname), (cls,), members)

    @classmethod
    def _get_init_signature(mcls, cls, parameters):
        # Overrides _ModelMeta version of this
        degree_attrs = tuple(attr[0] for attr in cls._get_degree_attrs())
        arg_types = super(_PolynomialModelMeta, mcls)._get_init_signature(
                cls, parameters)
        args = arg_types[0]
        kwargs = arg_types[1]

        # args[0] should be 'self'
        args = (args[0],) + degree_attrs + args[1:]

        # add domain/window keyword args--get the default values from the
        # base class
        base_cls = cls.__bases__[0]  # TODO: This could by flaky
        base_init_sig = funcsigs.signature(base_cls.__init__)
        if cls.n_inputs == 1:
            arg_prefixes = ('',)
        else:
            arg_prefixes = (in_ + '_' for in_ in cls.inputs[::-1])

        n_params = len(cls.param_names)

        for arg_prefix in arg_prefixes:
            for arg_base in ('window', 'domain'):
                # Just 'domain' for 1D; x_domain, y_domain, etc.
                # for >= 2D
                arg_name = arg_prefix + arg_base
                param = base_init_sig.parameters.get(arg_name)

                if param is None:
                    # Maybe this model doesn't take a window/domain argument?
                    continue

                # Insert into kwargs list after the parameter kwargs
                kwargs.insert(n_params, (arg_name, param.default))

        return (args,) + arg_types[1:]


@six.add_metaclass(_PolynomialModelMeta)
class PolynomialBase(FittableModel):
    """
    Base class for all polynomial-like models with an arbitrary number of
    parameters in the form of coefficients.

    In this case Parameter instances are returned through the class's
    ``__getattr__`` rather than through class descriptors.
    """

    linear = True
    col_fit_deriv = False

    default_coefficient = 0.0
    """Default value given to all coefficients."""

    _degrees = None
    """
    This internal attribute is set by
    _PolynomialModelMeta._get_subclass_for_degree, and just stores the degree
    of the polynomial (or degrees for multiple-indexed polynomials).  A None
    value for this attribute indicates an "abstract" polynomial in the sense
    that it does not have a specific degree, and can't be instantiated by
    itself.
    """

    @abc.abstractproperty
    def coefficient_template(self):
        """
        Template for coefficient parameter names.  Template should use
        positional template arguments corresponding to exponent on each of
        the input variables.

        For example for a 2D polynomial one might use::

            coefficient_template = 'c{0}_{1}'

        This will generate coefficient names like ``c1_2`` corresponding
        to the coefficient on the :math:`xy^2` term.
        """

    _degree_subclass_cache_maxsize = 16
    """
    When creating a Polynomial instance of a given degree, a subclass of that
    Polynomial class is created representing a polynomial of that degree.
    These subclasses are cached in order to speed up instance creation, but
    there is a cap on the number of classes that are cached.

    If, in the rare case that cap needs to be changed (say, we're creating
    lots of polynomials of a wide range of degrees) that cap can be lifted
    by creating a subclass of the polynomial type in question with an
    increased value for this attribute.

    The default value of 16 is probably reasonable, but chosen somewhat
    arbitrarily (though the cache gets best performance if the size is a power
    of two).
    """

    # TODO: The argument parsing in __new__ is awkward--might be good to have a
    # better way of specifying how to look for the 'degree' arguments

    def __new__(cls, *args, **kwargs):
        if cls._degrees is None:
            degrees = []
            for idx, arg in enumerate(cls._get_degree_attrs()):
                if arg[0] in kwargs:
                    degrees.append(kwargs[arg[0]])
                elif args and idx <= len(args) - 1:
                    # The first N positional arguments should be the N
                    # different degree arguments for n_degrees = N
                    degrees.append(args[idx])
                else:
                    raise TypeError("missing required argument '{0}'".format(
                        arg[0]))

            new_cls = cls._get_subclass_for_degree(*degrees)
        else:
            # Just pass through the existing cls, which already has some
            # specific degree--this can happen in the case of copying /
            # unpickling an existing polynomial model with fixed degree
            new_cls = cls
        return super(PolynomialBase, cls).__new__(new_cls)

    # TODO: Make this 'lazy'
    @classproperty
    def n_coefficients(cls):
        """
        The number of coefficients in this polynomial (equivalently the
        number of terms in the polynomial).

        In most cases this will be identical to ``len(model.param_names)``
        since each coefficient is a fittable parameter.  However, some classes
        may define additional parameters that are not coefficients.
        """

        return len(list(cls._term_powers(*cls._degrees)))

    # The following classmethods may be overridden by subclasses to determine
    # the number of terms in the polynomial and indexing / naming of
    # coefficients

    @classmethod
    def _get_degree_attrs(cls):
        """
        Returns the names of class attributes representing the degree of the
        polynomial (or per-variable degrees).  In principle these can be any
        attributes that define a polynomial of the given type (i.e. results in
        a new subclass being created when the polynomial is instantiated).

        Also returns docstrings for those attributes as ``(name, doc)`` tuples.
        """

        return [('degree', "Degree of the polynomial.")]

    @classmethod
    def _term_powers(cls, *degrees):
        """
        Returns an iterator over exponents for terms in the polynomial of given
        degrees (if one degree is given it is the maximum degree of terms in
        polynomial; if multiple degrees are given it is the maximum degree in
        each of the polynomial's variables).

        Here the term "exponents" may be interpreted as any set of indices
        used to index terms in a generalized polynomial.
        """

        return polynomial_exponents(cls.n_inputs, degrees[0])

    @classmethod
    def _coefficient_name(cls, *powers):
        """
        Generates the parameter name for the coefficient to the term with
        the given exponents (in order of the model's inputs).  By default
        this simply uses the `~ParameterBase.coefficient_template` class
        attribute as a template, but this may be overridden by subclasses
        in order to take complete control.

        An alternative implementation may choose to ignore
        ``coefficient_template`` entirely if necessary.
        """

        return cls.coefficient_template.format(*powers)

    @classmethod
    def _make_coefficient_parameter(cls, name, *powers):
        """
        Creates a `Parameter` instance to serve as the unbound parameter
        descriptor for the given coefficient, given the coefficient's name,
        and the exponents on the independent variables identifying which term
        this is the coefficient of.

        `PolynomialBase` subclasses may override this to control creation of
        `Parameter` descriptors (for example to add default constraints, etc.)

        The default just returns a `Parameter` with the given name and the
        default determined by the `~PolynomialBase.default_coefficient` class
        attribute.
        """

        return Parameter(name=name, default=cls.default_coefficient)

    @deprecated('1.1.0', alternative='model.n_coefficients')
    def get_num_coeff(self):
        return self.n_coefficients


class Polynomial1D(PolynomialBase):
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
    coefficient_template = 'c{0}'

    def __init__(self, degree, domain=[-1, 1], window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        self.domain = domain
        self.window = window
        super(Polynomial1D, self).__init__(n_models=n_models,
                                           model_set_axis=model_set_axis,
                                           name=name, meta=meta, **params)

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


class Chebyshev1D(Polynomial1D):
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

    def __init__(self, degree, domain=None, window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Chebyshev1D, self).__init__(
                degree, domain=domain, window=window, n_models=n_models,
                model_set_axis=model_set_axis, name=name, meta=meta, **params)

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


class Legendre1D(Polynomial1D):
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

    def __init__(self, degree, domain=None, window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Legendre1D, self).__init__(
                degree, domain=domain, window=window, n_models=n_models,
                model_set_axis=model_set_axis, name=name, meta=meta, **params)

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


class Polynomial2D(PolynomialBase):
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
    coefficient_template = 'c{0}_{1}'

    def __init__(self, degree, x_domain=[-1, 1], y_domain=[-1, 1],
                 x_window=[-1, 1], y_window=[-1, 1], n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(Polynomial2D, self).__init__(
                n_models=n_models, model_set_axis=model_set_axis,
                name=name, meta=meta, **params)
        # TODO: Perhaps some of these other parameters should be properties?
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
        return (x, y), format_info

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
        return invlex_coeffs[::-1]

    def multivariate_horner(self, x, y, coeffs):
        """
        Multivariate Horner's scheme

        Parameters
        ----------
        x, y : array
        coeff : array of coefficients in inverse lexical order
        """

        alpha = self._invlex()
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

    def _invlex(self):
        c = []
        lencoeff = self.degree + 1
        for i in range(lencoeff):
            for j in range(lencoeff):
                if i + j <= self.degree:
                    c.append((j, i))
        return c[::-1]


class OrthoPolynomialBase(Polynomial2D):
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

    def __init__(self, x_degree, y_degree, x_domain=None, x_window=None,
                 y_domain=None, y_window=None, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        super(OrthoPolynomialBase, self).__init__(x_degree + y_degree,
                x_domain=x_domain, x_window=x_window, y_domain=y_domain,
                y_window=y_window, n_models=n_models,
                model_set_axis=model_set_axis, name=name, meta=meta, **params)

    def __repr__(self):
        return self._format_repr([self.x_degree, self.y_degree])

    def __str__(self):
        return self._format_str(*(('{0}-degree'.format(var.upper()),
                                  getattr(self, '{0}_degree'.format(var)))
                                  for var in self.inputs))

    @classmethod
    def _get_degree_attrs(cls):
        degree_attrs = tuple('{0}_degree'.format(input_)
                             for input_ in cls.inputs)
        doc_templ = "Maximum degree in the '{0}' variable."
        degree_docs = tuple(doc_templ.format(input_)
                            for input_ in cls.inputs)

        return zip(degree_attrs, degree_docs)

    @classmethod
    def _term_powers(cls, *degrees):
        # The original implementation increments the x degree fastest,
        # while itertools.product increments the y degree fastest,
        # hence the reverses
        ranges = (range(d + 1) for d in degrees[::-1])
        for ind in itertools.product(*ranges):
            yield ind[::-1]

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

    def _fcache(self, x, y):
        # TODO: Write a docstring explaining the actual purpose of this method
        """To be implemented by subclasses"""

        raise NotImplementedError("Subclasses should implement this")

    def evaluate(self, x, y, *coeffs):
        if self.x_domain is not None:
            x = poly_map_domain(x, self.x_domain, self.x_window)
        if self.y_domain is not None:
            y = poly_map_domain(y, self.y_domain, self.y_window)
        invcoeff = self.invlex_coeff()
        return self.imhorner(x, y, invcoeff)

    def prepare_inputs(self, x, y, **kwargs):
        inputs, format_info = \
                super(OrthoPolynomialBase, self).prepare_inputs(x, y, **kwargs)

        x, y = inputs

        if x.shape != y.shape:
            raise ValueError("Expected input arrays to have the same shape")

        return (x, y), format_info


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
    coefficient_template = '{0}_{1}_{2}'
    # The first item in the above template is the coeff_prefix

    def __init__(self, order, coeff_prefix, n_models=None,
                 model_set_axis=None, name=None, meta=None, **params):
        if order < 2 or order > 9:
            raise ValueError("Degree of polynomial must be 2 < deg < 9")

        super(_SIP1D, self).__init__(n_models=n_models,
                                     model_set_axis=model_set_axis,
                                     name=name, meta=meta, **params)

    def __repr__(self):
        return self._format_repr(args=[self.order, self.coeff_prefix])

    def __str__(self):
        return self._format_str(
            [('Order', self.order),
             ('Coeff. Prefix', self.coeff_prefix)])

    @classmethod
    def _get_degree_attrs(cls):
        # Note: coeff_prefix isn't really a 'degree', but we still need it to
        # generate the coefficient names for this model so we include it here
        return [('order', 'Order/degree of the polynomial'),
                ('coeff_prefix', 'Alphabetic prefix to the coefficient names')]

    @classmethod
    def _term_powers(cls, order, coeff_prefix):
        for exponents in polynomial_exponents(2, order):
            # Note: SIP distortion polynomials exclude monomial terms
            if sum(exponents) >= 2:
                yield (coeff_prefix,) + exponents

    def evaluate(self, x, y, *coeffs):
        # TODO: Rewrite this so that it uses a simpler method of determining
        # the matrix based on the number of given coefficients.
        mcoef = self._coeff_matrix(self.coeff_prefix, coeffs)
        return self._eval_sip(x, y, mcoef)

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


copyreg.pickle(_PolynomialModelMeta, _PolynomialModelMeta.__reduce__)
