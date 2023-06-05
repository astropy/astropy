# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Spline models and fitters."""
# pylint: disable=line-too-long, too-many-lines, too-many-arguments, invalid-name

import abc
import functools
import warnings

import numpy as np

from astropy.utils import isiterable
from astropy.utils.exceptions import AstropyUserWarning

from .core import FittableModel, ModelDefinitionError
from .parameters import Parameter

__all__ = [
    "Spline1D",
    "SplineInterpolateFitter",
    "SplineSmoothingFitter",
    "SplineExactKnotsFitter",
    "SplineSplrepFitter",
]
__doctest_requires__ = {"Spline1D": ["scipy"]}


class _Spline(FittableModel):
    """Base class for spline models."""

    _knot_names = ()
    _coeff_names = ()

    optional_inputs = {}

    def __init__(
        self,
        knots=None,
        coeffs=None,
        degree=None,
        bounds=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
    ):
        super().__init__(
            n_models=n_models, model_set_axis=model_set_axis, name=name, meta=meta
        )

        self._user_knots = False
        self._init_tck(degree)

        # Hack to allow an optional model argument
        self._create_optional_inputs()

        if knots is not None:
            self._init_spline(knots, coeffs, bounds)
        elif coeffs is not None:
            raise ValueError(
                "If one passes a coeffs vector one needs to also pass knots!"
            )

    @property
    def param_names(self):
        """
        Coefficient names generated based on the spline's degree and
        number of knots.
        """
        return tuple(list(self._knot_names) + list(self._coeff_names))

    @staticmethod
    def _optional_arg(arg):
        return f"_{arg}"

    def _create_optional_inputs(self):
        for arg in self.optional_inputs:
            attribute = self._optional_arg(arg)
            if hasattr(self, attribute):
                raise ValueError(
                    f"Optional argument {arg} already exists in this class!"
                )
            else:
                setattr(self, attribute, None)

    def _intercept_optional_inputs(self, **kwargs):
        new_kwargs = kwargs
        for arg in self.optional_inputs:
            if arg in kwargs:
                attribute = self._optional_arg(arg)
                if getattr(self, attribute) is None:
                    setattr(self, attribute, kwargs[arg])
                    del new_kwargs[arg]
                else:
                    raise RuntimeError(
                        f"{arg} has already been set, something has gone wrong!"
                    )

        return new_kwargs

    def evaluate(self, *args, **kwargs):
        """Extract the optional kwargs passed to call."""
        optional_inputs = kwargs
        for arg in self.optional_inputs:
            attribute = self._optional_arg(arg)

            if arg in kwargs:
                # Options passed in
                optional_inputs[arg] = kwargs[arg]
            elif getattr(self, attribute) is not None:
                # No options passed in and Options set
                optional_inputs[arg] = getattr(self, attribute)
                setattr(self, attribute, None)
            else:
                # No options passed in and No options set
                optional_inputs[arg] = self.optional_inputs[arg]

        return optional_inputs

    def __call__(self, *args, **kwargs):
        """
        Make model callable to model evaluation.
        """
        # Hack to allow an optional model argument
        kwargs = self._intercept_optional_inputs(**kwargs)

        return super().__call__(*args, **kwargs)

    def _create_parameter(self, name: str, index: int, attr: str, fixed=False):
        """
        Create a spline parameter linked to an attribute array.

        Parameters
        ----------
        name : str
            Name for the parameter
        index : int
            The index of the parameter in the array
        attr : str
            The name for the attribute array
        fixed : optional, bool
            If the parameter should be fixed or not
        """
        # Hack to allow parameters and attribute array to freely exchange values
        #   _getter forces reading value from attribute array
        #   _setter forces setting value to attribute array

        def _getter(value, model: "_Spline", index: int, attr: str):
            return getattr(model, attr)[index]

        def _setter(value, model: "_Spline", index: int, attr: str):
            getattr(model, attr)[index] = value
            return value

        getter = functools.partial(_getter, index=index, attr=attr)
        setter = functools.partial(_setter, index=index, attr=attr)

        default = getattr(self, attr)
        param = Parameter(
            name=name, default=default[index], fixed=fixed, getter=getter, setter=setter
        )
        # setter/getter wrapper for parameters in this case require the
        # parameter to have a reference back to its parent model
        param.model = self
        param.value = default[index]

        # Add parameter to model
        self.__dict__[name] = param

    def _create_parameters(self, base_name: str, attr: str, fixed=False):
        """
        Create a spline parameters linked to an attribute array for all
        elements in that array.

        Parameters
        ----------
        base_name : str
            Base name for the parameters
        attr : str
            The name for the attribute array
        fixed : optional, bool
            If the parameters should be fixed or not
        """
        names = []
        for index in range(len(getattr(self, attr))):
            name = f"{base_name}{index}"
            names.append(name)

            self._create_parameter(name, index, attr, fixed)

        return tuple(names)

    @abc.abstractmethod
    def _init_parameters(self):
        raise NotImplementedError("This needs to be implemented")

    @abc.abstractmethod
    def _init_data(self, knots, coeffs, bounds=None):
        raise NotImplementedError("This needs to be implemented")

    def _init_spline(self, knots, coeffs, bounds=None):
        self._init_data(knots, coeffs, bounds)
        self._init_parameters()

        # fill _parameters and related attributes
        self._initialize_parameters((), {})
        self._initialize_slices()

        # Calling this will properly fill the _parameter vector, which is
        #   used directly sometimes without being properly filled.
        _ = self.parameters

    def _init_tck(self, degree):
        self._c = None
        self._t = None
        self._degree = degree

    def __getstate__(self):
        return {
            "t": self._t,
            "c": self._c,
            "k": self._degree,
        }

    def __setstate__(self, state):
        return self.__init__(knots=state["t"], coeffs=state["c"], degree=state["k"])


class Spline1D(_Spline):
    """
    One dimensional Spline Model.

    Parameters
    ----------
    knots :  optional
        Define the knots for the spline. Can be 1) the number of interior
        knots for the spline, 2) the array of all knots for the spline, or
        3) If both bounds are defined, the interior knots for the spline
    coeffs : optional
        The array of knot coefficients for the spline
    degree : optional
        The degree of the spline. It must be 1 <= degree <= 5, default is 3.
    bounds : optional
        The upper and lower bounds of the spline.

    Notes
    -----
    Much of the functionality of this model is provided by
    `scipy.interpolate.BSpline` which can be directly accessed via the
    bspline property.

    Fitting for this model is provided by wrappers for:
    `scipy.interpolate.UnivariateSpline`,
    `scipy.interpolate.InterpolatedUnivariateSpline`,
    and `scipy.interpolate.LSQUnivariateSpline`.

    If one fails to define any knots/coefficients, no parameters will
    be added to this model until a fitter is called. This is because
    some of the fitters for splines vary the number of parameters and so
    we cannot define the parameter set until after fitting in these cases.

    Since parameters are not necessarily known at model initialization,
    setting model parameters directly via the model interface has been
    disabled.

    Direct constructors are provided for this model which incorporate the
    fitting to data directly into model construction.

    Knot parameters are declared as "fixed" parameters by default to
    enable the use of other `astropy.modeling` fitters to be used to
    fit this model.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.modeling.models import Spline1D
    >>> from astropy.modeling import fitting
    >>> np.random.seed(42)
    >>> x = np.linspace(-3, 3, 50)
    >>> y = np.exp(-x**2) + 0.1 * np.random.randn(50)
    >>> xs = np.linspace(-3, 3, 1000)

    A 1D interpolating spline can be fit to data:

    >>> fitter = fitting.SplineInterpolateFitter()
    >>> spl = fitter(Spline1D(), x, y)

    Similarly, a smoothing spline can be fit to data:

    >>> fitter = fitting.SplineSmoothingFitter()
    >>> spl = fitter(Spline1D(), x, y, s=0.5)

    Similarly, a spline can be fit to data using an exact set of interior knots:

    >>> t = [-1, 0, 1]
    >>> fitter = fitting.SplineExactKnotsFitter()
    >>> spl = fitter(Spline1D(), x, y, t=t)
    """

    n_inputs = 1
    n_outputs = 1
    _separable = True

    optional_inputs = {"nu": 0}

    def __init__(
        self,
        knots=None,
        coeffs=None,
        degree=3,
        bounds=None,
        n_models=None,
        model_set_axis=None,
        name=None,
        meta=None,
    ):
        super().__init__(
            knots=knots,
            coeffs=coeffs,
            degree=degree,
            bounds=bounds,
            n_models=n_models,
            model_set_axis=model_set_axis,
            name=name,
            meta=meta,
        )

    @property
    def t(self):
        """
        The knots vector.
        """
        if self._t is None:
            return np.concatenate(
                (np.zeros(self._degree + 1), np.ones(self._degree + 1))
            )
        else:
            return self._t

    @t.setter
    def t(self, value):
        if self._t is None:
            raise ValueError(
                "The model parameters must be initialized before setting knots."
            )
        elif len(value) == len(self._t):
            self._t = value
        else:
            raise ValueError(
                "There must be exactly as many knots as previously defined."
            )

    @property
    def t_interior(self):
        """
        The interior knots.
        """
        return self.t[self.degree + 1 : -(self.degree + 1)]

    @property
    def c(self):
        """
        The coefficients vector.
        """
        if self._c is None:
            return np.zeros(len(self.t))
        else:
            return self._c

    @c.setter
    def c(self, value):
        if self._c is None:
            raise ValueError(
                "The model parameters must be initialized before setting coeffs."
            )
        elif len(value) == len(self._c):
            self._c = value
        else:
            raise ValueError(
                "There must be exactly as many coeffs as previously defined."
            )

    @property
    def degree(self):
        """
        The degree of the spline polynomials.
        """
        return self._degree

    @property
    def _initialized(self):
        return self._t is not None and self._c is not None

    @property
    def tck(self):
        """
        Scipy 'tck' tuple representation.
        """
        return (self.t, self.c, self.degree)

    @tck.setter
    def tck(self, value):
        if self._initialized:
            if value[2] != self.degree:
                raise ValueError("tck has incompatible degree!")

            self.t = value[0]
            self.c = value[1]
        else:
            self._init_spline(value[0], value[1])

        # Calling this will properly fill the _parameter vector, which is
        #   used directly sometimes without being properly filled.
        _ = self.parameters

    @property
    def bspline(self):
        """
        Scipy bspline object representation.
        """
        from scipy.interpolate import BSpline

        return BSpline(*self.tck)

    @bspline.setter
    def bspline(self, value):
        from scipy.interpolate import BSpline

        if isinstance(value, BSpline):
            self.tck = value.tck
        else:
            self.tck = value

    @property
    def knots(self):
        """
        Dictionary of knot parameters.
        """
        return [getattr(self, knot) for knot in self._knot_names]

    @property
    def user_knots(self):
        """If the knots have been supplied by the user."""
        return self._user_knots

    @user_knots.setter
    def user_knots(self, value):
        self._user_knots = value

    @property
    def coeffs(self):
        """
        Dictionary of coefficient parameters.
        """
        return [getattr(self, coeff) for coeff in self._coeff_names]

    def _init_parameters(self):
        self._knot_names = self._create_parameters("knot", "t", fixed=True)
        self._coeff_names = self._create_parameters("coeff", "c")

    def _init_bounds(self, bounds=None):
        if bounds is None:
            bounds = [None, None]

        if bounds[0] is None:
            lower = np.zeros(self._degree + 1)
        else:
            lower = np.array([bounds[0]] * (self._degree + 1))

        if bounds[1] is None:
            upper = np.ones(self._degree + 1)
        else:
            upper = np.array([bounds[1]] * (self._degree + 1))

        if bounds[0] is not None and bounds[1] is not None:
            self.bounding_box = bounds
            has_bounds = True
        else:
            has_bounds = False

        return has_bounds, lower, upper

    def _init_knots(self, knots, has_bounds, lower, upper):
        if np.issubdtype(type(knots), np.integer):
            self._t = np.concatenate((lower, np.zeros(knots), upper))
        elif isiterable(knots):
            self._user_knots = True
            if has_bounds:
                self._t = np.concatenate((lower, np.array(knots), upper))
            else:
                if len(knots) < 2 * (self._degree + 1):
                    raise ValueError(
                        f"Must have at least {2*(self._degree + 1)} knots."
                    )
                self._t = np.array(knots)
        else:
            raise ValueError(f"Knots: {knots} must be iterable or value")

        # check that knots form a viable spline
        self.bspline  # noqa: B018

    def _init_coeffs(self, coeffs=None):
        if coeffs is None:
            self._c = np.zeros(len(self._t))
        else:
            self._c = np.array(coeffs)

        # check that coeffs form a viable spline
        self.bspline  # noqa: B018

    def _init_data(self, knots, coeffs, bounds=None):
        self._init_knots(knots, *self._init_bounds(bounds))
        self._init_coeffs(coeffs)

    def evaluate(self, *args, **kwargs):
        """
        Evaluate the spline.

        Parameters
        ----------
        x :
            (positional) The points where the model is evaluating the spline at
        nu : optional
            (kwarg) The derivative of the spline for evaluation, 0 <= nu <= degree + 1.
            Default: 0.
        """
        kwargs = super().evaluate(*args, **kwargs)
        x = args[0]

        if "nu" in kwargs:
            if kwargs["nu"] > self.degree + 1:
                raise RuntimeError(
                    "Cannot evaluate a derivative of "
                    f"order higher than {self.degree + 1}"
                )

        return self.bspline(x, **kwargs)

    def derivative(self, nu=1):
        """
        Create a spline that is the derivative of this one.

        Parameters
        ----------
        nu : int, optional
            Derivative order, default is 1.
        """
        if nu <= self.degree:
            bspline = self.bspline.derivative(nu=nu)

            derivative = Spline1D(degree=bspline.k)
            derivative.bspline = bspline

            return derivative
        else:
            raise ValueError(f"Must have nu <= {self.degree}")

    def antiderivative(self, nu=1):
        """
        Create a spline that is an antiderivative of this one.

        Parameters
        ----------
        nu : int, optional
            Antiderivative order, default is 1.

        Notes
        -----
        Assumes constant of integration is 0
        """
        if (nu + self.degree) <= 5:
            bspline = self.bspline.antiderivative(nu=nu)

            antiderivative = Spline1D(degree=bspline.k)
            antiderivative.bspline = bspline

            return antiderivative
        else:
            raise ValueError(
                "Supported splines can have max degree 5, "
                f"antiderivative degree will be {nu + self.degree}"
            )


class _SplineFitter(abc.ABC):
    """
    Base Spline Fitter.
    """

    def __init__(self):
        self.fit_info = {"resid": None, "spline": None}

    def _set_fit_info(self, spline):
        self.fit_info["resid"] = spline.get_residual()
        self.fit_info["spline"] = spline

    @abc.abstractmethod
    def _fit_method(self, model, x, y, **kwargs):
        raise NotImplementedError("This has not been implemented for _SplineFitter.")

    def __call__(self, model, x, y, z=None, **kwargs):
        model_copy = model.copy()
        if isinstance(model_copy, Spline1D):
            if z is not None:
                raise ValueError("1D model can only have 2 data points.")

            spline = self._fit_method(model_copy, x, y, **kwargs)

        else:
            raise ModelDefinitionError(
                "Only spline models are compatible with this fitter."
            )

        self._set_fit_info(spline)

        return model_copy


class SplineInterpolateFitter(_SplineFitter):
    """
    Fit an interpolating spline.
    """

    def _fit_method(self, model, x, y, **kwargs):
        weights = kwargs.pop("weights", None)
        bbox = kwargs.pop("bbox", [None, None])

        if model.user_knots:
            warnings.warn(
                "The current user specified knots maybe ignored for interpolating data",
                AstropyUserWarning,
            )
            model.user_knots = False

        if bbox != [None, None]:
            model.bounding_box = bbox

        from scipy.interpolate import InterpolatedUnivariateSpline

        spline = InterpolatedUnivariateSpline(
            x, y, w=weights, bbox=bbox, k=model.degree
        )

        model.tck = spline._eval_args
        return spline


class SplineSmoothingFitter(_SplineFitter):
    """
    Fit a smoothing spline.
    """

    def _fit_method(self, model, x, y, **kwargs):
        s = kwargs.pop("s", None)
        weights = kwargs.pop("weights", None)
        bbox = kwargs.pop("bbox", [None, None])

        if model.user_knots:
            warnings.warn(
                "The current user specified knots maybe ignored for smoothing data",
                AstropyUserWarning,
            )
            model.user_knots = False

        if bbox != [None, None]:
            model.bounding_box = bbox

        from scipy.interpolate import UnivariateSpline

        spline = UnivariateSpline(x, y, w=weights, bbox=bbox, k=model.degree, s=s)

        model.tck = spline._eval_args
        return spline


class SplineExactKnotsFitter(_SplineFitter):
    """
    Fit a spline using least-squares regression.
    """

    def _fit_method(self, model, x, y, **kwargs):
        t = kwargs.pop("t", None)
        weights = kwargs.pop("weights", None)
        bbox = kwargs.pop("bbox", [None, None])

        if t is not None:
            if model.user_knots:
                warnings.warn(
                    "The current user specified knots will be "
                    "overwritten for by knots passed into this function",
                    AstropyUserWarning,
                )
        else:
            if model.user_knots:
                t = model.t_interior
            else:
                raise RuntimeError("No knots have been provided")

        if bbox != [None, None]:
            model.bounding_box = bbox

        from scipy.interpolate import LSQUnivariateSpline

        spline = LSQUnivariateSpline(x, y, t, w=weights, bbox=bbox, k=model.degree)

        model.tck = spline._eval_args
        return spline


class SplineSplrepFitter(_SplineFitter):
    """
    Fit a spline using the `scipy.interpolate.splrep` function interface.
    """

    def __init__(self):
        super().__init__()
        self.fit_info = {"fp": None, "ier": None, "msg": None}

    def _fit_method(self, model, x, y, **kwargs):
        t = kwargs.pop("t", None)
        s = kwargs.pop("s", None)
        task = kwargs.pop("task", 0)
        weights = kwargs.pop("weights", None)
        bbox = kwargs.pop("bbox", [None, None])

        if t is not None:
            if model.user_knots:
                warnings.warn(
                    "The current user specified knots will be "
                    "overwritten for by knots passed into this function",
                    AstropyUserWarning,
                )
        else:
            if model.user_knots:
                t = model.t_interior

        if bbox != [None, None]:
            model.bounding_box = bbox

        from scipy.interpolate import splrep

        tck, fp, ier, msg = splrep(
            x,
            y,
            w=weights,
            xb=bbox[0],
            xe=bbox[1],
            k=model.degree,
            s=s,
            t=t,
            task=task,
            full_output=1,
        )
        model.tck = tck
        return fp, ier, msg

    def _set_fit_info(self, spline):
        self.fit_info["fp"] = spline[0]
        self.fit_info["ier"] = spline[1]
        self.fit_info["msg"] = spline[2]
