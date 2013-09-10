# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module defines two classes that deal with parameters.

It is unlikely users will need to work with these classes directly, unless they
define their own models.
"""

from __future__ import division

import numbers

import numpy as np

from ..utils import misc
from .utils import InputParameterError


__all__ = ['Parameters', 'Parameter']


def _tofloat(value):
    """Convert a parameter to float or float array"""

    if misc.isiterable(value):
        try:
            _value = np.array(value, dtype=np.float)
            shape = _value.shape
        except (TypeError, ValueError):
            # catch arrays with strings or user errors like different
            # types of parameters in a parameter set
            raise InputParameterError(
                "Parameter of {0} could not be converted to "
                "float".format(type(value)))
    elif isinstance(value, bool):
        raise InputParameterError(
            "Expected parameter to be of numerical type, not boolean")
    elif isinstance(value, numbers.Number):
        _value = float(value)
        shape = ()
    else:
        raise InputParameterError(
            "Don't know how to convert parameter of {0} to "
            "float".format(type(value)))
    return _value, shape


def getval(self, name):
    return getattr(self, '_' + name).value


class Parameter(object):
    """
    Wraps individual parameters.

    This class represents a model's parameter (in a somewhat broad
    sense). To support multiple parameter sets, a parameter has a dimension
    (dim) and is a list-like object. It supports indexing so that individual
    parameters can be updated.
    To support some level of validation a parameter has a shape attribute.
    Parameter objects behave like numbers.

    Parameters
    ----------
    name : string
        parameter name
    val :  number or an iterable of numbers
    model : object
        an instance of a Model class
    dim : int
        parameter dimension
    fixed: boolean
        if True the parameter is not varied during fitting
    tied: callable or False
        if callable is suplplied it provides a way to link to another parameter
    min: float
        the lower bound of a parameter
    max: float
        the upper bound of a parameter
    """

    def __init__(self, name, val, model, dim, fixed=False, tied=False,
                 min=None, max=None):
        self._dim = dim
        self._name = name

        if self._dim == 1:
            val, parshape = _tofloat(val)
            self._value = val
        else:
            try:
                val0, parshape = _tofloat(val[0])
                val = [_tofloat(v)[0] for v in val]
                if len(val) != self._dim:
                    raise InputParameterError(
                        "Expected parameter {0} to be of dim {1}".format(
                            self._name, self._dim))
            except TypeError:
                raise InputParameterError("Expected a multivalued"
                                          " parameter {0}".format(self._name))
            for shape in [_tofloat(v)[1] for v in val]:
                assert shape == parshape, "Multiple values for the same" \
                    " parameters should have the same shape"

            self._value = val[:]

        self._shape = parshape
        self._model = model
        self._fixed = fixed
        self._tied = tied
        self._min = min
        self._max = max

    def __repr__(self):
        return repr(self.value)

    @property
    def dim(self):
        """Number of parameter sets"""

        return self._dim

    @dim.setter
    def dim(self, val):
        self._dim = val

    @property
    def value(self):
        """Parameter value"""

        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    @property
    def shape(self):
        """Parameter shape"""

        return self._shape

    @property
    def model(self):
        """An instance of `~astropy.modeling.core.ParametricModel`"""

        return self._model

    @model.setter
    def model(self, val):
        self._model = val

    @property
    def name(self):
        """Parameter name"""

        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def fixed(self):
        """
        Boolean indicating if the parameter is kept fixed during fitting.
        """

        return self._fixed

    @fixed.setter
    def fixed(self, val):
        """Fix a parameter."""

        assert isinstance(val, bool), "Fixed can be True or False"
        self._fixed = val

    @property
    def tied(self):
        """
        Indicates that this parameter is linked to another one.

        A callable which provides the relationship of the two parameters.
        """

        return self._tied

    @tied.setter
    def tied(self, val):
        """Tie a parameter."""

        assert callable(val) or val is False, "Tied must be a callable"
        self._tied = val

    @property
    def min(self):
        """A value used as a lower bound when fitting a parameter"""

        return self._min

    @min.setter
    def min(self, val):
        """Set a minimum value of a parameter."""

        if val is not None:
            assert isinstance(val, numbers.Number), "Min value must be a number"
            self._min = float(val)
        else:
            self._min = None

    @property
    def max(self):
        """A value used as an upper bound when fitting a parameter"""

        return self._max

    @max.setter
    def max(self, val):
        """Set a maximum value of a parameter."""

        if val is not None:
            assert isinstance(val, numbers.Number), "Max value must be a number"
            self._max = float(val)
        else:
            self._max = None

    def __getitem__(self, i):
        return self._value[i]

    def __setslice__(self, i, j, val):
        self._value[i:j] = _tofloat(val)[0]
        setattr(self.model, self.name, self._value)

    def __setitem__(self, i, val):
        val = _tofloat(val)[0]
        self._value[i] = val

        setattr(self.model, self.name, self._value)

    def __add__(self, val):
        return np.asarray(self._value) + val

    def __radd__(self, val):
        return np.asarray(self._value) + val

    def __sub__(self, val):
        return np.asarray(self._value) - val

    def __rsub__(self, val):
        return val - np.asarray(self._value)

    def __mul__(self, val):
        return np.asarray(self._value) * val

    def __rmul__(self, val):
        return np.asarray(self._value) * val

    def __pow__(self, val):
        return np.asarray(self._value) ** val

    def __div__(self, val):
        return np.asarray(self._value) / val

    def __rdiv__(self, val):
        return val / np.asarray(self._value)

    def __truediv__(self, val):
        return np.asarray(self._value) / val

    def __rtruediv__(self, val):
        return val / np.asarray(self._value)

    def __eq__(self, val):
        return (np.asarray(self._value) == np.asarray(val)).all()

    def __ne__(self, val):
        return not (np.asarray(self._value) == np.asarray(val)).all()

    def __lt__(self, val):
        return (np.asarray(self._value) < np.asarray(val)).all()

    def __gt__(self, val):
        return (np.asarray(self._value) > np.asarray(val)).all()

    def __le__(self, val):
        return (np.asarray(self._value) <= np.asarray(val)).all()

    def __ge__(self, val):
        return (np.asarray(self._value) >= np.asarray(val)).all()

    def __neg__(self):
        return np.asarray(self._value) * (-1)

    def __abs__(self):
        return np.abs(np.asarray(self._value))


class Parameters(list):
    """
    Store model parameters as a flat list of floats.

    This is a list-like object which stores model parameters. Only  instances
    of `~astropy.modeling.core.ParametricModel` keep an instance of this class
    as an attribute. The list of parameters can be modified by the user or by
    an instance of `~astropy.modeling.fitting.Fitter`.  This list of parameters
    is kept in sync with single model parameter attributes.  When more than one
    dimensional, a `~astropy.modeling.fitting.Fitter` treats each set of
    parameters as belonging to the same model but different set of data.

    Parameters
    ----------
    mobj : object
        an instance of a subclass of `~astropy.modeling.core.ParametricModel`
    param_names : list of strings
        parameter names
    dim : int
        Number of parameter sets
    """

    def __init__(self, mobj, param_names, dim=1):
        self.mobj = mobj
        self.dim = dim
        # A flag set to True by a fitter to indicate that the flat
        # list of parameters has been changed.
        self._changed = False
        self.parinfo = {}
        parlist = [getval(mobj, attr) for attr in param_names]
        flat = self._flatten(param_names, parlist)
        super(Parameters, self).__init__(flat)

    def __setitem__(self, ind, value):
        _val = _tofloat(value)[0]
        super(Parameters, self).__setitem__(ind, _val)
        self._changed = True
        self._update_model_params()

    def __setslice__(self, istart, istop, vlist):
        super(Parameters, self).__setslice__(istart, istop, vlist)
        self._changed = True
        self._update_model_params()

    def _update_model_params(self):
        """Update individual parameters"""

        for key in self.parinfo.keys():
            sl = self.parinfo[key][0]
            val = self[sl]
            if len(val) == 1:
                val = val[0]
            else:
                val = val
            par = getattr(self.mobj, key)
            setattr(par, 'value', val)
        self._changed = False

    def _is_same_length(self, newparams):
        """
        Checks if the user supplied value of
        `~astropy.modeling.core.ParametricModel.parameters` has the same length
        as the original parameters list.
        """

        parsize = _tofloat(newparams)[0].size
        return parsize == self.__len__()

    def _flatten(self, param_names, parlist):
        """Create a flat list of model parameters"""

        flatpars = []
        start = 0
        for (name, par) in zip(param_names, parlist):
            pararr = np.array(par)
            fpararr = pararr.flatten()

            stop = start + len(fpararr)
            self.parinfo[name] = slice(start, stop, 1), pararr.shape
            start = stop
            flatpars.extend(list(fpararr))

        return flatpars
