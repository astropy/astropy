# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The following are private functions, included here **FOR REFERENCE ONLY** since
the io registry cannot be displayed. These functions are registered into
:meth:`~astropy.cosmology.Cosmology.to_format` and
:meth:`~astropy.cosmology.Cosmology.from_format` and should only be accessed
via these methods.
"""  # this is shown in the docs.

import abc
import copy
import inspect

import numpy as np

from astropy.cosmology.connect import convert_registry
from astropy.cosmology.core import Cosmology
from astropy.modeling import FittableModel, Model
from astropy.utils.decorators import classproperty

from .utils import convert_parameter_to_model_parameter

__all__ = []  # nothing is publicly scoped


class _CosmologyModel(FittableModel):
    """Base class for Cosmology redshift-method Models.

    .. note::

        This class is not publicly scoped so should not be used directly.
        Instead, from a Cosmology instance use ``.to_format("astropy.model")``
        to create an instance of a subclass of this class.

    `_CosmologyModel` (subclasses) wrap a redshift-method of a
    :class:`~astropy.cosmology.Cosmology` class, converting each non-`None`
    |Cosmology| :class:`~astropy.cosmology.Parameter` to a
    :class:`astropy.modeling.Model` :class:`~astropy.modeling.Parameter`
    and the redshift-method to the model's ``__call__ / evaluate``.

    See Also
    --------
    astropy.cosmology.Cosmology.to_format
    """

    @abc.abstractmethod
    def _cosmology_class(self):
        """Cosmology class as a private attribute. Set in subclasses."""

    @abc.abstractmethod
    def _method_name(self):
        """Cosmology method name as a private attribute. Set in subclasses."""

    @classproperty
    def cosmology_class(cls):
        """|Cosmology| class."""
        return cls._cosmology_class

    @property
    def cosmology(self):
        """Return |Cosmology| using `~astropy.modeling.Parameter` values."""
        cosmo = self._cosmology_class(
            name=self.name,
            **{k: (v.value if not (v := getattr(self, k)).unit else v.quantity)
               for k in self.param_names})
        return cosmo

    @classproperty
    def method_name(self):
        """Redshift-method name on |Cosmology| instance."""
        return self._method_name

    # ---------------------------------------------------------------

    def evaluate(self, *args, **kwargs):
        """Evaluate method {method!r} of {cosmo_cls!r} Cosmology.

        The Model wraps the :class:`~astropy.cosmology.Cosmology` method,
        converting each |Cosmology| :class:`~astropy.cosmology.Parameter` to a
        :class:`astropy.modeling.Model` :class:`~astropy.modeling.Parameter`
        (unless the Parameter is None, in which case it is skipped).
        Here an instance of the cosmology is created using the current
        Parameter values and the method is evaluated given the input.

        Parameters
        ----------
        *args, **kwargs
            The first ``n_inputs`` of ``*args`` are for evaluating the method
            of the cosmology. The remaining args and kwargs are passed to the
            cosmology class constructor.
            Any unspecified Cosmology Parameter use the current value of the
            corresponding Model Parameter.

        Returns
        -------
        Any
            Results of evaluating the Cosmology method.
        """
        # create BoundArgument with all available inputs beyond the Parameters,
        # which will be filled in next
        ba = self.cosmology_class._init_signature.bind_partial(*args[self.n_inputs:], **kwargs)

        # fill in missing Parameters
        for k in self.param_names:
            if k not in ba.arguments:
                v = getattr(self, k)
                ba.arguments[k] = v.value if not v.unit else v.quantity

            # unvectorize, since Cosmology is not vectorized
            # TODO! remove when vectorized
            if np.shape(ba.arguments[k]):  # only in __call__
                # m_nu is a special case  # TODO! fix by making it 'structured'
                if k == "m_nu" and len(ba.arguments[k].shape) == 1:
                    continue
                ba.arguments[k] = ba.arguments[k][0]

        # make instance of cosmology
        cosmo = self._cosmology_class(**ba.arguments)
        # evaluate method
        result = getattr(cosmo, self._method_name)(*args[:self.n_inputs])

        return result


##############################################################################


def from_model(model):
    """Load |Cosmology| from `~astropy.modeling.Model` object.

    Parameters
    ----------
    model : `_CosmologyModel` subclass instance
        See ``Cosmology.to_format.help("astropy.model") for details.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    >>> from astropy.cosmology import Cosmology, Planck18
    >>> model = Planck18.to_format("astropy.model", method="lookback_time")
    >>> Cosmology.from_format(model)
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
    """
    cosmology = model.cosmology_class
    meta = copy.deepcopy(model.meta)

    # assemble the Parameters
    params = {}
    for n in model.param_names:
        p = getattr(model, n)
        params[p.name] = p.quantity if p.unit else p.value
        # put all attributes in a dict
        meta[p.name] = {n: getattr(p, n) for n in dir(p)
                        if not (n.startswith("_") or callable(getattr(p, n)))}

    ba = cosmology._init_signature.bind(name=model.name, **params, meta=meta)
    return cosmology(*ba.args, **ba.kwargs)


def to_model(cosmology, *_, method):
    """Convert a `~astropy.cosmology.Cosmology` to a `~astropy.modeling.Model`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    method : str, keyword-only
        The name of the method on the ``cosmology``.

    Returns
    -------
    `_CosmologyModel` subclass instance
        The Model wraps the |Cosmology| method, converting each non-`None`
        :class:`~astropy.cosmology.Parameter` to a
        :class:`astropy.modeling.Model` :class:`~astropy.modeling.Parameter`
        and the method to the model's ``__call__ / evaluate``.

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> model = Planck18.to_format("astropy.model", method="lookback_time")
    >>> model
    <FlatLambdaCDMCosmologyLookbackTimeModel(H0=67.66 km / (Mpc s), Om0=0.30966,
        Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.  , 0.  , 0.06] eV, Ob0=0.04897,
        name='Planck18')>
    """
    cosmo_cls = cosmology.__class__

    # get bound method & sig from cosmology (unbound if class).
    if not hasattr(cosmology, method):
        raise AttributeError(f"{method} is not a method on {cosmology.__class__}.")
    func = getattr(cosmology, method)
    if not callable(func):
        raise ValueError(f"{cosmology.__class__}.{method} is not callable.")
    msig = inspect.signature(func)

    # introspect for number of positional inputs, ignoring "self"
    n_inputs = len([p for p in tuple(msig.parameters.values()) if (p.kind in (0, 1))])

    attrs = {}  # class attributes
    attrs["_cosmology_class"] = cosmo_cls
    attrs["_method_name"] = method
    attrs["n_inputs"] = n_inputs
    attrs["n_outputs"] = 1

    params = {}  # Parameters (also class attributes)
    for n in cosmology.__parameters__:
        v = getattr(cosmology, n)  # parameter value

        if v is None:  # skip unspecified parameters
            continue

        # add as Model Parameter
        params[n] = convert_parameter_to_model_parameter(getattr(cosmo_cls, n), v,
                                                         cosmology.meta.get(n))

    # class name is cosmology name + Cosmology + method name + Model
    clsname = (cosmo_cls.__qualname__.replace(".", "_")
               + "Cosmology"
               + method.replace("_", " ").title().replace(" ", "")
               + "Model")

    # make Model class
    CosmoModel = type(clsname, (_CosmologyModel, ), {**attrs, **params})
    # override __signature__ and format the doc.
    setattr(CosmoModel.evaluate, "__signature__", msig)
    CosmoModel.evaluate.__doc__ = CosmoModel.evaluate.__doc__.format(
        cosmo_cls=cosmo_cls.__qualname__, method=method)

    # instantiate class using default values
    ps = {n: getattr(cosmology, n) for n in params.keys()}
    model = CosmoModel(**ps, name=cosmology.name, meta=copy.deepcopy(cosmology.meta))

    return model


def model_identify(origin, format, *args, **kwargs):
    """Identify if object uses the :class:`~astropy.modeling.Model` format.

    Returns
    -------
    bool
    """
    itis = False
    if origin == "read":
        itis = isinstance(args[1], Model) and (format in (None, "astropy.model"))

    return itis


# ===================================================================
# Register

convert_registry.register_reader("astropy.model", Cosmology, from_model)
convert_registry.register_writer("astropy.model", Cosmology, to_model)
convert_registry.register_identifier("astropy.model", Cosmology, model_identify)
