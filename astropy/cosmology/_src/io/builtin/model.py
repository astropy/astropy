# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""|Cosmology| <-> Model I/O, using |Cosmology.to_format| and |Cosmology.from_format|.

This module provides functions to transform a |Cosmology| object to and from a
:class:`~astropy.modeling.Model`. The functions are registered with ``convert_registry``
under the format name "astropy.model".

Using ``format="astropy.model"`` any redshift(s) method of a cosmology may be turned
into a :class:`astropy.modeling.Model`. Each |Cosmology|
:class:`~astropy.cosmology.Parameter` is converted to a :class:`astropy.modeling.Model`
:class:`~astropy.modeling.Parameter` and the redshift-method to the model's ``__call__ /
evaluate``. Now you can fit cosmologies with data!

.. code-block::

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> model = Planck18.to_format("astropy.model", method="lookback_time")
    >>> model
    <FlatLambdaCDMCosmologyLookbackTimeModel(H0=67.66 km / (Mpc s), Om0=0.30966,
        Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.  , 0.  , 0.06] eV, Ob0=0.04897,
        name='Planck18')>

The |Planck18| cosmology can be recovered with |Cosmology.from_format|.

    >>> print(Cosmology.from_format(model))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
"""

from __future__ import annotations

import abc
import copy
import inspect
from dataclasses import replace
from typing import Generic

import numpy as np

from astropy.modeling import FittableModel, Model
from astropy.utils.decorators import classproperty

# isort: split
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import convert_registry
from astropy.cosmology._src.typing import _CosmoT

from .utils import convert_parameter_to_model_parameter

__all__: list[str] = []  # nothing is publicly scoped


class _CosmologyModel(FittableModel, Generic[_CosmoT]):
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
    def _cosmology_class(self) -> type[_CosmoT]:
        """Cosmology class as a private attribute.

        Set in subclasses.
        """

    @abc.abstractmethod
    def _method_name(self) -> str:
        """Cosmology method name as a private attribute.

        Set in subclasses.
        """

    @classproperty
    def cosmology_class(cls) -> type[_CosmoT]:
        """|Cosmology| class."""
        return cls._cosmology_class

    @classproperty(lazy=True)
    def _cosmology_class_sig(cls):
        """Signature of |Cosmology| class."""
        return inspect.signature(cls._cosmology_class)

    @property
    def cosmology(self) -> _CosmoT:
        """Return |Cosmology| using `~astropy.modeling.Parameter` values."""
        cosmo = self._cosmology_class(
            name=self.name,
            **{
                k: (v.value if not (v := getattr(self, k)).unit else v.quantity)
                for k in self.param_names
            },
        )
        return cosmo

    @classproperty
    def method_name(self) -> str:
        """Redshift-method name on |Cosmology| instance."""
        return self._method_name

    # ---------------------------------------------------------------

    # NOTE: cannot add type annotations b/c of how Model introspects
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
        # TODO: speed up using ``replace``

        # create BoundArgument with all available inputs beyond the Parameters,
        # which will be filled in next
        ba = self._cosmology_class_sig.bind_partial(*args[self.n_inputs :], **kwargs)

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
        result = getattr(cosmo, self._method_name)(*args[: self.n_inputs])

        return result


##############################################################################


def from_model(model: _CosmologyModel[_CosmoT]) -> _CosmoT:
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
    >>> print(Cosmology.from_format(model))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
    """
    cosmo = model.cosmology

    # assemble the metadata
    meta = copy.deepcopy(model.meta)
    for n in model.param_names:
        p = getattr(model, n)
        meta[p.name] = {
            n: getattr(p, n)
            for n in dir(p)
            if not (n.startswith("_") or callable(getattr(p, n)))
        }
    return replace(cosmo, meta=meta)


def to_model(cosmology: _CosmoT, *_: object, method: str) -> _CosmologyModel[_CosmoT]:
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

    params = {
        k: convert_parameter_to_model_parameter(
            cosmo_cls.parameters[k], v, meta=cosmology.meta.get(k)
        )
        for k, v in cosmology.parameters.items()
        if v is not None
    }

    # class name is cosmology name + Cosmology + method name + Model
    clsname = (
        cosmo_cls.__qualname__.replace(".", "_")
        + "Cosmology"
        + method.replace("_", " ").title().replace(" ", "")
        + "Model"
    )

    # make Model class
    CosmoModel = type(clsname, (_CosmologyModel,), {**attrs, **params})
    # override __signature__ and format the doc.
    CosmoModel.evaluate.__signature__ = msig
    if CosmoModel.evaluate.__doc__ is not None:
        # guard against PYTHONOPTIMIZE mode
        CosmoModel.evaluate.__doc__ = CosmoModel.evaluate.__doc__.format(
            cosmo_cls=cosmo_cls.__qualname__, method=method
        )

    # instantiate class using default values
    model = CosmoModel(
        **cosmology.parameters, name=cosmology.name, meta=copy.deepcopy(cosmology.meta)
    )

    return model


def model_identify(
    origin: str, format: str | None, *args: object, **kwargs: object
) -> bool:
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
