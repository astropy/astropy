# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains dictionaries with sets of parameters for a
given cosmology.

Each cosmology has the following parameters defined:

    ==========  =====================================
    Oc0         Omega cold dark matter at z=0
    Ob0         Omega baryon at z=0
    Om0         Omega matter at z=0
    flat        Is this assumed flat?  If not, Ode0 must be specified
    Ode0        Omega dark energy at z=0 if flat is False
    H0          Hubble parameter at z=0 in km/s/Mpc
    n           Density perturbation spectral index
    Tcmb0       Current temperature of the CMB
    Neff        Effective number of neutrino species
    m_nu        Assumed mass of neutrino species, in eV.
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr
    ==========  =====================================

The list of cosmologies available are given by the tuple
``available``. Current cosmologies available:

Planck 2018 (Planck18_arXiv_v2) parameters from Planck Collaboration 2018,
 arXiv:1807.06209v2 (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)

Planck 2015 (Planck15) parameters from Planck Collaboration 2016, A&A, 594, A13
 (Paper XIII), Table 4 (TT, TE, EE + lowP + lensing + ext)

Planck 2013 (Planck13) parameters from Planck Collaboration 2014, A&A, 571, A16
 (Paper XVI), Table 5 (Planck + WP + highL + BAO)

WMAP 9 year (WMAP9) parameters from Hinshaw et al. 2013, ApJS, 208, 19,
doi: 10.1088/0067-0049/208/2/19. Table 4 (WMAP9 + eCMB + BAO + H0)

WMAP 7 year (WMAP7) parameters from Komatsu et al. 2011, ApJS, 192, 18,
doi: 10.1088/0067-0049/192/2/18. Table 1 (WMAP + BAO + H0 ML).

WMAP 5 year (WMAP5) parameters from Komatsu et al. 2009, ApJS, 180, 330,
doi: 10.1088/0067-0049/180/2/330. Table 1 (WMAP + BAO + SN ML).

"""
import copy
import inspect
import sys
import typing as T
import warnings
from types import MappingProxyType

from .core import Cosmology, LambdaCDM, FlatLambdaCDM
from .. import units as u
from ..utils.state import ScienceState, _StateProxy
from ..utils.decorators import classproperty

from ..utils.exceptions import AstropyUserWarning

# Check pkg_resources exists
try:
    from pkg_resources import iter_entry_points

    HAS_PKG = True
except ImportError:
    HAS_PKG = False

__all__ = [
    "default_cosmology",  # *available (added below)
]


# Note: if you add a new cosmology, please also update the table
# in the 'Built-in Cosmologies' section of astropy/docs/cosmology/index.rst
# in addition to the list above.

_parameter_registry = dict(
    # Planck 2018 paper VI v2
    # Unlike Planck 2015, the paper includes massive neutrinos in Om0,
    # which here are included in m_nu.  Hence, the Om0 value differs slightly
    # from the paper.
    Planck18_arXiv_v2=_StateProxy(
        dict(
            parameters=dict(
                name="Planck18_arXiv_v2",
                Oc0=0.2607,
                Ob0=0.04897,
                Om0=0.30966,
                H0=67.66,
                n=0.9665,
                sigma8=0.8102,
                tau=0.0561,
                z_reion=7.82,
                t0=13.787,
                Tcmb0=2.7255,
                Neff=3.046,
                flat=True,
                m_nu=[0.0, 0.0, 0.06],
            ),
            references=(
                "Planck Collaboration 2018, arXiv:1807.06209 v2 (Paper VI),"
                " Table 2 (TT, TE, EE + lowE + lensing + BAO)"
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
    # Planck 2015 paper XII Table 4 final column (best fit)
    Planck15=_StateProxy(
        dict(
            parameters=dict(
                name="Planck15",
                Oc0=0.2589,
                Ob0=0.04860,
                Om0=0.3075,
                H0=67.74,
                n=0.9667,
                sigma8=0.8159,
                tau=0.066,
                z_reion=8.8,
                t0=13.799,
                Tcmb0=2.7255,
                Neff=3.046,
                flat=True,
                m_nu=[0.0, 0.0, 0.06],
            ),
            references=(
                "Planck Collaboration 2016, A&A, 594, A13 (Paper XIII),"
                " Table 4 (TT, TE, EE + lowP + lensing + ext)"
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
    # Planck 2013 paper XVI Table 5 penultimate column (best fit)
    Planck13=_StateProxy(
        dict(
            parameters=dict(
                name="Planck13",
                Oc0=0.25886,
                Ob0=0.048252,
                Om0=0.30712,
                H0=67.77,
                n=0.9611,
                sigma8=0.8288,
                tau=0.0952,
                z_reion=11.52,
                t0=13.7965,
                Tcmb0=2.7255,
                Neff=3.046,
                flat=True,
                m_nu=[0.0, 0.0, 0.06],
            ),
            references=(
                "Planck Collaboration 2014, A&A, 571, A16 (Paper XVI),"
                " Table 5 (Planck + WP + highL + BAO)"
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
    WMAP9=_StateProxy(
        dict(
            parameters=dict(
                name="WMAP9",
                Oc0=0.2402,
                Ob0=0.04628,
                Om0=0.2865,
                H0=69.32,
                n=0.9608,
                sigma8=0.820,
                tau=0.081,
                z_reion=10.1,
                t0=13.772,
                Tcmb0=2.725,
                Neff=3.04,
                m_nu=0.0,
                flat=True,
            ),
            references=(
                "Hinshaw et al. 2013, ApJS, 208, 19, "
                "doi: 10.1088/0067-0049/208/2/19. "
                "Table 4 (WMAP9 + eCMB + BAO + H0, last column)"
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
    WMAP7=_StateProxy(
        dict(
            parameters=dict(
                name="WMAP7",
                Oc0=0.226,
                Ob0=0.0455,
                Om0=0.272,
                H0=70.4,
                n=0.967,
                sigma8=0.810,
                tau=0.085,
                z_reion=10.3,
                t0=13.76,
                Tcmb0=2.725,
                Neff=3.04,
                m_nu=0.0,
                flat=True,
            ),
            references=(
                "Komatsu et al. 2011, ApJS, 192, 18, "
                "doi: 10.1088/0067-0049/192/2/18. "
                "Table 1 (WMAP + BAO + H0 ML)."
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
    WMAP5=_StateProxy(
        dict(
            parameters=dict(
                name="WMAP5",
                Oc0=0.231,
                Ob0=0.0459,
                Om0=0.277,
                H0=70.2,
                n=0.962,
                sigma8=0.817,
                tau=0.088,
                z_reion=11.3,
                t0=13.72,
                Tcmb0=2.725,
                Neff=3.04,
                m_nu=0.0,
                flat=True,
            ),
            references=(
                "Komatsu et al. 2009, ApJS, 180, 330, "
                "doi: 10.1088/0067-0049/180/2/330. "
                "Table 1 (WMAP + BAO + SN ML)."
            ),
            cosmo=FlatLambdaCDM,
        )
    ),
)


available = _parameter_registry.keys()  # auto-updates to reflect registry
__all__ += list(available)  # add all available at import


def _parse_cosmology_in_registry(name: str):
    """
    Parse cosmology from string name

    Parameters
    ----------
    name : str
        registry name

    Returns
    -------
    cosmo : :class:`~astropy.cosmology.core.Cosmology` instance
        With docstring set from name and references

    """
    if name not in _parameter_registry:
        s = "Unknown cosmology '{}'. Valid cosmologies:\n{}".format(
            name, available
        )
        raise ValueError(s)

    state: dict = copy.deepcopy(_parameter_registry[name])

    params: dict = state["parameters"]
    references = state.get("references", None)
    cls = state.get("cosmo", None)

    # Determine class, if None
    if cls is None:
        if params["flat"]:  # else try to grok  # TODO, better.
            cls = FlatLambdaCDM
        else:
            cls = LambdaCDM

    # parse parameters
    if "m_nu" in params:
        params["m_nu"] = u.Quantity(params["m_nu"], u.eV)

    # Create class instance
    sig = inspect.signature(cls.__init__)
    sig = sig.replace(parameters=list(sig.parameters.values())[1:])
    ba = sig.bind_partial(
        **{k: p for k, p in params.items() if k in sig.parameters.keys()},
    )  # apply any arguments that can be applied

    cosmo = cls(*ba.args, **ba.kwargs)
    docstr: str = "{} instance of {} cosmology\n\n(from {})"
    cosmo.__doc__ = docstr.format(name, cls.__name__, references)

    return cosmo


#########################################################################
# The science state below contains the current cosmology.
#########################################################################


class default_cosmology(ScienceState):
    """
    This class controls the default cosmology to use.

    This class controls the parameter settings by specifying a string name,
    or :class:`~astropy.cosmology.core.Cosmology` instance with the
    pre-specified options listed in ``available``

    Alternatively, user-defined cosmologies may be registered, with
    ``default_cosmology.register_parameters`` or
    ``default_cosmology.register_cosmology_instance``,
    and used identically as pre-specified options.

    Examples
    --------
    First import two Cosmology classes, an instance, and the ScienceState.

        >>> from astropy.cosmology import (Cosmology, LambdaCDM,
        ...                                default_cosmology, WMAP7)

    To change the default cosmology::

        >>> default_cosmology.set(WMAP7)   # doctest: +FLOAT_CMP
        <ScienceState default_cosmology:
            FlatLambdaCDM(name="WMAP7", H0=70.4 km / (Mpc s), Om0=0.272,
                Tcmb0=2.725 K, Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.0455)>

    Or, you may use a string::

        >>> _ = default_cosmology.set('Planck15')

    The default cosmology can be retrieved by the ``get`` method

        >>> cosmo = default_cosmology.get()
        >>> print(cosmo.name)
        Planck15

    Or just by calling :class:`~astropy.cosmology.core.Cosmology`

        >>> cosmo = Cosmology()
        >>> print(cosmo.name)
        Planck15

    The default cosmology can also be updated by using this class as a context
    manager::

        >>> with default_cosmology.set('WMAP9'):
        ...    cosmo = Cosmology()
        ...    print(cosmo.name)
        WMAP9

    Again, changing the default parameter values will not affect frame
    attributes that are explicitly specified (except the name)::

        >>> with default_cosmology.set('Planck15'):
        ...     cosmo = Cosmology(H0=68)
        ...     print(cosmo.name, cosmo.H0, sep=", ")
        None, 68.0 km / (Mpc s)


    Additional cosmologies may be registered from the parameters::

        >>> parameters = dict(
        ...     Oc0=0.231, Ob0=0.0459, Om0=0.277, Ode0=0.4461, z_reion=11.3,
        ...     H0=70.2, n=0.962, sigma8=0.817, tau=0.088, t0=13.72,
        ...     Tcmb0=2.725, Neff=3.04, m_nu=0.0, flat=False)
        >>> default_cosmology.register_parameters(
        ...     "example", parameters=parameters, cosmo=LambdaCDM)

    Or from a class instance::

        >>> cosmo = LambdaCDM(name="ex_instance", H0=68, Om0=0.3, Ode0=0.5)
        >>> default_cosmology.register_cosmology_instance(cosmo)

    """

    _default_value: str = "Planck15"
    _value = _default_value[:]  # set to default (copy)
    _registry = MappingProxyType(_parameter_registry)  # view only

    @classproperty
    def available(cls):
        """Registered cosmologies."""
        return available

    @classmethod
    def get_from_registry(cls, name: str):
        """
        Return cosmology parameters and metadata given string name.
        This method ensures the returned state is mutable,
        but the registry state remains unchanged.
        This is the preferred means of accessing the registry.

        Returns
        -------
        state : dict
            Copy of the registry for the string name.
            At minimum, includes the fields "parameters", "references", "cosmo"
            which hold the cosmology parameters, class, and references.

            - "parameters": dict
                cosmology parameters
            - "references" : dict or str or None
                References for "parameters".
                Fields are str or sequence of str.
            - "cosmo" : :class:`~astropy.cosmology.core.Cosmology`

        Raises
        ------
        KeyError
            If invalid string input to registry

        """
        # Get the state from the registry.
        # Copy to ensure registry is immutable to modifications of "_value".
        # Raises KeyError if ``name`` is invalid string input to registry
        state = copy.deepcopy(cls._registry[name])  # decouple & mutable

        # Get references form the state
        references = state.get("references", None) or {}  # get, if has
        references = references or None  # empty references -> None
        state["references"] = references

        cosmo = state.get("cosmo", None)
        state["cosmo"] = cosmo

        return state

    @staticmethod
    def get_cosmology_from_string(name: str) -> T.Optional[Cosmology]:
        """ Return a cosmology instance from a string.

        Parameters
        ----------
        name : str
            Name of the cosmology in registry
            or importable from :mod:`~astropy.cosmology.builtin`.

        Returns
        -------
        cosmo : :class:`~astropy.cosmology.core.Cosmology` instance

        """
        if name == "no_default":
            cosmo = None
        else:
            try:  # get from "cache"
                cosmo = getattr(sys.modules[__name__], name)
            except AttributeError:  # not already defined
                # make from registry
                cosmo = _parse_cosmology_in_registry(name)
                # cache so importable
                setattr(sys.modules[__name__], name, cosmo)

        return cosmo

    @classmethod
    def validate(cls, value: T.Union[Cosmology, str, None]) -> Cosmology:
        """Validate a Cosmology.

        Parameters
        ----------
        value: :class:`~astropy.cosmology.core.Cosmology` or str or None
            None becomes default value
            str calls ``get_cosmology_from_string``
            Cosmology instance is passed through

        Returns
        -------
        cosmo : :class:`~astropy.cosmology.core.Cosmology` or None
            None if value was "no_default"
            returned cosmology is set to ScienceState value.

        Raises
        ------
        TypeError
            If value is not None, str, or Cosmology instance.

        """
        if value is None:
            value = cls._default_value

        if isinstance(value, str):
            cosmo = cls.get_cosmology_from_string(value)
        elif isinstance(value, Cosmology):
            cosmo = value
        else:
            raise TypeError(
                "default_cosmology must be a string or Cosmology instance."
            )

        return cosmo

    @classmethod
    def register_parameters(
        cls,
        name: str,
        parameters: dict,
        cosmo: T.Optional[Cosmology] = None,
        references=None,
        *,
        viewonly=True,
    ):
        """
        Register a set of parameters.

        Parameters
        ----------
        name : str
            The registration name for the parameter and metadata set.
        parameters : dict
            The cosomological parameter.
        reference : dict or None, optional
            References for contents of parameters.
        viewonly : bool, optional, kwarg only
            whether to make registrant read-only.

        Notes
        -----
        Parameters must have:

            ==========  =====================================
            Oc0         Omega cold dark matter at z=0
            Ob0         Omega baryon at z=0
            Om0         Omega matter at z=0
            flat        Is this assumed flat?  If not, Ode0 must be specified
            Ode0        Omega dark energy at z=0 if flat is False
            H0          Hubble parameter at z=0 in km/s/Mpc
            n           Density perturbation spectral index
            Tcmb0       Current temperature of the CMB
            Neff        Effective number of neutrino species
            m_nu        Assumed mass of neutrino species, in eV.
            sigma8      Density perturbation amplitude
            tau         Ionisation optical depth
            z_reion     Redshift of hydrogen reionisation
            t0          Age of the universe in Gyr
            ==========  =====================================

        """
        parameters["name"] = name  # ensure name in parameters

        # check on contents of parameters
        must_have: tuple = (
            ("Oc0", "Ob0", "Om0")  # densities
            + ("H0", "n", "Tcmb0", "Neff", "m_nu", "sigma8", "tau")  # ps
            + ("z_reion", "t0")  # times
            + ("flat",)
        )
        if parameters["flat"] is False:
            must_have += ("Ode0",)

        missing: list = [k for k in must_have if k not in parameters]
        if missing:
            raise ValueError(f"Missing parameters {missing}")

        state = dict(
            parameters=copy.deepcopy(parameters),
            references=references,
            cosmo=cosmo,
        )
        if viewonly:
            state = _StateProxy(state)

        _parameter_registry[name] = state

    @classmethod
    def register_cosmology_instance(
        cls,
        cosmo: Cosmology,
        # references=None,
        viewonly=True,
    ):
        """
        Register parameters from a cosmology instance.

        Parameters
        ----------
        name : str
            The registration name for the parameter and metadata set.
        cosmo : :class:`~astropy.cosmology.core.Cosmology` instance

        """
        if cosmo.name is None:
            raise ValueError("Need to name cosmology.")
        elif cosmo.name in cls._registry:
            warnings.warn(
                AstropyUserWarning("Overwriting existing cosmology.")
            )

        cls = cosmo.__class__  # get class from cosmology
        if cls is Cosmology:
            raise TypeError("Cannot register the base cosmology class.")
        # TODO get references
        references = None

        # Get parameters
        # inspect the class, dropping "self"
        sig = inspect.signature(cosmo.__class__.__init__)
        sig = sig.replace(parameters=list(sig.parameters.values())[1:])
        # unpack the parameters used to initialize the class
        ba = sig.bind_partial(**cosmo._init_params)
        parameters = ba.arguments

        parameters["name"] = cosmo.name

        state = dict(parameters=parameters, references=references, cosmo=cls,)
        if viewonly:
            state = _StateProxy(state)

        _parameter_registry[cosmo.name] = state


#########################################################################
# Create Pre-Defined Cosmologies
#########################################################################

# TODO merge this method with the one in modeling.fitting as a class
def populate_entry_points(entry_points):
    """
    Inject entry points into the :mod:`~astropy.cosmology.builtin` namespace.
    This provides a means of inserting a cosmology instance without requiring
    it being merged into astropy's core.

    Parameters
    ----------

    entry_points : a list of `~pkg_resources.EntryPoint`
                  entry_points are objects which encapsulate
                  importable objects and are defined on the
                  installation of a package.

    Notes
    -----
    An explanation of entry points can be found `here <http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`

    """
    for entry_point in entry_points:
        name = entry_point.name
        try:
            entry_point = entry_point.load()
        except Exception as e:
            # Stops choking if an entry_point produces an error.
            warnings.warn(
                AstropyUserWarning(
                    f"{type(e).__name__} error occurred in entry point {name}."
                )
            )
        else:
            if isinstance(entry_point, Cosmology):
                default_cosmology.register_cosmology_instance(entry_point)
            else:
                warnings.warn(
                    AstropyUserWarning(
                        f"Cosmology entry point {name} expected to extend "
                        "astropy.cosmology.builtin._parameter_registry"
                    )
                )


# this is so fitting doesn't choke if pkg_resources doesn't exist
if HAS_PKG:
    populate_entry_points(
        iter_entry_points(group="astropy.cosmology", name=None)
    )


key: str
for key in available:
    cosmo = _parse_cosmology_in_registry(key)
    setattr(sys.modules[__name__], key, cosmo)

del key  # don't leave variable floating around in the namespace
