from abc import abstractmethod
from dataclasses import dataclass
from typing import TypeVar

from astropy.cosmology.core_wrapper import (
    AstropyCosmologyWrapper,
    CosmologyWrapperBase,
    _GetFromCosmo,
)

from .base import FLRW, FLRWAPI

Self = TypeVar("Self", bound="FLRWCosmologyWrapperBase")


@dataclass(frozen=True)
class FLRWCosmologyWrapperBase(CosmologyWrapperBase, FLRWAPI):
    """A base class for wrappers for astropy FLRW cosmology objects.

    This is a base class for wrappers for astropy cosmology objects that
    implement the FLRWAPI interface. It is used to allow the astropy
    cosmology objects to be used in the same way as the FLRWAPI objects.

    Parameters
    ----------
    cosmo : `astropy.cosmology.FLRW`
        The astropy cosmology object to wrap.

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> from lsst.afw.math import FLRWWrapperBase
    >>> cosmo = FLRWWrapperBase(Planck18)
    >>> cosmo.h
    .6766
    """

    cosmo: FLRWAPI

    # ==================================
    # Abstract methods from CosmologyAPI

    @property
    @abstractmethod
    def meta(self):
        ...

    @property
    @abstractmethod
    def name(self):
        ...

    @property
    def is_flat(self):
        return FLRWAPI.is_flat(self)

    @property
    @abstractmethod
    def clone(self, *, meta=None, **kwargs):
        ...

    # ==================================
    # Abstract methods from FLRWAPI

    @property
    @abstractmethod
    def H0(self):
        ...

    @property
    @abstractmethod
    def Om0(self):
        ...

    @property
    @abstractmethod
    def Ode0(self):
        ...

    @property
    @abstractmethod
    def Tcmb0(self):
        ...

    @property
    @abstractmethod
    def Neff(self):
        ...

    @property
    @abstractmethod
    def m_nu(self):
        ...

    @property
    @abstractmethod
    def Ob0(self):
        ...

    # --- Derived parameters ---

    @property
    @abstractmethod
    def Odm0(self):
        ...

    @property
    @abstractmethod
    def Ok0(self):
        ...

    @property
    @abstractmethod
    def Tnu0(self):
        ...

    @property
    @abstractmethod
    def has_massive_nu(self):
        ...

    @property
    @abstractmethod
    def Ogamma0(self):
        ...

    @property
    @abstractmethod
    def Onu0(self):
        ...

    # --- Methods ---

    @abstractmethod
    def w(self, z, /):
        ...

    @abstractmethod
    def nu_relative_density(self, z, /):
        ...

    @abstractmethod
    def de_density_scale(self, z, /):
        ...

    @abstractmethod
    def efunc(self, z, /):
        ...

    @abstractmethod
    def inv_efunc(self, z, /):
        ...

    @abstractmethod
    def lookback_time(self, z, /):
        ...

    @abstractmethod
    def age(self, z, /):
        ...

    @abstractmethod
    def comoving_distance(self, z, /):
        ...

    @abstractmethod
    def angular_diameter_distance_z1z2(self, z1, z2, /):
        ...

    @abstractmethod
    def absorption_distance(self, z, /):
        ...


@dataclass(frozen=True)
class AstropyFLRWCosmologyWrapper(FLRWCosmologyWrapperBase, AstropyCosmologyWrapper):
    """A wrapper for an astropy FLRW cosmology object.

    This is a wrapper for an astropy cosmology object that implements the
    FLRWAPI interface. It is used to allow the astropy cosmology objects to
    be used in the same way as the FLRWAPI objects.

    Parameters
    ----------
    cosmo : `astropy.cosmology.FLRW`
        The astropy cosmology object to wrap.

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> from lsst.afw.math import AstropyFLRWWrapper
    >>> cosmo = AstropyFLRWWrapper(Planck18)
    >>> cosmo.h
    .6766
    """

    cosmo: FLRW
    name = _GetFromCosmo()

    H0 = _GetFromCosmo()
    Om0 = _GetFromCosmo()
    Ode0 = _GetFromCosmo()
    Tcmb0 = _GetFromCosmo()
    Neff = _GetFromCosmo()
    m_nu = _GetFromCosmo()
    Ob0 = _GetFromCosmo()
    Odm0 = _GetFromCosmo()
    Ok0 = _GetFromCosmo()
    Tnu0 = _GetFromCosmo()
    has_massive_nu = _GetFromCosmo()
    hubble_time = _GetFromCosmo()
    hubble_distance = _GetFromCosmo()
    critical_density0 = _GetFromCosmo()
    Ogamma0 = _GetFromCosmo()
    Onu0 = _GetFromCosmo()

    w = _GetFromCosmo()
    Otot = _GetFromCosmo()
    Om = _GetFromCosmo()
    Ob = _GetFromCosmo()
    Odm = _GetFromCosmo()
    Ok = _GetFromCosmo()
    Ode = _GetFromCosmo()
    Ogamma = _GetFromCosmo()
    Onu = _GetFromCosmo()
    Tcmb = _GetFromCosmo()
    Tnu = _GetFromCosmo()
    nu_relative_density = _GetFromCosmo()
    de_density_scale = _GetFromCosmo()
    efunc = _GetFromCosmo()
    inv_efunc = _GetFromCosmo()
    H = _GetFromCosmo()
    scale_factor = _GetFromCosmo()
    lookback_time = _GetFromCosmo()
    lookback_distance = _GetFromCosmo()
    age = _GetFromCosmo()
    critical_density = _GetFromCosmo()
    comoving_distance = _GetFromCosmo()
    comoving_transverse_distance = _GetFromCosmo()
    angular_diameter_distance = _GetFromCosmo()
    luminosity_distance = _GetFromCosmo()
    angular_diameter_distance_z1z2 = _GetFromCosmo()
    absorption_distance = _GetFromCosmo()
    distmod = _GetFromCosmo()
    comoving_volume = _GetFromCosmo()
    differential_comoving_volume = _GetFromCosmo()
