from __future__ import annotations

from abc import ABCMeta, abstractmethod
from dataclasses import dataclass

from .core import Cosmology, CosmologyAPI


@dataclass(frozen=True)
class CosmologyWrapperBase(CosmologyAPI, metaclass=ABCMeta):

    cosmo: CosmologyAPI

    def __getattr__(self, name):
        return getattr(self.cosmo, name)

    @property
    @abstractmethod
    def meta(self):
        ...

    @property
    @abstractmethod
    def name(self):
        ...

    @property
    @abstractmethod
    def is_flat(self):
        ...

    @property
    @abstractmethod
    def clone(self, *, meta=None, **kwargs):
        ...


class _GetFromCosmo(property):
    def __set_name__(self, objcls, name):
        self.name = name

    def __get__(self, obj, objcls=None):
        return getattr(obj.cosmo, self.name) if obj is not None else self


@dataclass(frozen=True)
class AstropyCosmologyWrapper(CosmologyWrapperBase):
    """A wrapper for an astropy Cosmology object.

    This is a wrapper for an astropy cosmology object that implements the
    CosmologyAPI interface. It is used to allow the astropy cosmology
    objects to be used in the same way as the CosmologyAPI objects.

    Note that you probably want to use the ``AstropyFLRWWrapper`` class
    instead of this one.

    Parameters
    ----------
    cosmo : `astropy.cosmology.Cosmology`
        The astropy cosmology object to wrap.

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> from lsst.afw.math import AstropyCosmologyWrapper
    >>> cosmo = AstropyCosmologyWrapper(Planck18)
    >>> cosmo.h
    .6766
    """

    cosmo: Cosmology

    meta = _GetFromCosmo()
    name = _GetFromCosmo()
    is_flat = _GetFromCosmo()

    def clone(self, *, meta=None, **kwargs):
        return type(self)(self.cosmo.clone(meta=meta, **kwargs))
