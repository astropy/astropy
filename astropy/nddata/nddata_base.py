# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDDataBase class.


from abc import ABCMeta, abstractmethod

__all__ = ["NDDataBase"]


class NDDataBase(metaclass=ABCMeta):
    """Base metaclass that defines the interface for N-dimensional datasets
    with associated meta information used in ``astropy``.

    All properties and ``__init__`` have to be overridden in subclasses. See
    `NDData` for a subclass that defines this interface on `numpy.ndarray`-like
    ``data``.

    See also: https://docs.astropy.org/en/stable/nddata/

    """

    @abstractmethod
    def __init__(self):
        pass

    @property
    @abstractmethod
    def data(self):
        """The stored dataset."""

    @property
    @abstractmethod
    def mask(self):
        """Mask for the dataset.

        Masks should follow the ``numpy`` convention that **valid** data points
        are marked by ``False`` and **invalid** ones with ``True``.
        """
        return None

    @property
    @abstractmethod
    def unit(self):
        """Unit for the dataset."""
        return None

    @property
    @abstractmethod
    def wcs(self):
        """World coordinate system (WCS) for the dataset."""
        return None

    # psf is concrete to avoid introducing a breaking change in release 5.2.
    @property
    def psf(self):
        """Image representation of the PSF for the dataset.

        Should be `ndarray`-like.
        """
        return None

    @property
    @abstractmethod
    def meta(self):
        """Additional meta information about the dataset.

        Should be `dict`-like.
        """
        return None

    @property
    @abstractmethod
    def uncertainty(self):
        """Uncertainty in the dataset.

        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, such as ``"std"`` for standard deviation or
        ``"var"`` for variance.
        """
        return None
