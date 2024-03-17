# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the classes and utility functions for distance and
cartesian coordinates.
"""

import warnings

import numpy as np

from astropy import units as u
from astropy.utils.compat import COPY_IF_NEEDED
from astropy.utils.exceptions import AstropyWarning

from .angles import Angle

__all__ = ["Distance"]


__doctest_requires__ = {"*": ["scipy"]}


class Distance(u.SpecificTypeQuantity):
    """
    A one-dimensional distance.

    This can be initialized by providing one of the following:

    * Distance ``value`` (array or float) and a ``unit``
    * |Quantity| object with dimensionality of length
    * Redshift and (optionally) a `~astropy.cosmology.Cosmology`
    * Distance modulus
    * Parallax

    Parameters
    ----------
    value : scalar or `~astropy.units.Quantity` ['length']
        The value of this distance.
    unit : `~astropy.units.UnitBase` ['length']
        The unit for this distance.
    z : float
        A redshift for this distance.  It will be converted to a distance
        by computing the luminosity distance for this redshift given the
        cosmology specified by ``cosmology``. Must be given as a keyword
        argument.
    cosmology : `~astropy.cosmology.Cosmology` or None
        A cosmology that will be used to compute the distance from ``z``.
        If `None`, the current cosmology will be used (see
        `astropy.cosmology` for details).
    distmod : float or `~astropy.units.Quantity`
        The distance modulus for this distance. Note that if ``unit`` is not
        provided, a guess will be made at the unit between AU, pc, kpc, and Mpc.
    parallax : angle-like
        The parallax in angular units.
    dtype : `~numpy.dtype`, optional
        See `~astropy.units.Quantity`.
    copy : bool, optional
        See `~astropy.units.Quantity`.
    order : {'C', 'F', 'A'}, optional
        See `~astropy.units.Quantity`.
    subok : bool, optional
        See `~astropy.units.Quantity`.
    ndmin : int, optional
        See `~astropy.units.Quantity`.
    allow_negative : bool, optional
        Whether to allow negative distances (which are possible in some
        cosmologies).  Default: `False`.

    Raises
    ------
    `~astropy.units.UnitsError`
        If the ``unit`` is not a length unit.
    ValueError
        If value specified is less than 0 and ``allow_negative=False``.

        If ``cosmology`` is provided when ``z`` is *not* given.

        If either none or more than one of ``value``, ``z``, ``distmod``,
        or ``parallax`` were given.


    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.cosmology import WMAP5
    >>> Distance(10, u.Mpc)
    <Distance 10. Mpc>
    >>> Distance(40*u.pc, unit=u.kpc)
    <Distance 0.04 kpc>
    >>> Distance(z=0.23)                      # doctest: +FLOAT_CMP
    <Distance 1184.01657566 Mpc>
    >>> Distance(z=0.23, cosmology=WMAP5)     # doctest: +FLOAT_CMP
    <Distance 1147.78831918 Mpc>
    >>> Distance(distmod=24.47*u.mag)         # doctest: +FLOAT_CMP
    <Distance 783.42964277 kpc>
    >>> Distance(parallax=21.34*u.mas)        # doctest: +FLOAT_CMP
    <Distance 46.86035614 pc>
    """

    _equivalent_unit = u.m
    _include_easy_conversion_members = True

    def __new__(
        cls,
        value=None,
        unit=None,
        z=None,
        cosmology=None,
        distmod=None,
        parallax=None,
        dtype=np.inexact,
        copy=True,
        order=None,
        subok=False,
        ndmin=0,
        allow_negative=False,
    ):
        n_not_none = sum(x is not None for x in [value, z, distmod, parallax])
        if n_not_none == 0:
            raise ValueError(
                "none of `value`, `z`, `distmod`, or `parallax` "
                "were given to Distance constructor"
            )
        if n_not_none > 1:
            raise ValueError(
                "more than one of `value`, `z`, `distmod`, or "
                "`parallax` were given to Distance constructor"
            )

        if value is None:
            # If something else but `value` was provided then a new array will
            # be created anyways and there is no need to copy that.
            copy = False

        if z is not None:
            if cosmology is None:
                from astropy.cosmology import default_cosmology

                cosmology = default_cosmology.get()

            value = cosmology.luminosity_distance(z)

        elif cosmology is not None:
            raise ValueError(
                "a `cosmology` was given but `z` was not "
                "provided in Distance constructor"
            )

        elif distmod is not None:
            value = cls._distmod_to_pc(distmod)
            if unit is None:
                # if the output unit is not specified, convert `value`
                # based on the mean of the log of the distance.
                # Leaving `unit=None` is fine for the `super().__new__()` call below.
                meanlogval = np.log10(value.value).mean()
                if meanlogval > 6:
                    value <<= u.Mpc
                elif meanlogval > 3:
                    value <<= u.kpc
                elif meanlogval < -3:  # ~200 AU
                    value <<= u.AU

        elif parallax is not None:
            parallax = u.Quantity(parallax, copy=COPY_IF_NEEDED, subok=True)
            value = parallax.to(unit or u.pc, equivalencies=u.parallax())

            if np.any(parallax < 0):
                if not allow_negative:
                    raise ValueError(
                        "some parallaxes are negative, which are not "
                        "interpretable as distances. See the discussion in "
                        "this paper: https://arxiv.org/abs/1507.02105 . You "
                        "can convert negative parallaxes to NaN distances by "
                        "providing the `allow_negative=True` argument."
                    )
                warnings.warn(
                    "negative parallaxes are converted to NaN distances even when"
                    " `allow_negative=True`, because negative parallaxes cannot be"
                    " transformed into distances. See the discussion in this paper:"
                    " https://arxiv.org/abs/1507.02105",
                    AstropyWarning,
                )
            allow_negative = True  # No need to check twice.

        # now we have arguments like for a Quantity, so let it do the work
        distance = super().__new__(
            cls,
            value,
            unit,
            dtype=dtype,
            copy=copy,
            order=order,
            subok=subok,
            ndmin=ndmin,
        )

        if not allow_negative and np.any(distance.value < 0):
            raise ValueError(
                "distance must be >= 0. Use the argument "
                "`allow_negative=True` to allow negative values."
            )

        return distance

    @property
    def z(self):
        """Short for ``self.compute_z()``."""
        return self.compute_z()

    def compute_z(self, cosmology=None, **atzkw):
        """
        The redshift for this distance assuming its physical distance is
        a luminosity distance.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` or None
            The cosmology to assume for this calculation, or `None` to use the
            current cosmology (see `astropy.cosmology` for details).
        **atzkw
            keyword arguments for :func:`~astropy.cosmology.z_at_value`

        Returns
        -------
        z : `~astropy.units.Quantity`
            The redshift of this distance given the provided ``cosmology``.

        Warnings
        --------
        This method can be slow for large arrays.
        The redshift is determined using :func:`astropy.cosmology.z_at_value`,
        which handles vector inputs (e.g. an array of distances) by
        element-wise calling of :func:`scipy.optimize.minimize_scalar`.
        For faster results consider using an interpolation table;
        :func:`astropy.cosmology.z_at_value` provides details.

        See Also
        --------
        :func:`astropy.cosmology.z_at_value` : Find the redshift corresponding to a
            :meth:`astropy.cosmology.FLRW.luminosity_distance`.
        """
        from astropy.cosmology import z_at_value

        if cosmology is None:
            from astropy.cosmology import default_cosmology

            cosmology = default_cosmology.get()

        atzkw.setdefault("ztol", 1.0e-10)
        return z_at_value(cosmology.luminosity_distance, self, **atzkw)

    @property
    def distmod(self):
        """The distance modulus as a `~astropy.units.Quantity`."""
        val = 5.0 * np.log10(self.to_value(u.pc)) - 5.0
        return u.Quantity(val, u.mag, copy=COPY_IF_NEEDED)

    @classmethod
    def _distmod_to_pc(cls, dm):
        dm = u.Quantity(dm, u.mag)
        return cls(10 ** ((dm.value + 5) / 5.0), u.pc, copy=COPY_IF_NEEDED)

    @property
    def parallax(self):
        """The parallax angle as an `~astropy.coordinates.Angle` object."""
        return Angle(self.to(u.milliarcsecond, u.parallax()))
