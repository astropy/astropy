# This file includes the definition of a mix-in class that provides the low-
# and high-level WCS API to the astropy.wcs.WCS object. We keep this code
# isolated in this mix-in class to avoid making the main wcs.py file too
# long.

import warnings

import numpy as np

from astropy import units as u
from astropy.constants import c
from astropy.coordinates import ICRS, Galactic, SpectralCoord
from astropy.coordinates.spectral_coordinate import (
    attach_zero_velocities,
    update_differentials_to_match,
)
from astropy.utils.exceptions import AstropyUserWarning

from .high_level_api import HighLevelWCSMixin
from .low_level_api import BaseLowLevelWCS
from .wrappers import SlicedLowLevelWCS

__all__ = ["custom_ctype_to_ucd_mapping", "SlicedFITSWCS", "FITSWCSAPIMixin"]

C_SI = c.si.value

VELOCITY_FRAMES = {
    "GEOCENT": "gcrs",
    "BARYCENT": "icrs",
    "HELIOCENT": "hcrs",
    "LSRK": "lsrk",
    "LSRD": "lsrd",
}

# The spectra velocity frames below are needed for FITS spectral WCS
#  (see Greisen 06 table 12) but aren't yet defined as real
# astropy.coordinates frames, so we instead define them here as instances
# of existing coordinate frames with offset velocities. In future we should
# make these real frames so that users can more easily recognize these
# velocity frames when used in SpectralCoord.

# This frame is defined as a velocity of 220 km/s in the
# direction of l=90, b=0. The rotation velocity is defined
# in:
#
#   Kerr and Lynden-Bell 1986, Review of galactic constants.
#
# NOTE: this may differ from the assumptions of galcen_v_sun
# in the Galactocentric frame - the value used here is
# the one adopted by the WCS standard for spectral
# transformations.

VELOCITY_FRAMES["GALACTOC"] = Galactic(
    u=0 * u.km,
    v=0 * u.km,
    w=0 * u.km,
    U=0 * u.km / u.s,
    V=-220 * u.km / u.s,
    W=0 * u.km / u.s,
    representation_type="cartesian",
    differential_type="cartesian",
)

# This frame is defined as a velocity of 300 km/s in the
# direction of l=90, b=0. This is defined in:
#
#   Transactions of the IAU Vol. XVI B Proceedings of the
#   16th General Assembly, Reports of Meetings of Commissions:
#   Comptes Rendus Des Séances Des Commissions, Commission 28,
#   p201.
#
# Note that these values differ from those used by CASA
# (308 km/s towards l=105, b=-7) but we use the above values
# since these are the ones defined in Greisen et al (2006).

VELOCITY_FRAMES["LOCALGRP"] = Galactic(
    u=0 * u.km,
    v=0 * u.km,
    w=0 * u.km,
    U=0 * u.km / u.s,
    V=-300 * u.km / u.s,
    W=0 * u.km / u.s,
    representation_type="cartesian",
    differential_type="cartesian",
)

# This frame is defined as a velocity of 368 km/s in the
# direction of l=263.85, b=48.25. This is defined in:
#
#   Bennett et al. (2003), First-Year Wilkinson Microwave
#   Anisotropy Probe (WMAP) Observations: Preliminary Maps
#   and Basic Results
#
# Note that in that paper, the dipole is expressed as a
# temperature (T=3.346 +/- 0.017mK)

VELOCITY_FRAMES["CMBDIPOL"] = Galactic(
    l=263.85 * u.deg,
    b=48.25 * u.deg,
    distance=0 * u.km,
    radial_velocity=-(3.346e-3 / 2.725 * c).to(u.km / u.s),
)


# Mapping from CTYPE axis name to UCD1

CTYPE_TO_UCD1 = {
    # Celestial coordinates
    "RA": "pos.eq.ra",
    "DEC": "pos.eq.dec",
    "GLON": "pos.galactic.lon",
    "GLAT": "pos.galactic.lat",
    "ELON": "pos.ecliptic.lon",
    "ELAT": "pos.ecliptic.lat",
    "TLON": "pos.bodyrc.lon",
    "TLAT": "pos.bodyrc.lat",
    "HPLT": "custom:pos.helioprojective.lat",
    "HPLN": "custom:pos.helioprojective.lon",
    "HPRZ": "custom:pos.helioprojective.z",
    "HGLN": "custom:pos.heliographic.stonyhurst.lon",
    "HGLT": "custom:pos.heliographic.stonyhurst.lat",
    "CRLN": "custom:pos.heliographic.carrington.lon",
    "CRLT": "custom:pos.heliographic.carrington.lat",
    "SOLX": "custom:pos.heliocentric.x",
    "SOLY": "custom:pos.heliocentric.y",
    "SOLZ": "custom:pos.heliocentric.z",
    # Spectral coordinates (WCS paper 3)
    "FREQ": "em.freq",  # Frequency
    "ENER": "em.energy",  # Energy
    "WAVN": "em.wavenumber",  # Wavenumber
    "WAVE": "em.wl",  # Vacuum wavelength
    "VRAD": "spect.dopplerVeloc.radio",  # Radio velocity
    "VOPT": "spect.dopplerVeloc.opt",  # Optical velocity
    "ZOPT": "src.redshift",  # Redshift
    "AWAV": "em.wl",  # Air wavelength
    "VELO": "spect.dopplerVeloc",  # Apparent radial velocity
    "BETA": "custom:spect.doplerVeloc.beta",  # Beta factor (v/c)
    "STOKES": "phys.polarization.stokes",  # STOKES parameters
    # Time coordinates (https://www.aanda.org/articles/aa/pdf/2015/02/aa24653-14.pdf)
    "TIME": "time",
    "TAI": "time",
    "TT": "time",
    "TDT": "time",
    "ET": "time",
    "IAT": "time",
    "UT1": "time",
    "UTC": "time",
    "GMT": "time",
    "GPS": "time",
    "TCG": "time",
    "TCB": "time",
    "TDB": "time",
    "LOCAL": "time",
    # Distance coordinates
    "DIST": "pos.distance",
    "DSUN": "custom:pos.distance.sunToObserver"
    # UT() and TT() are handled separately in world_axis_physical_types
}

# Keep a list of additional custom mappings that have been registered. This
# is kept as a list in case nested context managers are used
CTYPE_TO_UCD1_CUSTOM = []


class custom_ctype_to_ucd_mapping:
    """
    A context manager that makes it possible to temporarily add new CTYPE to
    UCD1+ mapping used by :attr:`FITSWCSAPIMixin.world_axis_physical_types`.

    Parameters
    ----------
    mapping : dict
        A dictionary mapping a CTYPE value to a UCD1+ value

    Examples
    --------
    Consider a WCS with the following CTYPE::

        >>> from astropy.wcs import WCS
        >>> wcs = WCS(naxis=1)
        >>> wcs.wcs.ctype = ['SPAM']

    By default, :attr:`FITSWCSAPIMixin.world_axis_physical_types` returns `None`,
    but this can be overridden::

        >>> wcs.world_axis_physical_types
        [None]
        >>> with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
        ...     wcs.world_axis_physical_types
        ['food.spam']
    """

    def __init__(self, mapping):
        CTYPE_TO_UCD1_CUSTOM.insert(0, mapping)
        self.mapping = mapping

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        CTYPE_TO_UCD1_CUSTOM.remove(self.mapping)


class SlicedFITSWCS(SlicedLowLevelWCS, HighLevelWCSMixin):
    pass


class FITSWCSAPIMixin(BaseLowLevelWCS, HighLevelWCSMixin):
    """
    A mix-in class that is intended to be inherited by the
    :class:`~astropy.wcs.WCS` class and provides the low- and high-level WCS API.
    """

    @property
    def pixel_n_dim(self):
        return self.naxis

    @property
    def world_n_dim(self):
        return len(self.wcs.ctype)

    @property
    def array_shape(self):
        if self.pixel_shape is None:
            return None
        else:
            return self.pixel_shape[::-1]

    @array_shape.setter
    def array_shape(self, value):
        if value is None:
            self.pixel_shape = None
        else:
            self.pixel_shape = value[::-1]

    @property
    def pixel_shape(self):
        if self._naxis == [0, 0]:
            return None
        else:
            return tuple(self._naxis)

    @pixel_shape.setter
    def pixel_shape(self, value):
        if value is None:
            self._naxis = [0, 0]
        else:
            if len(value) != self.naxis:
                raise ValueError(
                    f"The number of data axes, {self.naxis}, does not equal the shape"
                    f" {len(value)}."
                )
            self._naxis = list(value)

    @property
    def pixel_bounds(self):
        return self._pixel_bounds

    @pixel_bounds.setter
    def pixel_bounds(self, value):
        if value is None:
            self._pixel_bounds = value
        else:
            if len(value) != self.naxis:
                raise ValueError(
                    "The number of data axes, "
                    f"{self.naxis}, does not equal the number of "
                    f"pixel bounds {len(value)}."
                )
            self._pixel_bounds = list(value)

    @property
    def world_axis_physical_types(self):
        types = []
        # TODO: need to support e.g. TT(TAI)
        for ctype in self.wcs.ctype:
            if ctype.upper().startswith(("UT(", "TT(")):
                types.append("time")
            else:
                ctype_name = ctype.split("-")[0]
                for custom_mapping in CTYPE_TO_UCD1_CUSTOM:
                    if ctype_name in custom_mapping:
                        types.append(custom_mapping[ctype_name])
                        break
                else:
                    types.append(CTYPE_TO_UCD1.get(ctype_name.upper(), None))
        return types

    @property
    def world_axis_units(self):
        units = []
        for unit in self.wcs.cunit:
            if unit is None:
                unit = ""
            elif isinstance(unit, u.Unit):
                unit = unit.to_string(format="vounit")
            else:
                try:
                    unit = u.Unit(unit).to_string(format="vounit")
                except u.UnitsError:
                    unit = ""
            units.append(unit)
        return units

    @property
    def world_axis_names(self):
        return list(self.wcs.cname)

    @property
    def axis_correlation_matrix(self):
        # If there are any distortions present, we assume that there may be
        # correlations between all axes. Maybe if some distortions only apply
        # to the image plane we can improve this?
        if self.has_distortion:
            return np.ones((self.world_n_dim, self.pixel_n_dim), dtype=bool)

        # Assuming linear world coordinates along each axis, the correlation
        # matrix would be given by whether or not the PC matrix is zero
        matrix = self.wcs.get_pc() != 0

        # We now need to check specifically for celestial coordinates since
        # these can assume correlations because of spherical distortions. For
        # each celestial coordinate we copy over the pixel dependencies from
        # the other celestial coordinates.
        celestial = (self.wcs.axis_types // 1000) % 10 == 2
        celestial_indices = np.nonzero(celestial)[0]
        for world1 in celestial_indices:
            for world2 in celestial_indices:
                if world1 != world2:
                    matrix[world1] |= matrix[world2]
                    matrix[world2] |= matrix[world1]

        return matrix

    def pixel_to_world_values(self, *pixel_arrays):
        world = self.all_pix2world(*pixel_arrays, 0)
        return world[0] if self.world_n_dim == 1 else tuple(world)

    def world_to_pixel_values(self, *world_arrays):
        # avoid circular import
        from astropy.wcs.wcs import NoConvergence

        try:
            pixel = self.all_world2pix(*world_arrays, 0)
        except NoConvergence as e:
            warnings.warn(str(e))
            # use best_solution contained in the exception and format the same
            # way as all_world2pix does (using _array_converter)
            pixel = self._array_converter(
                lambda *args: e.best_solution, "input", *world_arrays, 0
            )

        return pixel[0] if self.pixel_n_dim == 1 else tuple(pixel)

    @property
    def world_axis_object_components(self):
        return self._get_components_and_classes()[0]

    @property
    def world_axis_object_classes(self):
        return self._get_components_and_classes()[1]

    @property
    def serialized_classes(self):
        return False

    def _get_components_and_classes(self):
        # The aim of this function is to return whatever is needed for
        # world_axis_object_components and world_axis_object_classes. It's easier
        # to figure it out in one go and then return the values and let the
        # properties return part of it.

        # Since this method might get called quite a few times, we need to cache
        # it. We start off by defining a hash based on the attributes of the
        # WCS that matter here (we can't just use the WCS object as a hash since
        # it is mutable)
        wcs_hash = (
            self.naxis,
            list(self.wcs.ctype),
            list(self.wcs.cunit),
            self.wcs.radesys,
            self.wcs.specsys,
            self.wcs.equinox,
            self.wcs.dateobs,
            self.wcs.lng,
            self.wcs.lat,
        )

        # If the cache is present, we need to check that the 'hash' matches.
        if getattr(self, "_components_and_classes_cache", None) is not None:
            cache = self._components_and_classes_cache
            if cache[0] == wcs_hash:
                return cache[1]
            else:
                self._components_and_classes_cache = None

        # Avoid circular imports by importing here
        from astropy.coordinates import EarthLocation, SkyCoord, StokesCoord
        from astropy.time import Time, TimeDelta
        from astropy.time.formats import FITS_DEPRECATED_SCALES
        from astropy.wcs.utils import wcs_to_celestial_frame

        components = [None] * self.naxis
        classes = {}

        # Let's start off by checking whether the WCS has a pair of celestial
        # components

        if self.has_celestial:
            try:
                celestial_frame = wcs_to_celestial_frame(self)
            except ValueError:
                # Some WCSes, e.g. solar, can be recognized by WCSLIB as being
                # celestial but we don't necessarily have frames for them.
                celestial_frame = None
            else:
                kwargs = {}
                kwargs["frame"] = celestial_frame
                # Very occasionally (i.e. with TAB) wcs does not convert the units to degrees
                kwargs["unit"] = (
                    u.Unit(self.wcs.cunit[self.wcs.lng]),
                    u.Unit(self.wcs.cunit[self.wcs.lat]),
                )

                classes["celestial"] = (SkyCoord, (), kwargs)

                components[self.wcs.lng] = ("celestial", 0, "spherical.lon.degree")
                components[self.wcs.lat] = ("celestial", 1, "spherical.lat.degree")

        # Next, we check for spectral components

        if self.has_spectral:
            # Find index of spectral coordinate
            ispec = self.wcs.spec
            ctype = self.wcs.ctype[ispec][:4]
            ctype = ctype.upper()

            kwargs = {}

            # Determine observer location and velocity

            # TODO: determine how WCS standard would deal with observer on a
            # spacecraft far from earth. For now assume the obsgeo parameters,
            # if present, give the geocentric observer location.

            if np.isnan(self.wcs.obsgeo[0]):
                observer = None
            else:
                earth_location = EarthLocation(*self.wcs.obsgeo[:3], unit=u.m)

                # Get the time scale from TIMESYS or fall back to 'utc'
                tscale = self.wcs.timesys.lower() or "utc"

                if np.isnan(self.wcs.mjdavg):
                    obstime = Time(
                        self.wcs.mjdobs,
                        format="mjd",
                        scale=tscale,
                        location=earth_location,
                    )
                else:
                    obstime = Time(
                        self.wcs.mjdavg,
                        format="mjd",
                        scale=tscale,
                        location=earth_location,
                    )
                observer_location = SkyCoord(earth_location.get_itrs(obstime=obstime))

                if self.wcs.specsys in VELOCITY_FRAMES:
                    frame = VELOCITY_FRAMES[self.wcs.specsys]
                    observer = observer_location.transform_to(frame)
                    if isinstance(frame, str):
                        observer = attach_zero_velocities(observer)
                    else:
                        observer = update_differentials_to_match(
                            observer_location,
                            VELOCITY_FRAMES[self.wcs.specsys],
                            preserve_observer_frame=True,
                        )
                elif self.wcs.specsys == "TOPOCENT":
                    observer = attach_zero_velocities(observer_location)
                else:
                    raise NotImplementedError(
                        f"SPECSYS={self.wcs.specsys} not yet supported"
                    )

            # Determine target

            # This is tricker. In principle the target for each pixel is the
            # celestial coordinates of the pixel, but we then need to be very
            # careful about SSYSOBS which is tricky. For now, we set the
            # target using the reference celestial coordinate in the WCS (if
            # any).

            if self.has_celestial and celestial_frame is not None:
                # NOTE: celestial_frame was defined higher up

                # NOTE: we set the distance explicitly to avoid warnings in SpectralCoord

                target = SkyCoord(
                    self.wcs.crval[self.wcs.lng] * self.wcs.cunit[self.wcs.lng],
                    self.wcs.crval[self.wcs.lat] * self.wcs.cunit[self.wcs.lat],
                    frame=celestial_frame,
                    distance=1000 * u.kpc,
                )

                target = attach_zero_velocities(target)

            else:
                target = None

            # SpectralCoord does not work properly if either observer or target
            # are not convertible to ICRS, so if this is the case, we (for now)
            # drop the observer and target from the SpectralCoord and warn the
            # user.

            if observer is not None:
                try:
                    observer.transform_to(ICRS())
                except Exception:
                    warnings.warn(
                        "observer cannot be converted to ICRS, so will "
                        "not be set on SpectralCoord",
                        AstropyUserWarning,
                    )
                    observer = None

            if target is not None:
                try:
                    target.transform_to(ICRS())
                except Exception:
                    warnings.warn(
                        "target cannot be converted to ICRS, so will "
                        "not be set on SpectralCoord",
                        AstropyUserWarning,
                    )
                    target = None

            # NOTE: below we include Quantity in classes['spectral'] instead
            # of SpectralCoord - this is because we want to also be able to
            # accept plain quantities.

            if ctype == "ZOPT":

                def spectralcoord_from_redshift(redshift):
                    if isinstance(redshift, SpectralCoord):
                        return redshift
                    return SpectralCoord(
                        (redshift + 1) * self.wcs.restwav,
                        unit=u.m,
                        observer=observer,
                        target=target,
                    )

                def redshift_from_spectralcoord(spectralcoord):
                    # TODO: check target is consistent between WCS and SpectralCoord,
                    # if they are not the transformation doesn't make conceptual sense.
                    if (
                        observer is None
                        or spectralcoord.observer is None
                        or spectralcoord.target is None
                    ):
                        if observer is None:
                            msg = "No observer defined on WCS"
                        elif spectralcoord.observer is None:
                            msg = "No observer defined on SpectralCoord"
                        else:
                            msg = "No target defined on SpectralCoord"
                        warnings.warn(
                            f"{msg}, SpectralCoord "
                            "will be converted without any velocity "
                            "frame change",
                            AstropyUserWarning,
                        )
                        return spectralcoord.to_value(u.m) / self.wcs.restwav - 1.0
                    else:
                        return (
                            spectralcoord.with_observer_stationary_relative_to(
                                observer
                            ).to_value(u.m)
                            / self.wcs.restwav
                            - 1.0
                        )

                classes["spectral"] = (u.Quantity, (), {}, spectralcoord_from_redshift)
                components[self.wcs.spec] = ("spectral", 0, redshift_from_spectralcoord)

            elif ctype == "BETA":

                def spectralcoord_from_beta(beta):
                    if isinstance(beta, SpectralCoord):
                        return beta
                    return SpectralCoord(
                        beta * C_SI,
                        unit=u.m / u.s,
                        doppler_convention="relativistic",
                        doppler_rest=self.wcs.restwav * u.m,
                        observer=observer,
                        target=target,
                    )

                def beta_from_spectralcoord(spectralcoord):
                    # TODO: check target is consistent between WCS and SpectralCoord,
                    # if they are not the transformation doesn't make conceptual sense.
                    doppler_equiv = u.doppler_relativistic(self.wcs.restwav * u.m)
                    if (
                        observer is None
                        or spectralcoord.observer is None
                        or spectralcoord.target is None
                    ):
                        if observer is None:
                            msg = "No observer defined on WCS"
                        elif spectralcoord.observer is None:
                            msg = "No observer defined on SpectralCoord"
                        else:
                            msg = "No target defined on SpectralCoord"
                        warnings.warn(
                            f"{msg}, SpectralCoord "
                            "will be converted without any velocity "
                            "frame change",
                            AstropyUserWarning,
                        )
                        return spectralcoord.to_value(u.m / u.s, doppler_equiv) / C_SI
                    else:
                        return (
                            spectralcoord.with_observer_stationary_relative_to(
                                observer
                            ).to_value(u.m / u.s, doppler_equiv)
                            / C_SI
                        )

                classes["spectral"] = (u.Quantity, (), {}, spectralcoord_from_beta)
                components[self.wcs.spec] = ("spectral", 0, beta_from_spectralcoord)

            else:
                kwargs["unit"] = self.wcs.cunit[ispec]

                if self.wcs.restfrq > 0:
                    if ctype == "VELO":
                        kwargs["doppler_convention"] = "relativistic"
                        kwargs["doppler_rest"] = self.wcs.restfrq * u.Hz
                    elif ctype == "VRAD":
                        kwargs["doppler_convention"] = "radio"
                        kwargs["doppler_rest"] = self.wcs.restfrq * u.Hz
                    elif ctype == "VOPT":
                        kwargs["doppler_convention"] = "optical"
                        kwargs["doppler_rest"] = self.wcs.restwav * u.m

                def spectralcoord_from_value(value):
                    if isinstance(value, SpectralCoord):
                        return value
                    return SpectralCoord(
                        value, observer=observer, target=target, **kwargs
                    )

                def value_from_spectralcoord(spectralcoord):
                    # TODO: check target is consistent between WCS and SpectralCoord,
                    # if they are not the transformation doesn't make conceptual sense.
                    if (
                        observer is None
                        or spectralcoord.observer is None
                        or spectralcoord.target is None
                    ):
                        if observer is None:
                            msg = "No observer defined on WCS"
                        elif spectralcoord.observer is None:
                            msg = "No observer defined on SpectralCoord"
                        else:
                            msg = "No target defined on SpectralCoord"
                        warnings.warn(
                            f"{msg}, SpectralCoord "
                            "will be converted without any velocity "
                            "frame change",
                            AstropyUserWarning,
                        )
                        return spectralcoord.to_value(**kwargs)
                    else:
                        return spectralcoord.with_observer_stationary_relative_to(
                            observer
                        ).to_value(**kwargs)

                classes["spectral"] = (u.Quantity, (), {}, spectralcoord_from_value)
                components[self.wcs.spec] = ("spectral", 0, value_from_spectralcoord)

        # We can then make sure we correctly return Time objects where appropriate
        # (https://www.aanda.org/articles/aa/pdf/2015/02/aa24653-14.pdf)

        if "time" in self.world_axis_physical_types:
            multiple_time = self.world_axis_physical_types.count("time") > 1

            for i in range(self.naxis):
                if self.world_axis_physical_types[i] == "time":
                    if multiple_time:
                        name = f"time.{i}"
                    else:
                        name = "time"

                    # Initialize delta
                    reference_time_delta = None

                    # Extract time scale, and remove any algorithm code
                    scale = self.wcs.ctype[i].split("-")[0].lower()

                    if scale == "time":
                        if self.wcs.timesys:
                            scale = self.wcs.timesys.lower()
                        else:
                            scale = "utc"

                    # Drop sub-scales
                    if "(" in scale:
                        pos = scale.index("(")
                        scale, subscale = scale[:pos], scale[pos + 1 : -1]
                        warnings.warn(
                            "Dropping unsupported sub-scale "
                            f"{subscale.upper()} from scale {scale.upper()}",
                            UserWarning,
                        )

                    # TODO: consider having GPS as a scale in Time
                    # For now GPS is not a scale, we approximate this by TAI - 19s
                    if scale == "gps":
                        reference_time_delta = TimeDelta(19, format="sec")
                        scale = "tai"

                    elif scale.upper() in FITS_DEPRECATED_SCALES:
                        scale = FITS_DEPRECATED_SCALES[scale.upper()]

                    elif scale not in Time.SCALES:
                        raise ValueError(f"Unrecognized time CTYPE={self.wcs.ctype[i]}")

                    # Determine location
                    trefpos = self.wcs.trefpos.lower()

                    if trefpos.startswith("topocent"):
                        # Note that some headers use TOPOCENT instead of TOPOCENTER
                        if np.any(np.isnan(self.wcs.obsgeo[:3])):
                            warnings.warn(
                                "Missing or incomplete observer location "
                                "information, setting location in Time to None",
                                UserWarning,
                            )
                            location = None
                        else:
                            location = EarthLocation(*self.wcs.obsgeo[:3], unit=u.m)
                    elif trefpos == "geocenter":
                        location = EarthLocation(0, 0, 0, unit=u.m)
                    elif trefpos == "":
                        location = None
                    else:
                        # TODO: implement support for more locations when Time supports it
                        warnings.warn(
                            f"Observation location '{trefpos}' is not "
                            "supported, setting location in Time to None",
                            UserWarning,
                        )
                        location = None

                    reference_time = Time(
                        np.nan_to_num(self.wcs.mjdref[0]),
                        np.nan_to_num(self.wcs.mjdref[1]),
                        format="mjd",
                        scale=scale,
                        location=location,
                    )

                    if reference_time_delta is not None:
                        reference_time = reference_time + reference_time_delta

                    def time_from_reference_and_offset(offset):
                        if isinstance(offset, Time):
                            return offset
                        return reference_time + TimeDelta(offset, format="sec")

                    def offset_from_time_and_reference(time):
                        return (time - reference_time).sec

                    classes[name] = (Time, (), {}, time_from_reference_and_offset)
                    components[i] = (name, 0, offset_from_time_and_reference)

        if "phys.polarization.stokes" in self.world_axis_physical_types:
            for i in range(self.naxis):
                if self.world_axis_physical_types[i] == "phys.polarization.stokes":
                    name = "stokes"
                    classes[name] = (StokesCoord, (), {})
                    components[i] = (name, 0, "value")

        # Fallback: for any remaining components that haven't been identified, just
        # return Quantity as the class to use

        for i in range(self.naxis):
            if components[i] is None:
                name = self.wcs.ctype[i].split("-")[0].lower()
                if name == "":
                    name = "world"
                while name in classes:
                    name += "_"
                classes[name] = (u.Quantity, (), {"unit": self.wcs.cunit[i]})
                components[i] = (name, 0, "value")

        # Keep a cached version of result
        self._components_and_classes_cache = wcs_hash, (components, classes)

        return components, classes
