import warnings
from collections import namedtuple

import astropy.units as u
import numpy as np
from astropy.constants import c
from astropy.coordinates import (ICRS,
                                 CartesianDifferential,
                                 CartesianRepresentation, SkyCoord,
                                 Galactic, FK4, HCRS, GCRS)
from astropy.coordinates.baseframe import (BaseCoordinateFrame, FrameMeta,
                                           frame_transform_graph)
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.compat import NUMPY_LT_1_17

DOPPLER_CONVENTIONS = {
    'radio': u.doppler_radio,
    'optical': u.doppler_optical,
    'relativistic': u.doppler_relativistic
}

RV_RS_EQUIV = [(u.cm / u.s, u.Unit(''),
                lambda x: x / c.cgs.value,
                lambda x: x * c.cgs.value)]

DEFAULT_DISTANCE = 1 * u.AU

# TODO: disallow redshift if observer and target are specified

# FIXME: there are currently numerical issues when transforming frames with
# velocities when the position is exactly at the origin. To avoid this, we use
# a very small offset for now.
EPS = 1e-10

DopplerConversion = namedtuple('DopplerConversion', ['rest', 'convention'])

__all__ = ['SpectralCoord']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralCoord.*']


def update_differentials_to_match(original, velocity_reference):
    """
    Given an original coordinate object, update the differentials so that
    the final coordinate is at the same location as the original coordinate
    but co-moving with the velocity reference object.
    """

    if not velocity_reference.data.differentials:
        raise ValueError("Reference frame has no velocities")

    # If the reference has an obstime already defined, we should ignore
    # it and stick with the original observer obstime.
    if 'obstime' in velocity_reference.frame_attributes and hasattr(original, 'obstime'):
        velocity_reference = velocity_reference.replicate(obstime=original.obstime)

    # We transform both coordinates to ICRS for simplicity

    original_icrs = original.transform_to(ICRS())
    velocity_reference_icrs = velocity_reference.transform_to(ICRS())

    differentials = velocity_reference_icrs.data.represent_as(CartesianRepresentation,
                                                              CartesianDifferential).differentials
    if original_icrs.data.differentials:
        data_with_differentials = original_icrs.data.represent_as(CartesianRepresentation,
                                                                  CartesianDifferential).with_differentials(differentials)
    else:
        data_with_differentials = original_icrs.data.represent_as(CartesianRepresentation).with_differentials(differentials)
    final_icrs = original_icrs.realize_frame(data_with_differentials)
    final = final_icrs.transform_to(original)

    return final


def attach_zero_velocities(coord):
    """
    Set the differentials to be stationary on a coordinate object.
    """
    coord_diffs = CartesianDifferential(u.Quantity([0, 0, 0] * u.km / u.s))
    new_data = coord.data.to_cartesian().with_differentials(coord_diffs)
    return coord.realize_frame(new_data)


class SpectralCoord(u.Quantity):
    """
    Coordinate object representing spectral values.

    The `SpectralCoord` class is new in Astropy v4.1 and should be considered
    experimental at this time. It is possible that there will be API changes
    in future versions of Astropy based on user feedback. If you
    have specific ideas for how it might be improved, please  let us know on the
    `astropy-dev mailing list`_ or at http://feedback.astropy.org.

    Parameters
    ----------
    value : ndarray or `~astropy.units.Quantity` or `SpectralCoord`
        Spectral axis data values.
    unit : str or `~astropy.units.Unit`
        Unit for the given data.
    observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    target : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    radial_velocity : `~astropy.units.Quantity`, optional
        The radial velocity of the target with respect to the observer.
    redshift : float, optional
        The redshift of the target with respect to the observer.
    doppler_rest : `~astropy.units.Quantity`, optional
        The rest value to use for velocity space transformations.
    doppler_convention : str, optional
        The convention to use when converting the spectral data to/from
        velocity space.
    """
    _quantity_class = u.Quantity

    @u.quantity_input(doppler_rest=['length', 'frequency', None],
                      radial_velocity=['speed'])
    def __new__(cls, value, unit=None, observer=None, target=None,
                radial_velocity=None, redshift=None, doppler_rest=None,
                doppler_convention=None, **kwargs):
        obj = super().__new__(cls, value, unit=unit, subok=True, **kwargs)

        # Make sure incompatible inputs can't be specified at the same time

        if radial_velocity is not None and redshift is not None:
            raise ValueError("Cannot set both a radial velocity and "
                             "redshift on spectral coordinate.")

        if target is not None and observer is not None:
            if radial_velocity is not None:
                raise ValueError("Cannot specify radial velocity if both target "
                                 "and observer are specified")
            if redshift is not None:
                raise ValueError("Cannot specify radial velocity if both target "
                                 "and observer are specified")

        # The quantity machinery will drop the unit because type(value) !=
        #  SpectralCoord when passing in a Quantity object. Reassign the unit
        #  here to avoid this.
        if isinstance(value, u.Quantity) and unit is None:
            obj._unit = value.unit

        # If we're initializing from an existing SpectralCoord, keep any
        # parameters that aren't being overridden
        if isinstance(value, SpectralCoord):
            if observer is None:
                observer = value.observer
            if target is None:
                target = value.target
            if radial_velocity is None and redshift is None:
                radial_velocity = value.radial_velocity
            if doppler_rest is None:
                doppler_rest = value.doppler_rest
            if doppler_convention is None:
                doppler_convention = value.doppler_convention

        # Store state about whether the observer and target were defined
        #  explicitly (True), or implicity from rv/redshift (False)
        obj._frames_state = dict(observer=observer is not None,
                                 target=target is not None)

        obj._doppler_conversion = DopplerConversion(
            rest=doppler_rest, convention=doppler_convention)

        for x in [y for y in [observer, target] if y is not None]:
            if not isinstance(x, (SkyCoord, BaseCoordinateFrame)):
                raise ValueError("Observer must be a sky coordinate or "
                                 "coordinate frame.")

        # If no observer is defined, create a default observer centered in the
        #  ICRS frame.
        if observer is None:
            if target is None:
                observer = ICRS(ra=0 * u.degree, dec=0 * u.degree,
                                pm_ra_cosdec=0 * u.mas/u.yr, pm_dec=0 * u.mas/u.yr,
                                distance=0 * u.pc, radial_velocity=0 * u.km/u.s)
            else:
                if radial_velocity is None:
                    radial_velocity = 0 * u.km/u.s

                    if redshift is not None:
                        radial_velocity = u.Quantity(redshift).to(
                            'km/s', equivalencies=RV_RS_EQUIV)

                observer = SpectralCoord._target_from_observer(
                    target, -radial_velocity)

        # If no target is defined, create a default target with any provided
        #  redshift/radial velocities.
        if target is None:
            if radial_velocity is None:
                radial_velocity = 0 * u.km/u.s

                if redshift is not None:
                    radial_velocity = u.Quantity(redshift).to(
                        'km/s', equivalencies=RV_RS_EQUIV)

            target = SpectralCoord._target_from_observer(
                observer, radial_velocity)

        obj._observer = cls._validate_coordinate(observer)
        obj._target = cls._validate_coordinate(target)

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)

        self._frames_state = getattr(obj, '_frames_state', None)
        self._doppler_conversion = getattr(obj, '_doppler_conversion', None)

        self._observer = getattr(obj, '_observer', None)
        self._target = getattr(obj, '_target', None)

    def __quantity_subclass__(self, unit):
        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.
        """
        return SpectralCoord, True

    @staticmethod
    def _target_from_observer(observer, radial_velocity):
        """
        Generates a default target from a provided observer with an offset
        defined such as to create the provided radial velocity.

        Parameters
        ----------
        observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            Observer frame off which to base the target frame.
        radial_velocity : `~astropy.units.Quantity`
            Radial velocity used to calculate appropriate offsets between
            provided observer and generated target.

        Returns
        -------
        target : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            Generated target frame.
        """
        observer = SpectralCoord._validate_coordinate(observer)
        observer_icrs = observer.transform_to(ICRS)

        d = observer_icrs.cartesian.norm()
        drep = CartesianRepresentation([DEFAULT_DISTANCE.to(d.unit),
                                        0 * d.unit, 0 * d.unit])

        obs_vel = observer_icrs.cartesian.differentials['s']
        tot_rv = radial_velocity  # + observer_icrs.radial_velocity

        target = (observer_icrs.cartesian.without_differentials() + drep).with_differentials(
            CartesianDifferential([obs_vel.d_x + tot_rv,
                                   obs_vel.d_y.to(tot_rv.unit),
                                   obs_vel.d_z.to(tot_rv.unit)]))

        target = observer_icrs.realize_frame(target)

        return target

    @staticmethod
    def _validate_coordinate(coord):
        """
        Checks the type of the frame and whether a velocity differential and a
        distance has been defined on the frame object.

        If no distance is defined, the target is assumed to be "really far
        away", and the observer is assumed to be "in the solar system".

        Parameters
        ----------
        coord : `~astropy.coordinates.BaseCoordinateFrame`
            The new frame to be used for target or observer.
        """
        if not issubclass(coord.__class__, (BaseCoordinateFrame, FrameMeta)):
            if isinstance(coord, SkyCoord):
                coord = coord.frame
            else:
                raise ValueError("`{}` is not a subclass of "
                                 "`~astropy.coordinates.BaseCoordinateFrame` or "
                                 "`~astropy.coordinates.SkyCoord`.".format(coord))

        # If the distance is not well-defined, ensure that it works properly
        # for generating differentials
        # TODO: change this to not set the distance and yield a warning once there's a good way to address this in astropy.coordinates
        if hasattr(coord, 'distance') and \
                coord.distance.unit.physical_type == 'dimensionless':
            coord = SkyCoord(coord, distance=1e6 * u.kpc)
            warnings.warn(
                "Distance on coordinate object is dimensionless, an "
                "abritrary distance value of 1e6 kpc will be set instead.",
                AstropyUserWarning)

        # If the observer frame does not contain information about the
        # velocity of the system, assume that the velocity is zero in the
        # system.
        if 's' not in coord.data.differentials:
            warnings.warn(
                "No velocity defined on frame, assuming {}.".format(
                    u.Quantity([0, 0, 0], unit=u.km/u.s)),
                AstropyUserWarning)

            coord = attach_zero_velocities(coord)

        return coord

    def _copy(self, **kwargs):

        default_kwargs = {
            'value': self.value,
            'unit': self.unit,
            'doppler_rest': self.doppler_rest,
            'doppler_convention': self.doppler_convention,
            'observer': self.observer,
            'target': self.target,
        }

        # Only include radial_velocity if it is not auto-computed from the
        # observer and target.
        if self.observer is None or self.target is None:
            default_kwargs['radial_velocity'] = self.radial_velocity

        # If the new kwargs dict contains a value argument and it is a
        #  quantity, use the unit information provided by that quantity.
        # TODO: if a new quantity value is provided *and* the unit is set,
        #  do we implicitly try and convert to the provided unit?
        new_value = kwargs.get('value')

        if isinstance(new_value, u.Quantity):
            kwargs['unit'] = None

        default_kwargs.update(kwargs)

        return self.__class__(**default_kwargs)

    @property
    def quantity(self):
        """
        Convert the ``SpectralCoord`` to a `~astropy.units.Quantity`.
        Equivalent to ``self.view(u.Quantity)``.

        Returns
        -------
        `~astropy.units.Quantity`
            This object viewed as a `~astropy.units.Quantity`.

        """
        return self.view(u.Quantity)

    @property
    def observer(self):
        """
        The coordinate frame from which the observation was taken.

        Returns
        -------
        `~astropy.coordinates.BaseCoordinateFrame`
            The astropy coordinate frame representing the observation.
        """
        if self._frames_state['observer']:
            return self._observer

    @observer.setter
    def observer(self, value):
        if self.observer is not None:
            raise ValueError("Spectral coordinate already has a defined "
                             "observer.")

        self._frames_state['observer'] = value is not None

        value = self._validate_coordinate(value)

        # The default target is based off the observer frame. In the case
        #  where both observer/target are initialized to defaults, and then
        #  the user sets a new observer, we need to create a new target based
        #  on the input frame to maintain rv/redshift continuity.
        if self.target is None:
            self._target = self._target_from_observer(
                value, self.radial_velocity)

        self._observer = value

    @property
    def target(self):
        """
        The coordinate frame of the object being observed.

        Returns
        -------
        `~astropy.coordinates.BaseCoordinateFrame`
            The astropy coordinate frame representing the target.
        """
        if self._frames_state['target']:
            return self._target

    @target.setter
    def target(self, value):
        if self.target is not None:
            raise ValueError("Spectral coordinate already has a defined "
                             "target.")

        self._frames_state['target'] = value is not None

        value = self._validate_coordinate(value)

        self._target = value

    @property
    def doppler_rest(self):
        """
        The rest value of the spectrum used for transformations to/from
        velocity space.

        Returns
        -------
        `~astropy.units.Quantity`
            Rest value as an astropy `~astropy.units.Quantity` object.
        """
        return self._doppler_conversion.rest

    @doppler_rest.setter
    @u.quantity_input(value=['length', 'frequency', 'energy', 'speed', None])
    def doppler_rest(self, value):
        """
        New rest value needed for velocity-space conversions.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Rest value.
        """
        if self._doppler_conversion.rest is not None:
            raise ValueError("Doppler rest value has already been set. Use "
                             "the `to` method to update the stored value.")

        self._doppler_conversion = self._doppler_conversion._replace(
            rest=value)

    @property
    def doppler_convention(self):
        """
        The defined convention for conversions to/from velocity space.

        Returns
        -------
        str
            One of 'optical', 'radio', or 'relativistic' representing the
            equivalency used in the unit conversions.
        """
        return self._doppler_conversion.convention

    @doppler_convention.setter
    def doppler_convention(self, value):
        """
        New velocity convention used for velocity space conversions.

        Parameters
        ----------
        value

        Notes
        -----
        More information on the equations dictating the transformations can be
        found in the astropy documentation [1]_.

        References
        ----------
        .. [1] Astropy documentation: https://docs.astropy.org/en/stable/units/equivalencies.html#spectral-doppler-equivalencies

        """
        if value is not None and value not in DOPPLER_CONVENTIONS:
            raise ValueError("Unrecognized velocity convention: {}.".format(
                value))

        if self._doppler_conversion.convention is not None:
            raise ValueError("Doppler convention has already been set. Use "
                             "the `to` method to update the stored value.")

        self._doppler_conversion = self._doppler_conversion._replace(
            convention=value)

    @property
    def radial_velocity(self):
        """
        Radial velocity of target relative to the observer.

        Returns
        -------
        `~astropy.units.Quantity`
            Radial velocity of target.

        Notes
        -----
        This is different from the ``.radial_velocity`` property of a
        coordinate frame in that this calculates the radial velocity with
        respect to the *observer*, not the origin of the frame.
        """
        return self._calculate_radial_velocity(self._observer, self._target,
                                               as_scalar=True)

    @property
    def redshift(self):
        """
        Redshift of target relative to observer. Calculated from the radial
        velocity.

        Returns
        -------
        float
            Redshift of target.
        """
        try:
            return self.radial_velocity.to('', equivalencies=RV_RS_EQUIV)
        except Exception as exc:
            print(exc)
            raise

    @staticmethod
    def _calculate_radial_velocity(observer, target, as_scalar=False):
        """
        Compute the line-of-sight velocity from the observer to the target.

        Parameters
        ----------
        observer : `~astropy.coordinates.BaseCoordinateFrame`
            The frame of the observer.
        target : `~astropy.coordinates.BaseCoordinateFrame`
            The frame of the target.

        Returns
        -------
        `~astropy.units.Quantity`
            The radial velocity of the target with respect to the observer.
        """

        # Convert observer and target to ICRS to avoid finite differencing
        #  calculations that lack numerical precision.
        observer_icrs = observer.transform_to(ICRS)
        target_icrs = target.transform_to(ICRS)

        pos_hat = SpectralCoord._norm_d_pos(observer_icrs, target_icrs)

        d_vel = target_icrs.velocity - observer_icrs.velocity

        vel_mag = np.dot(d_vel.d_xyz, pos_hat.xyz)

        if NUMPY_LT_1_17:
            vel_mag *= d_vel.d_xyz.unit

        if as_scalar:
            return vel_mag
        else:
            return vel_mag * pos_hat.xyz

    @staticmethod
    def _norm_d_pos(observer, target):
        """
        Calculate the normalized position vector between two frames.

        Parameters
        ----------
        observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The observation frame or coordinate.
        target : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The target frame or coordinate.

        Returns
        -------
        pos_hat : `BaseRepresentation`
            Position representation.
        """
        d_pos = (target.data.without_differentials() -
                 observer.data.without_differentials()).to_cartesian()

        dp_norm = d_pos.norm()

        # Reset any that are 0 to 1 to avoid nans from 0/0
        dp_norm.ravel()[dp_norm.ravel() == 0] = 1 * dp_norm.unit

        pos_hat = d_pos / dp_norm

        return pos_hat

    def _change_observer_to(self, observer, target=None):
        """
        Moves the observer to the provided coordinate/frame.

        Parameters
        ----------
        observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The new observation frame or coordinate.
        target : `~astropy.coordinates.SkyCoord`, optional
            The `~astropy.coordinates.SkyCoord` object representing the target of the observation.
            If none given, defaults to currently defined target.

        Returns
        -------
        new_coord : `SpectralCoord`
            The new coordinate object representing the spectral data
            transformed to the new observer frame.
        """
        if self.observer is None:
            raise ValueError("No observer has been set, cannot change "
                             "observer.")

        target = self.target if target is None else target

        # Check velocities and distance values on the new frame. This is
        #  handled in the frame validation.
        observer = self._validate_coordinate(observer)

        # Calculate the initial and final los velocity
        init_obs_vel = self._calculate_radial_velocity(self.observer, target)
        fin_obs_vel = self._calculate_radial_velocity(observer, target)

        line_of_sight_unit_vec = self._norm_d_pos(observer, target).xyz

        new_data = self._project_velocity_and_shift(init_obs_vel, fin_obs_vel, line_of_sight_unit_vec)

        new_coord = self._copy(value=new_data,
                               observer=observer,
                               target=target)

        return new_coord

    def _project_velocity_and_shift(self, init_vel, fin_vel, line_of_sight_unit_vec):
        """
        Calculated the velocity projection given two vectors.

        Parameters
        ----------
        init_vel : `u.Quantity`
            Initial velocity vector.
        fin_vel : `u.Quantity`
            Final velocity vector.

        Returns
        -------
        new_data : `u.Quantity`
            Spectral axis data with velocity shift applied.
        """

        # Project the velocity shift vector onto the line-of-sight vector
        # between the target and the new observation frame.
        init_proj_vel = np.dot(init_vel, line_of_sight_unit_vec) * line_of_sight_unit_vec

        if NUMPY_LT_1_17:
            init_proj_vel *= init_vel.unit

        # Calculate the magnitude of the velocity shift between the two vectors.
        # The vectors are aligned but may be in opposite directions, so we use
        # the dot product to determine this.

        diff_vel = fin_vel - init_proj_vel
        delta_vel = np.dot(diff_vel, line_of_sight_unit_vec)

        if NUMPY_LT_1_17:
            delta_vel *= diff_vel.unit

        # In the case where the projected velocity is nan, we can assume that
        #  the final velocity different is zero, and thus the actual velocity
        #  delta is equal to the original radial velocity.

        # TODO: Due to lack of precision in some coordinate transformations,
        #  we may end up with a final velocity very close to, but not quite at,
        #  zero. In this case, set a tolerance for the final velocity; if it's
        #  below this tolerance, assume the delta velocity is essentially nan.

        fin_vel_mag = np.linalg.norm(fin_vel)

        if NUMPY_LT_1_17:
            fin_vel_mag *= fin_vel.unit

        if np.isnan(delta_vel) or fin_vel_mag < 1e-7 * fin_vel.unit:
            delta_vel = -self.radial_velocity

        if self.unit.is_equivalent(u.m):  # wavelength
            new_data = self * (1 + delta_vel / c.cgs)
        elif self.unit.is_equivalent(u.Hz) or self.unit.is_equivalent(u.eV):  # frequency or energy
            new_data = self / (1 + delta_vel / c.cgs)
        elif self.unit.is_equivalent(u.km / u.s):  # velocity
            new_data = self + delta_vel
        else:
            raise TypeError(f"Unexpected units in velocity shift: {self.unit}")

        return new_data

    def in_observer_velocity_frame(self, frame):
        """
        Alters the velocity frame of the observer, but not the position.

        Parameters
        ----------
        frame : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The observation frame containing the new velocity for the observer.

        Returns
        -------
        new_coord : `SpectralCoord`
            The new coordinate object representing the spectral data
            transformed based on the observer's new velocity frame.
        """
        if hasattr(frame, 'frame'):
            frame = frame.frame
        elif isinstance(frame, str):
            frame_cls = frame_transform_graph.lookup_name(frame)
            frame = frame_cls(0 * u.m, 0 * u.m, 0 * u.m,
                              0 * u.m / u.s, 0 * u.m / u.s, 0 * u.m / u.s,
                              representation_type='cartesian',
                              differential_type='cartesian')

        observer = update_differentials_to_match(self.observer, frame)

        new_coord = self._change_observer_to(observer)

        return new_coord

    def with_los_shift(self, target_shift=None, observer_shift=None):
        """
        Apply a velocity shift to this spectral coordinate. The shift
        can be provided as a redshift (float value) or radial velocity
        (quantity with physical type of 'speed').

        Parameters
        ----------
        target_shift : float or `~astropy.units.Quantity`
            Shift value to apply to current target.
        observer_shift : float or `~astropy.units.Quantity`
            Shift value to apply to current observer.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate with the target/observer velocity changed
            to incorporate the shift.
        """
        if observer_shift is not None and (self.target is None or
                                           self.observer is None):
            raise ValueError("Both an observer and target must be defined "
                             "before applying a velocity shift.")

        for arg in [x for x in [target_shift, observer_shift] if x is not None]:
            if isinstance(arg, u.Quantity) and \
                    arg.unit.physical_type not in ['speed', 'dimensionless']:
                raise u.UnitsError("Argument must have unit physical type "
                                   "'speed' for radial velocty or "
                                   "'dimesionless' for redshift.")

        # The target or observer value is defined but is not a quantity object,
        #  assume it's a redshift float value and convert to velocity
        if isinstance(target_shift, (float, int)) or \
                isinstance(target_shift, u.Quantity) and \
                target_shift.unit.physical_type == 'dimensionless':
            target_shift = u.Quantity(target_shift).to(
                'km/s', equivalencies=RV_RS_EQUIV)

        if isinstance(observer_shift, (float, int)) or \
                isinstance(observer_shift, u.Quantity) and \
                observer_shift.unit.physical_type == 'dimensionless':
            observer_shift = u.Quantity(observer_shift).to(
                'km/s', equivalencies=RV_RS_EQUIV)

        target_icrs = self._target.transform_to(ICRS)
        observer_icrs = self._observer.transform_to(ICRS)

        target_shift = 0 * u.km / u.s if target_shift is None else target_shift
        observer_shift = 0 * u.km / u.s if observer_shift is None else observer_shift

        pos_hat = SpectralCoord._norm_d_pos(observer_icrs, target_icrs)

        target_velocity = target_icrs.velocity + target_shift * pos_hat
        observer_velocity = observer_icrs.velocity + observer_shift * pos_hat

        new_target = target_icrs.realize_frame(
            target_icrs.cartesian.with_differentials(
                CartesianDifferential(target_velocity.xyz))).transform_to(
            self._target)

        new_observer = observer_icrs.realize_frame(
            observer_icrs.cartesian.with_differentials(
                CartesianDifferential(observer_velocity.xyz))).transform_to(
            self._observer)

        init_obs_vel = self._calculate_radial_velocity(observer_icrs, target_icrs)
        fin_obs_vel = self._calculate_radial_velocity(new_observer, new_target)

        line_of_sight_unit_vec = self._norm_d_pos(observer_icrs, target_icrs).xyz

        new_data = self._project_velocity_and_shift(init_obs_vel, fin_obs_vel, line_of_sight_unit_vec)

        # If an observer/target pair were not defined already, we want to avoid
        #  providing an explicit pair, so create a new spectral coord object
        #  with just the radial velocity set so new implicit observer/target
        #  will be created. Otherwise, use the available observer/target
        #  instances to create the pair.
        return self._copy(value=new_data,
                          observer=new_observer if self.observer is not None else None,
                          target=new_target if self.target is not None else None,
                          radial_velocity=np.sum(fin_obs_vel, axis=0) if self.observer is None and self.target is None else None)

    @u.quantity_input(rv=['speed'])
    def with_radial_velocity(self, rv):
        """
        Creates a new `SpectralCoord` object with the updated radial
        velocity value.

        Parameters
        ----------
        rv : `~astropy.units.Quantity`
            New radial velocity to a store in the `SpectralCoord` object.

        Returns
        -------
        `SpectralCoord`
            A new instance with the updated radial velocity value.
        """
        if self.observer is not None and self.target is not None:
            raise ValueError("Radial velocity cannot be set explicitly when "
                             "providing both an observer and target.")

        return self._copy(radial_velocity=rv, redshift=None)

    def with_redshift(self, rs):
        """
        Creates a new `SpectralCoord` object with the updated redshift value.

        Parameters
        ----------
        rs : float
            New redshift to a store in the `SpectralCoord` object.

        Returns
        -------
        `SpectralCoord`
            A new instance with the updated redshift value.
        """
        if self.observer is not None and self.target is not None:
            raise ValueError("Redshift cannot be set explicitly when "
                             "providing both an observer and target.")

        return self._copy(redshift=rs, radial_velocity=None)

    def to_rest(self):
        """
        Transforms the spectral axis to the rest frame.
        """
        rest_frame_value = self / (1 + self.redshift)

        return self._copy(value=rest_frame_value, radial_velocity=0 * u.km / u.s)

    def to(self, unit, equivalencies=[], doppler_rest=None,
           doppler_convention=None):
        """
        Overloaded parent ``to`` method to provide parameters for defining
        rest value and pre-defined conventions for unit transformations.

        Parameters
        ----------
        doppler_rest : `~astropy.units.Quantity`, optional
            The rest value used in the velocity space conversions. Providing
            the value here will set the value stored on the `SpectralCoord`
            instance.
        doppler_convention : {'relativistic', 'optical', 'radio'}, optional
            The velocity convention to use during conversions. Providing the
            value here will set the value stored on the `SpectralCoord`
            instance.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate object with data converted to the new unit.
        """
        if doppler_rest is not None:
            self._doppler_conversion = self._doppler_conversion._replace(
                rest=doppler_rest)

        if doppler_convention is not None:
            if doppler_convention not in DOPPLER_CONVENTIONS:
                raise ValueError(
                    "Unrecognized doppler convention: {}.".format(
                        doppler_convention))

            self._doppler_conversion = self._doppler_conversion._replace(
                convention=doppler_convention)

        equivs = u.spectral()

        if self.doppler_rest is not None and self.doppler_convention is not None:
            vel_equiv = DOPPLER_CONVENTIONS[self.doppler_convention](self.doppler_rest)
            equivs += vel_equiv

        # Compose the equivalencies for spectral conversions including the
        # appropriate velocity handling.
        equivalencies += equivs

        return super().to(unit, equivalencies=equivalencies)

    def to_value(self, *args, **kwargs):
        return self.to(*args, **kwargs).value

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        sep = ', '
        arrstr = np.array2string(self.view(np.ndarray), separator=sep,
                                 prefix=prefixstr)
        obs_frame = self.observer.__class__.__name__ \
            if self.observer is not None else 'None'
        tar_frame = self.target.__class__.__name__ \
            if self.target is not None else 'None'

        try:
            radial_velocity = self.radial_velocity
            redshift = self.redshift
        except ValueError:
            radial_velocity = redshift = 'Undefined'

        return f'{prefixstr}{arrstr}{self._unitstr:s}, \n' \
            f'\tradial_velocity={radial_velocity}, \n' \
            f'\tredshift={redshift}, \n' \
            f'\tdoppler_rest={self.doppler_rest}, \n' \
            f'\tdoppler_convention={self.doppler_convention}, \n' \
            f'\tobserver={obs_frame}, \n' \
            f'\ttarget={tar_frame}>'

    # Now define a set of standard velocity frames for use with the
    # in_observer_velocity_frame method. Note that for frames that define a
    # non-zero velocity, the velocities given are the reverse of the actual
    # velocity we want for the frame - for example to get a velocity frame
    # moving towards a star, we need to set the coordinate object to have a
    # negative velocity so that the star is moving towards the frame origin.

    # First off, we consider velocity frames that are stationaly relative
    # to already defined 3-d celestial frames.

    GEOCENTRIC = GCRS(x=EPS * u.km, y=EPS * u.km, z=EPS * u.km,
                      v_x=0 * u.km / u.s, v_y=0 * u.km / u.s, v_z=0 * u.km / u.s,
                      representation_type='cartesian',
                      differential_type='cartesian')
    """
    Geocentric velocity frame (stationary relative to GCRS origin).
    """

    BARYCENTRIC = ICRS(x=EPS * u.km, y=EPS * u.km, z=EPS * u.km,
                       v_x=0 * u.km / u.s, v_y=0 * u.km / u.s, v_z=0 * u.km / u.s,
                       representation_type='cartesian',
                       differential_type='cartesian')
    """
    Barycentric velocity frame (stationary relative to ICRS origin).
    """

    HELIOCENTRIC = HCRS(x=EPS * u.km, y=EPS * u.km, z=EPS * u.km,
                        v_x=0 * u.km / u.s, v_y=0 * u.km / u.s, v_z=0 * u.km / u.s,
                        representation_type='cartesian',
                        differential_type='cartesian')
    """
    Heliocentric velocity frame (stationary relative to HCRS origin).
    """

    # The LSRK velocity frame is defined as having a velocity
    # of 20 km/s towards RA=270 Dec=30 (B1900) relative to the
    # solar system Barycenter. This is defined in:
    #
    #   Gordon 1975, Methods of Experimental Physics: Volume 12:
    #   Astrophysics, Part C: Radio Observations - Section 6.1.5.
    #
    # We use the astropy.coordinates FK4 class here since it is
    # defined with the solar system barycenter as the origin.

    LSRK_GORDON1975 = FK4(x=EPS * u.km, y=EPS * u.km, z=EPS * u.km,
                          v_x=-20 * u.km / u.s * np.cos(270 * u.deg) * np.cos(30 * u.deg),
                          v_y=-20 * u.km / u.s * np.sin(270 * u.deg) * np.cos(30 * u.deg),
                          v_z=-20 * u.km / u.s * np.sin(30 * u.deg),
                          representation_type='cartesian',
                          differential_type='cartesian',
                          equinox='B1900')
    """
    Kinematic Local Standard of Rest (as defined by Gordon 1975).
    """

    # The LSRD velocity frame is defined as a velocity of
    # U=9 km/s, V=12 km/s, and W=7 km/s or 16.552945 km/s
    # towards l=53.13 b=25.02. This is defined in:
    #
    #   Delhaye 1965, Solar Motion and Velocity Distribution of
    #   Common Stars.

    LSRD_DELHAYE1965 = Galactic(u=EPS * u.km, v=EPS * u.km, w=EPS * u.km,
                                U=-9 * u.km / u.s, V=-12 * u.km / u.s, W=-7 * u.km / u.s,
                                representation_type='cartesian',
                                differential_type='cartesian')
    """
    Kinematic Local Standard of Rest (as defined by Delhaye 1965).
    """

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

    GALACTOCENTRIC_KLB1986 = Galactic(u=EPS * u.km, v=EPS * u.km, w=EPS * u.km,
                                      U=0 * u.km / u.s, V=-220 * u.km / u.s, W=0 * u.km / u.s,
                                      representation_type='cartesian',
                                      differential_type='cartesian')
    """
    Galactocentric velocity frame (as defined by Kerr and Lynden-Bell 1986).
    """

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

    LOCALGROUP_IAU1976 = Galactic(u=EPS * u.km, v=EPS * u.km, w=EPS * u.km,
                                  U=0 * u.km / u.s, V=-300 * u.km / u.s, W=0 * u.km / u.s,
                                  representation_type='cartesian',
                                  differential_type='cartesian')
    """
    Local group velocity frame (as defined in the IAU Proceedings of the 16th General Assembly).
    """

    # This frame is defined as a velocity of 368 km/s in the
    # direction of l=263.85, b=48.25. This is defined in:
    #
    #   Bennett et al. (2003), First-Year Wilkinson Microwave
    #   Anisotropy Probe (WMAP) Observations: Preliminary Maps
    #   and Basic Results
    #
    # Note that in that paper, the dipole is expressed as a
    # temperature (T=3.346 +/- 0.017mK)

    CMBDIPOL_WMAP1 = Galactic(l=263.85 * u.deg, b=48.25 * u.deg, distance=EPS * u.km,
                              radial_velocity=-(3.346e-3 / 2.725 * c).to(u.km/u.s))
    """
    CMB Dipole velocity frame (using parameters from WMAP1 in Bennett et al., 2003).
    """
