import warnings
from collections import namedtuple

import astropy.units as u
import numpy as np
from astropy.constants import c
from astropy.coordinates import (ICRS, BaseRepresentation,
                                 CartesianDifferential,
                                 CartesianRepresentation, SkyCoord)
from astropy.coordinates.baseframe import BaseCoordinateFrame, FrameMeta
from astropy.utils.exceptions import AstropyUserWarning

DOPPLER_CONVENTIONS = {
    'radio': u.doppler_radio,
    'optical': u.doppler_optical,
    'relativistic': u.doppler_relativistic
}

RV_RS_EQUIV = [(u.cm / u.s, u.Unit(''),
                lambda x: x / c.cgs.value,
                lambda x: x * c.cgs.value)]

DEFAULT_DISTANCE = 1 * u.AU

DopplerConversion = namedtuple('DopplerConversion', ['rest', 'convention'])

__all__ = ['SpectralCoord']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralCoord.*']


class SpectralCoord(u.Quantity):
    """
    Coordinate object representing spectral values.

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
                elif redshift is not None:
                    raise ValueError("Cannot set both a radial velocity and "
                                     "redshift on spectral coordinate.")

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
            elif redshift is not None:
                raise ValueError("Cannot set both a radial velocity and "
                                 "redshift on spectral coordinate.")

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

        # If the observer frame does not contain information about the
        # velocity of the system, assume that the velocity is zero in the
        # system.
        if 's' not in coord.data.differentials:
            warnings.warn(
                "No velocity defined on frame, assuming {}.".format(
                    u.Quantity([0, 0, 0], unit=u.km/u.s)),
                AstropyUserWarning)

            coord_diffs = CartesianDifferential(u.Quantity([0, 0, 0] * u.km / u.s))

            new_data = coord.data.to_cartesian().with_differentials(coord_diffs)
            coord = coord.realize_frame(new_data)

        return coord

    def _copy(self, **kwargs):
        default_kwargs = {
            'value': self.value,
            'unit': self.unit,
            'doppler_rest': self.doppler_rest,
            'doppler_convention': self.doppler_convention,
            'observer': self.observer,
            'target': self.target,
            'radial_velocity': self.radial_velocity
        }

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
        return np.sum(self._calculate_radial_velocity(
            self._observer, self._target), axis=0)

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
        return self.radial_velocity.to('', equivalencies=RV_RS_EQUIV)

    @staticmethod
    def _calculate_radial_velocity(observer, target):
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

        return d_vel.d_xyz * pos_hat.xyz

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

        # TODO: getting proper blueshifted results in the test involving
        #  shifting an observer to the target velocity frame requires that we
        #  negate the initial velocity between observer/target. Unsure why.
        new_data = self._project_velocity_and_shift(-init_obs_vel, fin_obs_vel)

        new_coord = self._copy(value=new_data,
                               observer=observer,
                               target=target)

        return new_coord

    def _project_velocity_and_shift(self, init_vel, fin_vel):
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
        # Project the velocity shift vector onto the line-on-sight vector
        #  between the target and the new observation frame.
        init_proj_vel = np.dot(init_vel, fin_vel / np.linalg.norm(fin_vel))

        # Calculate the velocity shift between the two vectors
        delta_vel = np.sum(fin_vel - init_proj_vel, axis=0)

        # In the case where the projected velocity is nan, we can assume that
        #  the final velocity different is zero, and thus the actual velocity
        #  delta is equal to the original radial velocity.

        # TODO: Due to lack of precision in some coordinate transformations,
        #  we may end up with a final velocity very close to, but not quite at,
        #  zero. In this case, set a tolerance for the final velocity; if it's
        #  below this tolerance, assume the delta velocity is essentially nan.
        if np.isnan(delta_vel) or np.abs(np.sum(fin_vel, axis=0)) < 1e-7 * fin_vel.unit:
            delta_vel = -self.radial_velocity

        # Apply the velocity shift to the stored spectral axis data
        new_data = self * (1 + delta_vel / c.cgs)

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

        if not frame.data.differentials:
            raise ValueError("Frame has no velocities, cannot transform "
                             "velocity frame.")

        observer_icrs = self.observer.transform_to(ICRS)
        frames_icrs = frame.transform_to(ICRS)

        data_with_rv = observer_icrs.data.with_differentials(
            frames_icrs.data.differentials)

        observer_icrs = observer_icrs.realize_frame(data_with_rv)

        observer = observer_icrs.transform_to(self.observer)

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

        new_data = self._project_velocity_and_shift(init_obs_vel, fin_obs_vel)

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

        return self._copy(value=rest_frame_value)

    def to_observed(self):
        """
        Transforms the spectral axis to the observed frame.
        """
        observed_frame_value = self * (1 + self.redshift)

        return self._copy(value=observed_frame_value)

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
