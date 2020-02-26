import warnings

import numpy as np
from collections import namedtuple

import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord, ICRS, Distance, GCRS, RadialDifferential, CartesianDifferential
from astropy.coordinates.baseframe import BaseCoordinateFrame, FrameMeta
from astropy.utils.exceptions import AstropyUserWarning

DOPPLER_CONVENTIONS = {
    'radio': u.doppler_radio,
    'optical': u.doppler_optical,
    'relativistic': u.doppler_relativistic
}

RV_RS_EQUIV = [(u.cm / u.s, u.Unit(''),
                lambda x: np.sqrt((1 + x/c.cgs.value)/(1 - x/c.cgs.value)) - 1,
                lambda x: ((x + 1) ** 2 - 1) / ((x + 1) ** 2 + 1) * c.cgs.value)]

DopplerConversion = namedtuple('DopplerConversion', ['rest', 'convention'])

__all__ = ['SpectralCoord']


class SpectralCoord(u.Quantity):
    """
    Coordinate object representing spectral values.

    Attributes
    ----------
    value : ndarray or `Quantity` or `SpectralCoord`
        Spectral axis data values.
    unit : str or `Unit`
        Unit for the given data.
    observer : `BaseCoordinateFrame` or `SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    target : `BaseCoordinateFrame` or `SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    doppler_rest : `Quantity`, optional
        The rest value to use for velocity space transformations.
    doppler_convention : str, optional
        The convention to use when converting the spectral data to/from
        velocity space.
    """
    _quantity_class = u.Quantity

    @u.quantity_input(doppler_rest=['length', 'frequency', None])
    def __new__(cls, value, unit=None, observer=None, target=None,
                doppler_rest=None, doppler_convention=None, **kwargs):
        obj = super().__new__(cls, value, unit=unit, subok=True, **kwargs)

        # The quantity machinery will drop the unit because type(value) !=
        #  SpectralCoord when passing in a Quantity object. Reassign the unit
        #  here to avoid this.
        if isinstance(value, u.Quantity) and unit is None:
            obj._unit = value.unit

        obj._doppler_conversion = DopplerConversion(
            rest=doppler_rest, convention=doppler_convention)

        for x in [y for y in [observer, target] if y is not None]:
            if not isinstance(x, (SkyCoord, BaseCoordinateFrame)):
                raise ValueError("Observer must be a sky coordinate or "
                                 "coordinate frame.")

        obj._observer = cls._validate_coordinate(observer) if observer is not None else None
        obj._target = cls._validate_coordinate(target) if target is not None else None

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)

        self._doppler_conversion = getattr(obj, '_doppler_conversion', None)

        self._observer = getattr(obj, 'observer', None)
        self._target = getattr(obj, 'target', None)

    def __quantity_subclass__(self, unit):
        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.
        """
        return SpectralCoord, True

    @staticmethod
    def _validate_coordinate(coord):
        """
        Checks the type of the frame and whether a velocity differential and a
        distance has been defined on the frame object.

        If no distance is defined, the target is assumed to be "really far
        away", and the observer is assumed to be "in the solar system".

        Parameters
        ----------
        coord : `BaseCoordinateFrame`
            The new frame to be used for target or observer.
        """
        if not issubclass(coord.__class__, (BaseCoordinateFrame, FrameMeta)):
            if isinstance(coord, SkyCoord):
                coord = coord.frame
            else:
                raise ValueError("`{}` is not a subclass of "
                                 "`BaseCoordinateFrame` or "
                                 "`SkyCoord`.".format(coord))

        # If the observer frame does not contain information about the
        # velocity of the system, assume that the velocity is zero in the
        # system.
        if 's' not in coord.data.differentials:
            warnings.warn(
                "No velocity defined on frame, assuming {}.".format(
                   u.Quantity([0, 0, 0], unit=u.km/u.s)),
                AstropyUserWarning)

            vel_to_add = CartesianDifferential(
                0 * u.km / u.s, 0 * u.km / u.s, 0 * u.km/u.s)
            new_data = coord.data.to_cartesian().with_differentials(vel_to_add)
            coord = coord.realize_frame(new_data)

        return coord

    def _copy(self, **kwargs):
        default_kwargs = {
            'value': self.value,
            'unit': self.unit,
            'doppler_rest': self.doppler_rest,
            'doppler_convention': self.doppler_convention,
            'observer': self.observer,
            'target': self.target
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
    def observer(self):
        """
        The coordinate frame from which the observation was taken.

        Returns
        -------
        `BaseCoordinateFrame`
            The astropy coordinate frame representing the observation.
        """
        return self._observer

    @observer.setter
    def observer(self, value):
        if self.observer is not None:
            raise ValueError("Spectral coordinate already has a defined "
                             "observer.")

        value = self._validate_coordinate(value)

        self._observer = value

    @property
    def target(self):
        """
        The coordinate frame of the object being observed.

        Returns
        -------
        `BaseCoordinateFrame`
            The astropy coordinate frame representing the target.
        """
        return self._target

    @target.setter
    def target(self, value):
        if self.target is not None:
            raise ValueError("Spectral coordinate already has a defined "
                             "target.")

        value = self._validate_coordinate(value)

        self._target = value

    @property
    def doppler_rest(self):
        """
        The rest value of the spectrum used for transformations to/from
        velocity space.

        Returns
        -------
        `Quantity`
            Rest value as an astropy `Quantity` object.
        """
        return self._doppler_conversion.rest

    @doppler_rest.setter
    @u.quantity_input(value=['length', 'frequency', 'energy', 'speed', None])
    def doppler_rest(self, value):
        """
        New rest value needed for velocity-space conversions.

        Parameters
        ----------
        value : `Quantity`
            Rest value.
        """
        if self._doppler_conversion.rest is not None:
            raise ValueError("Doppler rest value has already been set. Use "
                             "the `to` method to update the stored value.")

        self._doppler_conversion = self._doppler_conversion._replace(rest=value)

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

        self._doppler_conversion = self._doppler_conversion._replace(convention=value)

    @property
    def radial_velocity(self):
        """
        Radial velocity of target relative to the observer.

        Returns
        -------
        `u.Quantity`
            Radial velocity of target.

        Notes
        -----
        This is different from the ``.radial_velocity`` property of a
        coordinate frame in that this calculates the radial velocity with
        respect to the *observer*, not the origin of the frame.
        """
        if self.observer is None or self.target is None:
            return 0 * u.km / u.s

        return self._calculate_radial_velocity(self.observer, self.target)

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
        observer : `BaseCoordinateFrame`
            The frame of the observer.
        target : `BaseCoordinateFrame`
            The frame of the target.

        Returns
        -------
        `Quantity`
            The radial velocity of the target with respect to the observer.
        """
        # Convert observer and target to ICRS to avoid finite
        # differencing calculations that lack numerical precision.
        if observer is None:
            raise ValueError("No observer has been specified; cannot "
                             "calculate radial velocity.")

        observer_icrs = observer.transform_to(ICRS)

        if target is None:
            raise ValueError("No target has been specified; cannot calculate "
                             "radial velocity.")

        target_icrs = target.transform_to(ICRS)

        d_pos = (target_icrs.data.without_differentials() -
                 observer_icrs.data.without_differentials()).to_cartesian()

        pos_hat = d_pos / (d_pos.norm() or 1)
        d_vel = target_icrs.velocity - observer_icrs.velocity

        return np.sum(d_vel.d_xyz * pos_hat.xyz, axis=0)

    def _change_observer_to(self, observer, target=None):
        """
        Moves the observer to the provided coordinate/frame.

        Parameters
        ----------
        observer : `BaseCoordinateFrame` or `SkyCoord`
            The new observation frame or coordinate.
        target : `SkyCoord`, optional
            The `SkyCoord` object representing
            the target of the observation.

        Returns
        -------
        `SpectralCoord`
            The new coordinate object representing the spectral data
            transformed to the new observer frame.
        """
        if self.observer is None:
            raise ValueError("No observer has been set, cannot transform "
                             "observer frame.")

        target = self.target if target is None else target

        # Check velocities and distance values on the new frame. This is
        #  handled in the frame validation.
        observer = self._validate_coordinate(observer)

        # Calculate the initial and final los velocity
        init_obs_vel = self._calculate_radial_velocity(self.observer, target)
        fin_obs_vel = self._calculate_radial_velocity(observer, target)

        # Project the velocity shift vector onto the the line-on-sight vector
        #  between the target and the new observation frame.
        init_proj_vel = np.dot(init_obs_vel, fin_obs_vel) / np.linalg.norm(
            fin_obs_vel)

        # Calculate the velocity shift between the two vectors
        delta_vel = fin_obs_vel - init_proj_vel

        # Apply the velocity shift to the stored spectral data
        new_data = self * (1 + delta_vel / c.cgs)
        new_coord = self._copy(value=new_data,
                               observer=observer,
                               target=target)

        return new_coord

    def in_observer_velocity_frame(self, frame):
        """
        Alters the velocity frame of the observer, but not the position.

        Parameters
        ----------
        frame : `BaseCoordinateFrame` or `SkyCoord`
            The observation frame containing the new velocity for the observer.

        Returns
        -------
        `SpectralCoord`
            The new coordinate object representing the spectral data
            transformed based on the observer's new velocity frame.
        """
        if hasattr(frame, 'frame'):
            frame = frame.frame

        if not frame.data.differentials:
            raise ValueError("Frame has no velocities, cannot transform "
                             "velocity frame.")

        data_with_rv = self.observer.data.with_differentials(
            frame.data.differentials)

        observer = self.observer.realize_frame(data_with_rv)

        return self._change_observer_to(observer)

    def with_los_shift(self, target=None, observer=None):
        """
        Apply a velocity shift to either the target or the observer. The shift
        can be provided as a redshift (float value) or radial velocity
        (quantity with physical type of 'speed').

        Parameters
        ----------
        target : float or `Quantity`
            Shift value to apply to current target.
        observer : float or `Quantity`
            Shift value to apply to current observer.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate with the target/observer velocity changed
            to incorporate the shift.
        """
        if self.target is None or self.observer is None:
            raise ValueError("Both an observer and target must be defined "
                             "before applying a velocity shift.")

        for arg in [x for x in [target, observer] if x is not None]:
            if isinstance(arg, u.Quantity) and arg.unit.physical_type != 'speed':
                raise u.UnitsError("Argument must has unit physical type "
                                   "'speed'.")

        target_icrs = self.target.transform_to(ICRS)
        observer_icrs = self.observer.transform_to(ICRS)

        target_shift = 0 * u.km / u.s if target is None else target
        observer_shift = 0 * u.km / u.s if observer is None else observer

        d_pos = (target_icrs.data.without_differentials() -
                 observer_icrs.data.without_differentials()).to_cartesian()

        pos_hat = d_pos / (d_pos.norm() or 1)

        target_velocity = target_icrs.velocity + target_shift * pos_hat
        observer_velocity = observer_icrs.velocity + observer_shift * pos_hat

        new_target = self.target.realize_frame(
            self.target.data.with_differentials(
                CartesianDifferential(target_velocity)))

        new_observer = self.observer.realize_frame(
            self.observer.data.with_differentials(
                CartesianDifferential(observer_velocity)))

        return self._copy(observer=new_observer,
                          target=new_target)

    def to(self, *args, doppler_rest=None, doppler_convention=None, **kwargs):
        """
        Overloaded parent ``to`` method to provide parameters for defining
        rest value and pre-defined conventions for unit transformations.

        Parameters
        ----------
        doppler_rest : `Quantity`, optional
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
            self.doppler_rest = doppler_rest

        if doppler_convention is not None:
            if doppler_convention not in DOPPLER_CONVENTIONS:
                raise ValueError(
                    "Unrecognized doppler convention: {}.".format(
                        doppler_convention))

            self.doppler_convention = doppler_convention

        equivs = u.spectral()

        if self.doppler_rest is not None and self.doppler_convention is not None:
            vel_equiv = DOPPLER_CONVENTIONS[self.doppler_convention](self.rest)
            equivs += vel_equiv

        # Compose the equivalencies for spectral conversions including the
        # appropriate velocity handling.
        kwargs.setdefault('equivalencies',
                          kwargs.pop('equivalencies', []) + equivs)

        return super().to(*args, **kwargs)

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
            f'radial_velocity={radial_velocity}, ' \
            f'redshift={redshift}, ' \
            f'\trest_value={self.doppler_rest}, ' \
            f'velocity_convention={self.doppler_convention}, ' \
            f'observer={obs_frame}, target={tar_frame}>'
