import warnings

import numpy as np

import astropy.units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord, ICRS, Distance, GCRS
from astropy.coordinates.baseframe import BaseCoordinateFrame, FrameMeta
from astropy.utils.compat import NUMPY_LT_1_14
from astropy.utils.exceptions import AstropyUserWarning

DOPPLER_CONVENTIONS = {
    'radio': u.doppler_radio,
    'optical': u.doppler_optical,
    'relativistic': u.doppler_relativistic
}

RV_RS_EQUIV = [(u.cm / u.s, u.Unit(''),
                lambda x: x/c.cgs.value - 1,
                lambda x: (x + 1) * c.cgs.value)]

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
    radial_velocity : `Quantity`, optional
        The radial velocity of the target with respect to the observer.
    redshift : float, optional
        The redshift of the target with respect to the observer.
    rest : `Quantity`, optional
        The rest value to use for velocity space transformations.
    velocity_convention : str
        The convention to use when converting the spectral data to/from
        velocity space.
    observer : `BaseCoordinateFrame` or `SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    target : `BaseCoordinateFrame` or `SkyCoord`, optional
        The coordinate (position and velocity) of observer.
    """
    def __new__(cls, value, unit=None, radial_velocity=None, redshift=None,
                rest=None, velocity_convention=None, observer=None,
                target=None, **kwargs):
        obj = super().__new__(cls, value, unit=unit, **kwargs)

        obj.radial_velocity = radial_velocity
        obj.redshift = redshift

        obj.rest = rest
        obj.velocity_convention = velocity_convention

        obj.observer = observer
        obj.target = target

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)

        # If we're a new object or viewing an ndarray, nothing has to be done.
        if obj is None or obj.__class__ is np.ndarray:
            return

        self.radial_velocity = getattr(obj, 'radial_velocity', None)
        self.redshift = getattr(obj, 'redshift', None)

        self.rest = getattr(obj, 'rest', None)
        self.velocity_convention = getattr(obj, 'velocity_convention', None)

        self.observer = getattr(obj, 'observer', None)
        self.target = getattr(obj, 'target', None)

    def __quantity_subclass__(self, unit):
        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.
        """
        return SpectralCoord, True

    @staticmethod
    def _validate_frame(frame, is_observer):
        """
        Checks the type of the frame and whether a velocity differential and a
        distance has been defined on the frame object.

        If not distance is defined, the target is assumed to be "really far
        away", and the observer is assumed to be "in the solar system".

        Parameters
        ----------
        frame : `BaseCoordinateFrame`
            The new frame to be used for target or observer.
        is_observer : bool
            Whether the frame represents the observer or target.
        """
        if frame is not None:
            frame_loc = 'observer' if is_observer else 'target'

            if not issubclass(frame.__class__, (BaseCoordinateFrame, FrameMeta)):
                if isinstance(frame, SkyCoord):
                    frame = frame.frame
                else:
                    raise ValueError("`{}` is not a subclass of "
                                     "`BaseCoordinateFrame` or "
                                     "`SkyCoord`.".format(frame))

            # If no distance value is defined on the frame, assume a default
            # distance of either zero (for an observer), or 1000 kpc for target
            if not isinstance(frame.distance, Distance):
                auto_dist = 0 * u.AU if is_observer else 1000 * u.kpc

                warnings.warn(
                    "No distance defined on the {} frame, assuming {}.".format(
                        frame_loc, auto_dist),
                    AstropyUserWarning)

                frame = SkyCoord(frame, distance=auto_dist).frame

            # If the observer frame does not contain information about the
            # velocity of the system, assume that the velocity is zero in the
            # system.
            if 's' not in frame.data.differentials:
                warnings.warn(
                    "No velocity defined on {} frame, assuming {}.".format(
                       frame_loc, u.Quantity([0, 0, 0], unit=u.km/u.s)),
                    AstropyUserWarning)

                frame = SkyCoord(frame, radial_velocity=0 * u.km / u.s).frame

        return frame

    def _copy(self, **kwargs):
        default_kwargs = {
            'value': self.value,
            'unit': self.unit,
            'radial_velocity': self.radial_velocity,
            'rest': self.rest,
            'velocity_convention': self.velocity_convention,
            'observer': self.observer,
            'target': self.target
        }

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
        value = self._validate_frame(value, True)
        self._observer = value

        # Reset radial velocity and redshift values
        self.radial_velocity = None

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
        value = self._validate_frame(value, False)
        self._target = value

        # Reset radial velocity and redshift values
        self.radial_velocity = None

    @property
    def rest(self):
        """
        The rest value of the spectrum used for transformations to/from
        velocity space.

        Returns
        -------
        `Quantity`
            Rest value as an astropy `Quantity` object.
        """
        return self._rest

    @rest.setter
    @u.quantity_input(value=['length', 'frequency', 'energy', 'speed', None])
    def rest(self, value):
        """
        New rest value needed for velocity-space conversions.

        Parameters
        ----------
        value : `Quantity`
            Rest value.
        """
        self._rest = value

    @property
    def velocity_convention(self):
        """
        The defined convention for conversions to/from velocity space.

        Returns
        -------
        str
            One of 'optical', 'radio', or 'relativistic' representing the
            equivalency used in the unit conversions.
        """
        return self._velocity_convention

    @velocity_convention.setter
    def velocity_convention(self, value):
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

        self._velocity_convention = value

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
        if self.velocity_convention is not None and self.rest is not None:
            self._radial_velocity = DOPPLER_CONVENTIONS.get(
                self.velocity_convention)(self.rest)(self)

        return self._radial_velocity

    @radial_velocity.setter
    @u.quantity_input(value=['speed', None])
    def radial_velocity(self, value):
        """
        New radial velocity of target.

        Parameters
        ----------
        value : `u.Quantity`
            New radial velocity of target relative to observer.
        """
        self._radial_velocity = value

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
        if self.radial_velocity is not None:
            return self.radial_velocity.to('', equivalencies=RV_RS_EQUIV).value

    @redshift.setter
    def redshift(self, value):
        """
        New redshift of target, updates the stored radial velocity value.

        Parameters
        ----------
        value : float
            New redshift of target.
        """
        if value is not None:
            self.radial_velocity = u.Quantity(value).to(
                'km/s', equivalencies=RV_RS_EQUIV)

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
        observer_icrs = observer.transform_to(ICRS)
        target_icrs = target.transform_to(observer_icrs)

        d_pos = (target_icrs.data.without_differentials() -
                 observer_icrs.data.without_differentials()).to_cartesian()
        pos_hat = d_pos / d_pos.norm()
        d_vel = target_icrs.velocity - observer_icrs.velocity

        return d_vel.d_xyz * pos_hat.xyz

    def transform_frame(self, frame, target=None):
        """
        Transform the frame of the observation.

        Parameters
        ----------
        frame : `BaseCoordinateFrame` or `SkyCoord`
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
        if target is not None:
            self.target = target

        # Check velocities and distance values on the new frame. This is
        # handled in the frame validation.
        # frame = self._validate_frame(frame, True)
        new_obs_frame = self.observer.transform_to(frame)

        # In the case the user has defined a custom radial velocity, use that
        # instead of re-calculating the radial velocity based on frames
        init_obs_vel = self._calculate_radial_velocity(self.observer,
                                                       self.target)
        fin_obs_vel = self._calculate_radial_velocity(new_obs_frame,
                                                      self.target)

        # Calculate the velocity shift between the two vectors
        obs_vel_shift = init_obs_vel - fin_obs_vel

        # Project the velocity shift vector onto the the line-on-sight vector
        # between the target and the new observation frame.
        delta_vel = np.dot(obs_vel_shift, fin_obs_vel) / np.linalg.norm(
            fin_obs_vel) * obs_vel_shift.unit

        # Store the new frame as the current observer frame
        self.observer = new_obs_frame

        # Apply the velocity shift to the stored spectral data
        new_data = (self.to('Hz') * (1 + delta_vel / c.cgs)).to(self.unit)
        new_coord = self._copy(value=new_data)

        return new_coord

    def to(self, *args, rest=None, velocity_convention='relativistic',
           **kwargs):
        """
        Overloaded parent ``to`` method to provide parameters for defining
        rest value and pre-defined conventions for unit transformations.

        Parameters
        ----------
        rest : `Quantity`, optional
            The rest value used in the velocity space conversions. Providing
            the value here will set the value stored on the `SpectralCoord`
            instance.
        velocity_convention : {'relativistic', 'optical', 'radio'}, optional
            The velocity convention to use during conversions. Providing the
            value here will set the value stored on the `SpectralCoord`
            instance.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate object with data converted to the new unit.
        """
        if rest is not None:
            self.rest = rest

        if velocity_convention is not None:
            self.velocity_convention = velocity_convention

        equivs = u.spectral()

        if self.rest is not None and self.velocity_convention is not None:
            vel_equiv = DOPPLER_CONVENTIONS[self.velocity_convention](self.rest)
            equivs += vel_equiv

        # Compose the equivalencies for spectral conversions including the
        # appropriate velocity handling.
        kwargs.setdefault('equivalencies',
                          kwargs.pop('equivalencies', []) + equivs)

        return super().to(*args, **kwargs)

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        sep = ',' if NUMPY_LT_1_14 else ', '
        arrstr = np.array2string(self.view(np.ndarray), separator=sep,
                                 prefix=prefixstr)
        obs_frame = self.observer.__class__.__name__ \
            if self.observer is not None else 'None'
        tar_frame = self.target.__class__.__name__ \
            if self.target is not None else 'None'
        return f'{prefixstr}{arrstr}{self._unitstr:s}, \n' \
            f'radial_velocity={self.radial_velocity}, ' \
            f'redshift={self.redshift}, ' \
            f'\trest_value={self.rest}, ' \
            f'velocity_convention={self.velocity_convention}, ' \
            f'observer={obs_frame}, target={tar_frame}>'

