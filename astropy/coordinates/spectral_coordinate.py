import warnings
from textwrap import indent

import astropy.units as u
import numpy as np
from astropy.constants import c
from astropy.coordinates import (ICRS,
                                 CartesianDifferential,
                                 CartesianRepresentation, SkyCoord)
from astropy.coordinates.spectral_quantity import SpectralQuantity
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           frame_transform_graph)
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ['SpectralCoord']


class NoVelocityWarning(AstropyUserWarning):
    pass


class NoDistanceWarning(AstropyUserWarning):
    pass


KMS = u.km / u.s
ZERO_VELOCITIES = CartesianDifferential([0, 0, 0] * KMS)

# Default distance to use for target when none is provided
DEFAULT_DISTANCE = 1e6 * u.kpc

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralCoord.*']


def _apply_relativistic_doppler_shift(scoord, velocity):
    """
    Given a `SpectralQuantity` and a velocity, return a new `SpectralQuantity`
    that is Doppler shifted by this amount.

    Note that the Doppler shift applied is the full relativistic one, so
    `SpectralQuantity` currently expressed in velocity and not using the
    relativistic convention will temporarily be converted to use the
    relativistic convention while the shift is applied.

    Positive velocities are assumed to redshift the spectral quantity,
    while negative velocities blueshift the spectral quantity.
    """

    # NOTE: we deliberately don't keep sub-classes of SpectralQuantity intact
    # since we can't guarantee that their metadata would be correct/consistent.
    squantity = scoord.view(SpectralQuantity)

    beta = velocity / c
    doppler_factor = np.sqrt((1 + beta) / (1 - beta))

    if squantity.unit.is_equivalent(u.m):  # wavelength
        return squantity * doppler_factor
    elif (squantity.unit.is_equivalent(u.Hz) or
          squantity.unit.is_equivalent(u.eV) or
          squantity.unit.is_equivalent(1 / u.m)):
        return squantity / doppler_factor
    elif squantity.unit.is_equivalent(KMS):  # velocity
        return (squantity.to(u.Hz) / doppler_factor).to(squantity.unit)
    else:  # pragma: no cover
        raise RuntimeError(f"Unexpected units in velocity shift: {squantity.unit}. "
                           "This should not happen, so please report this in the "
                           "astropy issue tracker!")


def update_differentials_to_match(original, velocity_reference, preserve_observer_frame=False):
    """
    Given an original coordinate object, update the differentials so that
    the final coordinate is at the same location as the original coordinate
    but co-moving with the velocity reference object.

    If preserve_original_frame is set to True, the resulting object will be in
    the frame of the original coordinate, otherwise it will be in the frame of
    the velocity reference.
    """

    if not velocity_reference.data.differentials:
        raise ValueError("Reference frame has no velocities")

    # If the reference has an obstime already defined, we should ignore
    # it and stick with the original observer obstime.
    if 'obstime' in velocity_reference.frame_attributes and hasattr(original, 'obstime'):
        velocity_reference = velocity_reference.replicate(obstime=original.obstime)

    # We transform both coordinates to ICRS for simplicity and because we know
    # it's a simple frame that is not time-dependent (it could be that both
    # the original and velocity_reference frame are time-dependent)

    original_icrs = original.transform_to(ICRS())
    velocity_reference_icrs = velocity_reference.transform_to(ICRS())

    differentials = velocity_reference_icrs.data.represent_as(CartesianRepresentation,
                                                              CartesianDifferential).differentials

    data_with_differentials = (original_icrs.data.represent_as(CartesianRepresentation)
                               .with_differentials(differentials))

    final_icrs = original_icrs.realize_frame(data_with_differentials)

    if preserve_observer_frame:
        final = final_icrs.transform_to(original)
    else:
        final = final_icrs.transform_to(velocity_reference)

    return final.replicate(representation_type=CartesianRepresentation,
                           differential_type=CartesianDifferential)


def attach_zero_velocities(coord):
    """
    Set the differentials to be stationary on a coordinate object.
    """
    new_data = coord.cartesian.with_differentials(ZERO_VELOCITIES)
    return coord.realize_frame(new_data)


def _get_velocities(coord):
    if 's' in coord.data.differentials:
        return coord.velocity
    else:
        return ZERO_VELOCITIES


class SpectralCoord(SpectralQuantity):
    """
    A spectral coordinate with its corresponding unit.

    .. note:: The |SpectralCoord| class is new in Astropy v4.1 and should be
              considered experimental at this time. Note that we do not fully
              support cases where the observer and target are moving
              relativistically relative to each other, so care should be taken
              in those cases. It is possible that there will be API changes in
              future versions of Astropy based on user feedback. If you have
              specific ideas for how it might be improved, please  let us know
              on the `astropy-dev mailing list`_ or at
              http://feedback.astropy.org.

    Parameters
    ----------
    value : ndarray or `~astropy.units.Quantity` or `SpectralCoord`
        Spectral values, which should be either wavelength, frequency,
        energy, wavenumber, or velocity values.
    unit : unit-like
        Unit for the given spectral values.
    observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
        The coordinate (position and velocity) of observer. If no velocities
        are present on this object, the observer is assumed to be stationary
        relative to the frame origin.
    target : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
        The coordinate (position and velocity) of target. If no velocities
        are present on this object, the target is assumed to be stationary
        relative to the frame origin.
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional
        The radial velocity of the target with respect to the observer. This
        can only be specified if ``redshift`` is not specified.
    redshift : float, optional
        The relativistic redshift of the target with respect to the observer.
        This can only be specified if ``radial_velocity`` cannot be specified.
    doppler_rest : `~astropy.units.Quantity`, optional
        The rest value to use when expressing the spectral value as a velocity.
    doppler_convention : str, optional
        The Doppler convention to use when expressing the spectral value as a velocity.
    """

    @u.quantity_input(radial_velocity=u.km/u.s)
    def __new__(cls, value, unit=None,
                observer=None, target=None,
                radial_velocity=None, redshift=None,
                **kwargs):

        obj = super().__new__(cls, value, unit=unit, **kwargs)

        # There are two main modes of operation in this class. Either the
        # observer and target are both defined, in which case the radial
        # velocity and redshift are automatically computed from these, or
        # only one of the observer and target are specified, along with a
        # manually specified radial velocity or redshift. So if a target and
        # observer are both specified, we can't also accept a radial velocity
        # or redshift.
        if target is not None and observer is not None:
            if radial_velocity is not None or redshift is not None:
                raise ValueError("Cannot specify radial velocity or redshift if both "
                                 "target and observer are specified")

        # We only deal with redshifts here and in the redshift property.
        # Otherwise internally we always deal with velocities.
        if redshift is not None:
            if radial_velocity is not None:
                raise ValueError("Cannot set both a radial velocity and redshift")
            redshift = u.Quantity(redshift)
            # For now, we can't specify redshift=u.one in quantity_input above
            # and have it work with plain floats, but if that is fixed, for
            # example as in https://github.com/astropy/astropy/pull/10232, we
            # can remove the check here and add redshift=u.one to the decorator
            if not redshift.unit.is_equivalent(u.one):
                raise u.UnitsError('redshift should be dimensionless')
            radial_velocity = redshift.to(u.km / u.s, u.doppler_redshift())

        # If we're initializing from an existing SpectralCoord, keep any
        # parameters that aren't being overridden
        if observer is None:
            observer = getattr(value, 'observer', None)
        if target is None:
            target = getattr(value, 'target', None)

        # As mentioned above, we should only specify the radial velocity
        # manually if either or both the observer and target are not
        # specified.
        if observer is None or target is None:
            if radial_velocity is None:
                radial_velocity = getattr(value, 'radial_velocity', None)

        obj._radial_velocity = radial_velocity
        obj._observer = cls._validate_coordinate(observer, label='observer')
        obj._target = cls._validate_coordinate(target, label='target')

        return obj

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        self._radial_velocity = getattr(obj, '_radial_velocity', None)
        self._observer = getattr(obj, '_observer', None)
        self._target = getattr(obj, '_target', None)

    @staticmethod
    def _validate_coordinate(coord, label=''):
        """
        Checks the type of the frame and whether a velocity differential and a
        distance has been defined on the frame object.

        If no distance is defined, the target is assumed to be "really far
        away", and the observer is assumed to be "in the solar system".

        Parameters
        ----------
        coord : `~astropy.coordinates.BaseCoordinateFrame`
            The new frame to be used for target or observer.
        label : str, optional
            The name of the object being validated (e.g. 'target' or 'observer'),
            which is then used in error messages.
        """

        if coord is None:
            return

        if not issubclass(coord.__class__, BaseCoordinateFrame):
            if isinstance(coord, SkyCoord):
                coord = coord.frame
            else:
                raise TypeError(f"{label} must be a SkyCoord or coordinate frame instance")

        # If the distance is not well-defined, ensure that it works properly
        # for generating differentials
        # TODO: change this to not set the distance and yield a warning once
        # there's a good way to address this in astropy.coordinates
        # https://github.com/astropy/astropy/issues/10247
        with np.errstate(all='ignore'):
            distance = getattr(coord, 'distance', None)
        if distance is not None and distance.unit.physical_type == 'dimensionless':
            coord = SkyCoord(coord, distance=DEFAULT_DISTANCE)
            warnings.warn(
                "Distance on coordinate object is dimensionless, an "
                f"arbitrary distance value of {DEFAULT_DISTANCE} will be set instead.",
                NoDistanceWarning)

        # If the observer frame does not contain information about the
        # velocity of the system, assume that the velocity is zero in the
        # system.
        if 's' not in coord.data.differentials:
            warnings.warn(
                f"No velocity defined on frame, assuming {ZERO_VELOCITIES}.",
                NoVelocityWarning)

            coord = attach_zero_velocities(coord)

        return coord

    def replicate(self, value=None, unit=None,
                  observer=None, target=None,
                  radial_velocity=None, redshift=None,
                  doppler_convention=None, doppler_rest=None,
                  copy=False):
        """
        Return a replica of the `SpectralCoord`, optionally changing the
        values or attributes.

        Note that no conversion is carried out by this method - this keeps
        all the values and attributes the same, except for the ones explicitly
        passed to this method which are changed.

        If ``copy`` is set to `True` then a full copy of the internal arrays
        will be made.  By default the replica will use a reference to the
        original arrays when possible to save memory.

        Parameters
        ----------
        value : ndarray or `~astropy.units.Quantity` or `SpectralCoord`, optional
            Spectral values, which should be either wavelength, frequency,
            energy, wavenumber, or velocity values.
        unit : unit-like
            Unit for the given spectral values.
        observer : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
            The coordinate (position and velocity) of observer.
        target : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`, optional
            The coordinate (position and velocity) of target.
        radial_velocity : `~astropy.units.Quantity` ['speed'], optional
            The radial velocity of the target with respect to the observer.
        redshift : float, optional
            The relativistic redshift of the target with respect to the observer.
        doppler_rest : `~astropy.units.Quantity`, optional
            The rest value to use when expressing the spectral value as a velocity.
        doppler_convention : str, optional
            The Doppler convention to use when expressing the spectral value as a velocity.
        copy : bool, optional
            If `True`, and ``value`` is not specified, the values are copied to
            the new `SkyCoord` - otherwise a reference to the same values is used.

        Returns
        -------
        sc : `SpectralCoord` object
            Replica of this object
        """

        if isinstance(value, u.Quantity):
            if unit is not None:
                raise ValueError("Cannot specify value as a Quantity and also specify unit")
            else:
                value, unit = value.value, value.unit

        value = value if value is not None else self.value
        unit = unit or self.unit
        observer = self._validate_coordinate(observer) or self.observer
        target = self._validate_coordinate(target) or self.target
        doppler_convention = doppler_convention or self.doppler_convention
        doppler_rest = doppler_rest or self.doppler_rest

        # If value is being taken from self and copy is Tru
        if copy:
            value = value.copy()

        # Only include radial_velocity if it is not auto-computed from the
        # observer and target.
        if (self.observer is None or self.target is None) and radial_velocity is None and redshift is None:
            radial_velocity = self.radial_velocity

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', NoVelocityWarning)
            return self.__class__(value=value, unit=unit,
                                  observer=observer, target=target,
                                  radial_velocity=radial_velocity, redshift=redshift,
                                  doppler_convention=doppler_convention, doppler_rest=doppler_rest, copy=False)

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
        The coordinates of the observer.

        If set, and a target is set as well, this will override any explicit
        radial velocity passed in.

        Returns
        -------
        `~astropy.coordinates.BaseCoordinateFrame`
            The astropy coordinate frame representing the observation.
        """
        return self._observer

    @observer.setter
    def observer(self, value):

        if self.observer is not None:
            raise ValueError("observer has already been set")

        self._observer = self._validate_coordinate(value, label='observer')

        # Switch to auto-computing radial velocity
        if self._target is not None:
            self._radial_velocity = None

    @property
    def target(self):
        """
        The coordinates of the target being observed.

        If set, and an observer is set as well, this will override any explicit
        radial velocity passed in.

        Returns
        -------
        `~astropy.coordinates.BaseCoordinateFrame`
            The astropy coordinate frame representing the target.
        """
        return self._target

    @target.setter
    def target(self, value):

        if self.target is not None:
            raise ValueError("target has already been set")

        self._target = self._validate_coordinate(value, label='target')

        # Switch to auto-computing radial velocity
        if self._observer is not None:
            self._radial_velocity = None

    @property
    def radial_velocity(self):
        """
        Radial velocity of target relative to the observer.

        Returns
        -------
        `~astropy.units.Quantity` ['speed']
            Radial velocity of target.

        Notes
        -----
        This is different from the ``.radial_velocity`` property of a
        coordinate frame in that this calculates the radial velocity with
        respect to the *observer*, not the origin of the frame.
        """
        if self._observer is None or self._target is None:
            if self._radial_velocity is None:
                return 0 * KMS
            else:
                return self._radial_velocity
        else:
            return self._calculate_radial_velocity(self._observer, self._target,
                                                   as_scalar=True)

    @property
    def redshift(self):
        """
        Redshift of target relative to observer. Calculated from the radial
        velocity.

        Returns
        -------
        `astropy.units.Quantity`
            Redshift of target.
        """
        return self.radial_velocity.to(u.dimensionless_unscaled, u.doppler_redshift())

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
        as_scalar : bool
            If `True`, the magnitude of the velocity vector will be returned,
            otherwise the full vector will be returned.

        Returns
        -------
        `~astropy.units.Quantity` ['speed']
            The radial velocity of the target with respect to the observer.
        """

        # Convert observer and target to ICRS to avoid finite differencing
        # calculations that lack numerical precision.
        observer_icrs = observer.transform_to(ICRS())
        target_icrs = target.transform_to(ICRS())

        pos_hat = SpectralCoord._normalized_position_vector(observer_icrs, target_icrs)

        d_vel = target_icrs.velocity - observer_icrs.velocity

        vel_mag = pos_hat.dot(d_vel)

        if as_scalar:
            return vel_mag
        else:
            return vel_mag * pos_hat

    @staticmethod
    def _normalized_position_vector(observer, target):
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
        d_pos = (target.cartesian.without_differentials() -
                 observer.cartesian.without_differentials())

        dp_norm = d_pos.norm()

        # Reset any that are 0 to 1 to avoid nans from 0/0
        dp_norm[dp_norm == 0] = 1 * dp_norm.unit

        pos_hat = d_pos / dp_norm

        return pos_hat

    @u.quantity_input(velocity=u.km/u.s)
    def with_observer_stationary_relative_to(self, frame, velocity=None, preserve_observer_frame=False):
        """
        A new  `SpectralCoord` with the velocity of the observer altered,
        but not the position.

        If a coordinate frame is specified, the observer velocities will be
        modified to be stationary in the specified frame. If a coordinate
        instance is specified, optionally with non-zero velocities, the
        observer velocities will be updated so that the observer is co-moving
        with the specified coordinates.

        Parameters
        ----------
        frame : str, `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The observation frame in which the observer will be stationary. This
            can be the name of a frame (e.g. 'icrs'), a frame class, frame instance
            with no data, or instance with data. This can optionally include
            velocities.
        velocity : `~astropy.units.Quantity` or `~astropy.coordinates.CartesianDifferential`, optional
            If ``frame`` does not contain velocities, these can be specified as
            a 3-element `~astropy.units.Quantity`. In the case where this is
            also not specified, the velocities default to zero.
        preserve_observer_frame : bool
            If `True`, the final observer frame class will be the same as the
            original one, and if `False` it will be the frame of the velocity
            reference class.

        Returns
        -------
        new_coord : `SpectralCoord`
            The new coordinate object representing the spectral data
            transformed based on the observer's new velocity frame.
        """

        if self.observer is None or self.target is None:
            raise ValueError("This method can only be used if both observer "
                             "and target are defined on the SpectralCoord.")

        # Start off by extracting frame if a SkyCoord was passed in
        if isinstance(frame, SkyCoord):
            frame = frame.frame

        if isinstance(frame, BaseCoordinateFrame):

            if not frame.has_data:
                frame = frame.realize_frame(CartesianRepresentation(0 * u.km, 0 * u.km, 0 * u.km))

            if frame.data.differentials:
                if velocity is not None:
                    raise ValueError('frame already has differentials, cannot also specify velocity')
                # otherwise frame is ready to go
            else:
                if velocity is None:
                    differentials = ZERO_VELOCITIES
                else:
                    differentials = CartesianDifferential(velocity)
                frame = frame.realize_frame(frame.data.with_differentials(differentials))

        if isinstance(frame, (type, str)):
            if isinstance(frame, type):
                frame_cls = frame
            elif isinstance(frame, str):
                frame_cls = frame_transform_graph.lookup_name(frame)
            if velocity is None:
                velocity = 0 * u.m / u.s, 0 * u.m / u.s, 0 * u.m / u.s
            elif velocity.shape != (3,):
                raise ValueError('velocity should be a Quantity vector with 3 elements')
            frame = frame_cls(0 * u.m, 0 * u.m, 0 * u.m,
                              *velocity,
                              representation_type='cartesian',
                              differential_type='cartesian')

        observer = update_differentials_to_match(self.observer, frame,
                                                 preserve_observer_frame=preserve_observer_frame)

        # Calculate the initial and final los velocity
        init_obs_vel = self._calculate_radial_velocity(self.observer, self.target, as_scalar=True)
        fin_obs_vel = self._calculate_radial_velocity(observer, self.target, as_scalar=True)

        # Apply transformation to data
        new_data = _apply_relativistic_doppler_shift(self, fin_obs_vel - init_obs_vel)

        new_coord = self.replicate(value=new_data, observer=observer)

        return new_coord

    def with_radial_velocity_shift(self, target_shift=None, observer_shift=None):
        """
        Apply a velocity shift to this spectral coordinate.

        The shift can be provided as a redshift (float value) or radial
        velocity (`~astropy.units.Quantity` with physical type of 'speed').

        Parameters
        ----------
        target_shift : float or `~astropy.units.Quantity` ['speed']
            Shift value to apply to current target.
        observer_shift : float or `~astropy.units.Quantity` ['speed']
            Shift value to apply to current observer.

        Returns
        -------
        `SpectralCoord`
            New spectral coordinate with the target/observer velocity changed
            to incorporate the shift. This is always a new object even if
            ``target_shift`` and ``observer_shift`` are both `None`.
        """

        if observer_shift is not None and (self.target is None or
                                           self.observer is None):
            raise ValueError("Both an observer and target must be defined "
                             "before applying a velocity shift.")

        for arg in [x for x in [target_shift, observer_shift] if x is not None]:
            if isinstance(arg, u.Quantity) and not arg.unit.is_equivalent((u.one, KMS)):
                raise u.UnitsError("Argument must have unit physical type "
                                   "'speed' for radial velocty or "
                                   "'dimensionless' for redshift.")

        # The target or observer value is defined but is not a quantity object,
        #  assume it's a redshift float value and convert to velocity

        if target_shift is None:
            if self._observer is None or self._target is None:
                return self.replicate()
            target_shift = 0 * KMS
        else:
            target_shift = u.Quantity(target_shift)
            if target_shift.unit.physical_type == 'dimensionless':
                target_shift = target_shift.to(u.km / u.s, u.doppler_redshift())
            if self._observer is None or self._target is None:
                return self.replicate(value=_apply_relativistic_doppler_shift(self, target_shift),
                                      radial_velocity=self.radial_velocity + target_shift)

        if observer_shift is None:
            observer_shift = 0 * KMS
        else:
            observer_shift = u.Quantity(observer_shift)
            if observer_shift.unit.physical_type == 'dimensionless':
                observer_shift = observer_shift.to(u.km / u.s, u.doppler_redshift())

        target_icrs = self._target.transform_to(ICRS())
        observer_icrs = self._observer.transform_to(ICRS())

        pos_hat = SpectralCoord._normalized_position_vector(observer_icrs, target_icrs)

        target_velocity = _get_velocities(target_icrs) + target_shift * pos_hat
        observer_velocity = _get_velocities(observer_icrs) + observer_shift * pos_hat

        target_velocity = CartesianDifferential(target_velocity.xyz)
        observer_velocity = CartesianDifferential(observer_velocity.xyz)

        new_target = (target_icrs
                      .realize_frame(target_icrs.cartesian.with_differentials(target_velocity))
                      .transform_to(self._target))

        new_observer = (observer_icrs
                        .realize_frame(observer_icrs.cartesian.with_differentials(observer_velocity))
                        .transform_to(self._observer))

        init_obs_vel = self._calculate_radial_velocity(observer_icrs, target_icrs, as_scalar=True)
        fin_obs_vel = self._calculate_radial_velocity(new_observer, new_target, as_scalar=True)

        new_data = _apply_relativistic_doppler_shift(self, fin_obs_vel - init_obs_vel)

        return self.replicate(value=new_data,
                              observer=new_observer,
                              target=new_target)

    def to_rest(self):
        """
        Transforms the spectral axis to the rest frame.
        """

        if self.observer is not None and self.target is not None:
            return self.with_observer_stationary_relative_to(self.target)

        result = _apply_relativistic_doppler_shift(self, -self.radial_velocity)

        return self.replicate(value=result, radial_velocity=0. * KMS, redshift=None)

    def __repr__(self):

        prefixstr = '<' + self.__class__.__name__ + ' '

        try:
            radial_velocity = self.radial_velocity
            redshift = self.redshift
        except ValueError:
            radial_velocity = redshift = 'Undefined'

        repr_items = [f'{prefixstr}']

        if self.observer is not None:
            observer_repr = indent(repr(self.observer), 14 * ' ').lstrip()
            repr_items.append(f'    observer: {observer_repr}')

        if self.target is not None:
            target_repr = indent(repr(self.target), 12 * ' ').lstrip()
            repr_items.append(f'    target: {target_repr}')

        if (self._observer is not None and self._target is not None) or self._radial_velocity is not None:
            if self.observer is not None and self.target is not None:
                repr_items.append('    observer to target (computed from above):')
            else:
                repr_items.append('    observer to target:')
            repr_items.append(f'      radial_velocity={radial_velocity}')
            repr_items.append(f'      redshift={redshift}')

        if self.doppler_rest is not None or self.doppler_convention is not None:
            repr_items.append(f'    doppler_rest={self.doppler_rest}')
            repr_items.append(f'    doppler_convention={self.doppler_convention}')

        arrstr = np.array2string(self.view(np.ndarray), separator=', ',
                                 prefix='  ')

        if len(repr_items) == 1:
            repr_items[0] += f'{arrstr}{self._unitstr:s}'
        else:
            repr_items[1] = '   (' + repr_items[1].lstrip()
            repr_items[-1] += ')'
            repr_items.append(f'  {arrstr}{self._unitstr:s}')

        return '\n'.join(repr_items) + '>'
