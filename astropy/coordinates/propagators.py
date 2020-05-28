"""
Propagators that define the ICRS<->BCRS transformations
"""
from abc import ABC, abstractmethod

import numpy as np

from astropy import log
from astropy.constants import c as speed_of_light
from astropy.time import Time
from astropy.utils.state import ScienceState
import astropy.units as u

from .builtin_frames import BCRS, ICRS
from .sky_coordinate import SkyCoord


__all__ = ['transform_propagator', 'BasePropagator',
           'LinearPropagator', 'SolarSystemLinearPropagator']


DEFAULT_REF_EPOCH = Time('J2000', scale='tt')


class transform_propagator(ScienceState):
    _value = None

    @classmethod
    def validate(cls, value):
        if value is not None and not isinstance(value, BasePropagator):
            raise ValueError("Must be a `Propagator` object or `None`.")
        return value


class BasePropagator(ABC):
    attributes = []

    def __repr__(self):
        stubs = []
        for name in self.attributes:
            stubs.append(f"{name}={getattr(self, name)}")
        return f"<{self.__class__.__name__} ({', '.join(stubs)})>"

    def __getattr__(self, name):
        if name not in self.attributes:
            raise AttributeError("Invalid attribute name")
        return getattr(self, f"_{name}", None)

    @abstractmethod
    def propagate(self, coord, to_time=None):
        pass

    @abstractmethod
    def depropagate(self, bcrs_coord):
        pass


class LinearPropagator(BasePropagator):
    """
    ICRS<->BCRS assuming linear motion
    """
    attributes = ['ref_epoch']

    def __init__(self, ref_epoch=DEFAULT_REF_EPOCH):
        self._ref_epoch = Time(ref_epoch)

    def propagate(self, icrs_coord, to_time=None):
        if not isinstance(icrs_coord, ICRS):
            raise ValueError("Coordinate must be in ICRS")
        if to_time is None:
            to_time = icrs_coord.obstime

        if 's' in icrs_coord.data.differentials:
            dt = (to_time - self.ref_epoch).to(u.yr)
            log.debug(f"Propagating forward by {dt} after {self.ref_epoch}")
            # TODO: code for SkyCoord.apply_space_motion() needs to moved out of that class
            result = SkyCoord(icrs_coord).apply_space_motion(dt=dt).frame
            return BCRS(result.data, obstime=to_time)
        else:
            log.debug("The coordinate does not have velocity information.")
            return BCRS(icrs_coord.data, obstime=to_time)

    def depropagate(self, bcrs_coord):
        if not isinstance(bcrs_coord, BCRS):
            raise ValueError("Coordinate must be in BCRS")
        from_time = bcrs_coord.obstime

        if 's' in bcrs_coord.data.differentials:
            dt = (from_time - self.ref_epoch).to(u.yr)
            log.debug(f"Propagating backward by {dt} to {self.ref_epoch}")
            # TODO: code for SkyCoord.apply_space_motion() needs to moved out of that class
            result = SkyCoord(ICRS(bcrs_coord.data)).apply_space_motion(dt=-dt).frame
            return ICRS(result.data)
        else:
            log.debug("The coordinate does not have velocity information.")
            return ICRS(bcrs_coord.data)


class SolarSystemLinearPropagator(BasePropagator):
    """
    ICRS<->BCRS assuming linear motion, with special handling of solar-system
    bodies

    Requires 3D coordinates because it uses the barycentric distance to
    distinguish between solar-system bodies and cosmic objects.
    """
    attributes = ['observer', 'distance_threshold', 'ref_epoch']

    @u.quantity_input
    def __init__(self, observer, distance_threshold: u.lyr = 1*u.lyr, ref_epoch=DEFAULT_REF_EPOCH):
        self._observer = observer.transform_to(ICRS)
        if hasattr(self._observer, "frame"):
            self._observer = self._observer.frame
        self._distance_threshold = distance_threshold
        self._ref_epoch = Time(ref_epoch)

        self._linear_propagator = LinearPropagator(ref_epoch=self.ref_epoch)

    def propagate(self, icrs_coord, to_time=None):
        if not isinstance(icrs_coord, ICRS):
            raise ValueError("Coordinate must be in ICRS")
        if to_time is None:
            to_time = icrs_coord.obstime

        posvel = icrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < self.distance_threshold):  # solar-system body
                pos = posvel.without_differentials()
                vel = posvel.differentials['s'].to_cartesian()

                # The light travel time (dt) for a body in linear motion satisfies the equation:
                #   norm(D - V * dt) = c * dt
                # where D is the observer-body vector, V is the body velocity vector,
                # and c is the speed of light.  This turns into a quadratic equation in dt, and
                # only the positive solution is valid:
                #   dt = (sqrt((D dot V)^2 + (D dot D) * c2_V2) - D dot V) / c2_V2
                # where c2_V2 = c^2 - V dot V
                D = pos - self.observer.cartesian
                c2_V2 = speed_of_light ** 2 - vel.dot(vel)
                DdotV = D.dot(vel)
                dt = (np.sqrt(DdotV ** 2 + D.dot(D) * c2_V2) - DdotV) / c2_V2
                log.debug(f"Retarding for {dt.to(u.s)} of light travel time")

                pos -= vel * dt
                posvel = pos.with_differentials(posvel.differentials)
                return BCRS(posvel, obstime=to_time)
            elif np.all(distance >= self.distance_threshold):  # cosmic object
                return self._linear_propagator.propagate(icrs_coord, to_time)
            else:  # a mix
                raise ConvertError("The propagation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")
        else:
            log.debug("The coordinate does not have velocity information.")
            return BCRS(posvel, obstime=to_time)

    def depropagate(self, bcrs_coord):
        if not isinstance(bcrs_coord, BCRS):
            raise ValueError("Coordinate must be in BCRS")
        from_time = bcrs_coord.obstime

        posvel = bcrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < self.distance_threshold):  # solar-system body
                pos = posvel.without_differentials()

                # The light travel time is straightforward to compute for an astrometric location
                dt = (pos - self.observer.cartesian).norm() / speed_of_light
                log.debug(f"Advancing for {dt.to(u.s)} of light travel time")

                pos += posvel.differentials['s'] * dt
                posvel = pos.with_differentials(posvel.differentials)
                return ICRS(posvel)
            elif np.all(distance >= self.distance_threshold):  # cosmic object
                return self._linear_propagator.depropagate(bcrs_coord)
            else:  # a mix
                raise ConvertError("The propagation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")
        else:
            log.debug("The coordinate does not have velocity information.")
            return ICRS(posvel)
