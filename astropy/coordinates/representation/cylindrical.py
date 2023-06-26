# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Cylindrical representations and differentials."""

import operator

import numpy as np

import astropy.units as u
from astropy.coordinates.angles import Angle

from .base import BaseDifferential, BaseRepresentation
from .cartesian import CartesianRepresentation
from .spherical import _spherical_op_funcs


class CylindricalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cylindrical coordinates.

    Parameters
    ----------
    rho : `~astropy.units.Quantity`
        The distance from the z axis to the point(s).

    phi : `~astropy.units.Quantity` or str
        The azimuth of the point(s), in angular units, which will be wrapped
        to an angle between 0 and 360 degrees. This can also be instances of
        `~astropy.coordinates.Angle`,

    z : `~astropy.units.Quantity`
        The z coordinate(s) of the point(s)

    differentials : dict, `~astropy.coordinates.CylindricalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `~astropy.coordinates.CylindricalDifferential` instance, or a dictionary of of differential
        instances with keys set to a string representation of the SI unit with
        which the differential (derivative) is taken. For example, for a
        velocity differential on a positional representation, the key would be
        ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"rho": u.Quantity, "phi": Angle, "z": u.Quantity}

    def __init__(self, rho, phi=None, z=None, differentials=None, copy=True):
        super().__init__(rho, phi, z, copy=copy, differentials=differentials)

        if not self._rho.unit.is_equivalent(self._z.unit):
            raise u.UnitsError("rho and z should have matching physical types")

    @property
    def rho(self):
        """
        The distance of the point(s) from the z-axis.
        """
        return self._rho

    @property
    def phi(self):
        """
        The azimuth of the point(s).
        """
        return self._phi

    @property
    def z(self):
        """
        The height of the point(s).
        """
        return self._z

    def unit_vectors(self):
        sinphi, cosphi = np.sin(self.phi), np.cos(self.phi)
        l = np.broadcast_to(1.0, self.shape)
        return {
            "rho": CartesianRepresentation(cosphi, sinphi, 0, copy=False),
            "phi": CartesianRepresentation(-sinphi, cosphi, 0, copy=False),
            "z": CartesianRepresentation(0, 0, l, unit=u.one, copy=False),
        }

    def scale_factors(self):
        rho = self.rho / u.radian
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"rho": l, "phi": rho, "z": l}

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to cylindrical polar
        coordinates.
        """
        rho = np.hypot(cart.x, cart.y)
        phi = np.arctan2(cart.y, cart.x)
        z = cart.z

        return cls(rho=rho, phi=phi, z=z, copy=False)

    def to_cartesian(self):
        """
        Converts cylindrical polar coordinates to 3D rectangular cartesian
        coordinates.
        """
        x = self.rho * np.cos(self.phi)
        y = self.rho * np.sin(self.phi)
        z = self.z

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    def _scale_operation(self, op, *args):
        if any(
            differential.base_representation is not self.__class__
            for differential in self.differentials.values()
        ):
            return super()._scale_operation(op, *args)

        phi_op, _, rho_op = _spherical_op_funcs(op, *args)
        z_op = lambda x: op(x, *args)

        result = self.__class__(
            rho_op(self.rho), phi_op(self.phi), z_op(self.z), copy=False
        )
        for key, differential in self.differentials.items():
            new_comps = (
                op(getattr(differential, comp))
                for op, comp in zip(
                    (rho_op, operator.pos, z_op), differential.components
                )
            )
            result.differentials[key] = differential.__class__(*new_comps, copy=False)
        return result


class CylindricalDifferential(BaseDifferential):
    """Differential(s) of points in cylindrical coordinates.

    Parameters
    ----------
    d_rho : `~astropy.units.Quantity` ['speed']
        The differential cylindrical radius.
    d_phi : `~astropy.units.Quantity` ['angular speed']
        The differential azimuth.
    d_z : `~astropy.units.Quantity` ['speed']
        The differential height.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = CylindricalRepresentation

    def __init__(self, d_rho, d_phi=None, d_z=None, copy=False):
        super().__init__(d_rho, d_phi, d_z, copy=copy)
        if not self._d_rho.unit.is_equivalent(self._d_z.unit):
            raise u.UnitsError("d_rho and d_z should have equivalent units.")
