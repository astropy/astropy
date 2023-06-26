# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Cartesian representations and differentials."""

import numpy as np
from erfa import ufunc as erfa_ufunc

import astropy.units as u

from .base import BaseDifferential, BaseRepresentation


class CartesianRepresentation(BaseRepresentation):
    """
    Representation of points in 3D cartesian coordinates.

    Parameters
    ----------
    x, y, z : `~astropy.units.Quantity` or array
        The x, y, and z coordinates of the point(s). If ``x``, ``y``, and ``z``
        have different shapes, they should be broadcastable. If not quantity,
        ``unit`` should be set.  If only ``x`` is given, it is assumed that it
        contains an array with the 3 coordinates stored along ``xyz_axis``.
    unit : unit-like
        If given, the coordinates will be converted to this unit (or taken to
        be in this unit if not given.
    xyz_axis : int, optional
        The axis along which the coordinates are stored when a single array is
        provided rather than distinct ``x``, ``y``, and ``z`` (default: 0).

    differentials : dict, `~astropy.coordinates.CartesianDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `~astropy.coordinates.CartesianDifferential` instance, or a dictionary of
        `~astropy.coordinates.CartesianDifferential` s with keys set to a string representation of
        the SI unit with which the differential (derivative) is taken. For
        example, for a velocity differential on a positional representation, the
        key would be ``'s'`` for seconds, indicating that the derivative is a
        time derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    attr_classes = {"x": u.Quantity, "y": u.Quantity, "z": u.Quantity}

    _xyz = None

    def __init__(
        self, x, y=None, z=None, unit=None, xyz_axis=None, differentials=None, copy=True
    ):
        if y is None and z is None:
            if isinstance(x, np.ndarray) and x.dtype.kind not in "OV":
                # Short-cut for 3-D array input.
                x = u.Quantity(x, unit, copy=copy, subok=True)
                # Keep a link to the array with all three coordinates
                # so that we can return it quickly if needed in get_xyz.
                self._xyz = x
                if xyz_axis:
                    x = np.moveaxis(x, xyz_axis, 0)
                    self._xyz_axis = xyz_axis
                else:
                    self._xyz_axis = 0

                self._x, self._y, self._z = x
                self._differentials = self._validate_differentials(differentials)
                return

            elif (
                isinstance(x, CartesianRepresentation)
                and unit is None
                and xyz_axis is None
            ):
                if differentials is None:
                    differentials = x._differentials

                return super().__init__(x, differentials=differentials, copy=copy)

            else:
                x, y, z = x

        if xyz_axis is not None:
            raise ValueError(
                "xyz_axis should only be set if x, y, and z are in a single array"
                " passed in through x, i.e., y and z should not be not given."
            )

        if y is None or z is None:
            raise ValueError(
                f"x, y, and z are required to instantiate {self.__class__.__name__}"
            )

        if unit is not None:
            x = u.Quantity(x, unit, copy=copy, subok=True)
            y = u.Quantity(y, unit, copy=copy, subok=True)
            z = u.Quantity(z, unit, copy=copy, subok=True)
            copy = False

        super().__init__(x, y, z, copy=copy, differentials=differentials)
        if not (
            self._x.unit.is_equivalent(self._y.unit)
            and self._x.unit.is_equivalent(self._z.unit)
        ):
            raise u.UnitsError("x, y, and z should have matching physical types")

    def unit_vectors(self):
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        o = np.broadcast_to(0.0 * u.one, self.shape, subok=True)
        return {
            "x": CartesianRepresentation(l, o, o, copy=False),
            "y": CartesianRepresentation(o, l, o, copy=False),
            "z": CartesianRepresentation(o, o, l, copy=False),
        }

    def scale_factors(self):
        l = np.broadcast_to(1.0 * u.one, self.shape, subok=True)
        return {"x": l, "y": l, "z": l}

    def get_xyz(self, xyz_axis=0):
        """Return a vector array of the x, y, and z coordinates.

        Parameters
        ----------
        xyz_axis : int, optional
            The axis in the final array along which the x, y, z components
            should be stored (default: 0).

        Returns
        -------
        xyz : `~astropy.units.Quantity`
            With dimension 3 along ``xyz_axis``.  Note that, if possible,
            this will be a view.
        """
        if self._xyz is not None:
            if self._xyz_axis == xyz_axis:
                return self._xyz
            else:
                return np.moveaxis(self._xyz, self._xyz_axis, xyz_axis)

        # Create combined array.  TO DO: keep it in _xyz for repeated use?
        # But then in-place changes have to cancel it. Likely best to
        # also update components.
        return np.stack([self._x, self._y, self._z], axis=xyz_axis)

    xyz = property(get_xyz)

    @classmethod
    def from_cartesian(cls, other):
        return other

    def to_cartesian(self):
        return self

    def transform(self, matrix):
        """
        Transform the cartesian coordinates using a 3x3 matrix.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : ndarray
            A 3x3 transformation matrix, such as a rotation matrix.

        Examples
        --------
        We can start off by creating a cartesian representation object:

            >>> from astropy import units as u
            >>> from astropy.coordinates import CartesianRepresentation
            >>> rep = CartesianRepresentation([1, 2] * u.pc,
            ...                               [2, 3] * u.pc,
            ...                               [3, 4] * u.pc)

        We now create a rotation matrix around the z axis:

            >>> from astropy.coordinates.matrix_utilities import rotation_matrix
            >>> rotation = rotation_matrix(30 * u.deg, axis='z')

        Finally, we can apply this transformation:

            >>> rep_new = rep.transform(rotation)
            >>> rep_new.xyz  # doctest: +FLOAT_CMP
            <Quantity [[ 1.8660254 , 3.23205081],
                       [ 1.23205081, 1.59807621],
                       [ 3.        , 4.        ]] pc>
        """
        # erfa rxp: Multiply a p-vector by an r-matrix.
        p = erfa_ufunc.rxp(matrix, self.get_xyz(xyz_axis=-1))
        # transformed representation
        rep = self.__class__(p, xyz_axis=-1, copy=False)
        # Handle differentials attached to this representation
        new_diffs = {
            k: d.transform(matrix, self, rep) for k, d in self.differentials.items()
        }
        return rep.with_differentials(new_diffs)

    def _combine_operation(self, op, other, reverse=False):
        self._raise_if_has_differentials(op.__name__)

        try:
            other_c = other.to_cartesian()
        except Exception:
            return NotImplemented

        first, second = (self, other_c) if not reverse else (other_c, self)
        return self.__class__(
            *(
                op(getattr(first, component), getattr(second, component))
                for component in first.components
            )
        )

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Note that any associated differentials will be dropped during this
        operation.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        # erfa pm: Modulus of p-vector.
        return erfa_ufunc.pm(self.get_xyz(xyz_axis=-1))

    def mean(self, *args, **kwargs):
        """Vector mean.

        Returns a new CartesianRepresentation instance with the means of the
        x, y, and z components.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials("mean")
        return self._apply("mean", *args, **kwargs)

    def sum(self, *args, **kwargs):
        """Vector sum.

        Returns a new CartesianRepresentation instance with the sums of the
        x, y, and z components.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.
        """
        self._raise_if_has_differentials("sum")
        return self._apply("sum", *args, **kwargs)

    def dot(self, other):
        """Dot product of two representations.

        Note that any associated differentials will be dropped during this
        operation.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            If not already cartesian, it is converted.

        Returns
        -------
        dot_product : `~astropy.units.Quantity`
            The sum of the product of the x, y, and z components of ``self``
            and ``other``.
        """
        try:
            other_c = other.to_cartesian()
        except Exception as err:
            raise TypeError(
                "can only take dot product with another "
                f"representation, not a {type(other)} instance."
            ) from err
        # erfa pdp: p-vector inner (=scalar=dot) product.
        return erfa_ufunc.pdp(self.get_xyz(xyz_axis=-1), other_c.get_xyz(xyz_axis=-1))

    def cross(self, other):
        """Cross product of two representations.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            If not already cartesian, it is converted.

        Returns
        -------
        cross_product : `~astropy.coordinates.CartesianRepresentation`
            With vectors perpendicular to both ``self`` and ``other``.
        """
        self._raise_if_has_differentials("cross")
        try:
            other_c = other.to_cartesian()
        except Exception as err:
            raise TypeError(
                "cannot only take cross product with another "
                f"representation, not a {type(other)} instance."
            ) from err
        # erfa pxp: p-vector outer (=vector=cross) product.
        sxo = erfa_ufunc.pxp(self.get_xyz(xyz_axis=-1), other_c.get_xyz(xyz_axis=-1))
        return self.__class__(sxo, xyz_axis=-1)


class CartesianDifferential(BaseDifferential):
    """Differentials in of points in 3D cartesian coordinates.

    Parameters
    ----------
    d_x, d_y, d_z : `~astropy.units.Quantity` or array
        The x, y, and z coordinates of the differentials. If ``d_x``, ``d_y``,
        and ``d_z`` have different shapes, they should be broadcastable. If not
        quantities, ``unit`` should be set.  If only ``d_x`` is given, it is
        assumed that it contains an array with the 3 coordinates stored along
        ``xyz_axis``.
    unit : `~astropy.units.Unit` or str
        If given, the differentials will be converted to this unit (or taken to
        be in this unit if not given.
    xyz_axis : int, optional
        The axis along which the coordinates are stored when a single array is
        provided instead of distinct ``d_x``, ``d_y``, and ``d_z`` (default: 0).
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.
    """

    base_representation = CartesianRepresentation
    _d_xyz = None

    def __init__(self, d_x, d_y=None, d_z=None, unit=None, xyz_axis=None, copy=True):
        if d_y is None and d_z is None:
            if isinstance(d_x, np.ndarray) and d_x.dtype.kind not in "OV":
                # Short-cut for 3-D array input.
                d_x = u.Quantity(d_x, unit, copy=copy, subok=True)
                # Keep a link to the array with all three coordinates
                # so that we can return it quickly if needed in get_xyz.
                self._d_xyz = d_x
                if xyz_axis:
                    d_x = np.moveaxis(d_x, xyz_axis, 0)
                    self._xyz_axis = xyz_axis
                else:
                    self._xyz_axis = 0

                self._d_x, self._d_y, self._d_z = d_x
                return

            else:
                d_x, d_y, d_z = d_x

        if xyz_axis is not None:
            raise ValueError(
                "xyz_axis should only be set if d_x, d_y, and d_z are in a single array"
                " passed in through d_x, i.e., d_y and d_z should not be not given."
            )

        if d_y is None or d_z is None:
            raise ValueError(
                "d_x, d_y, and d_z are required to instantiate"
                f" {self.__class__.__name__}"
            )

        if unit is not None:
            d_x = u.Quantity(d_x, unit, copy=copy, subok=True)
            d_y = u.Quantity(d_y, unit, copy=copy, subok=True)
            d_z = u.Quantity(d_z, unit, copy=copy, subok=True)
            copy = False

        super().__init__(d_x, d_y, d_z, copy=copy)
        if not (
            self._d_x.unit.is_equivalent(self._d_y.unit)
            and self._d_x.unit.is_equivalent(self._d_z.unit)
        ):
            raise u.UnitsError("d_x, d_y and d_z should have equivalent units.")

    def to_cartesian(self, base=None):
        return CartesianRepresentation(*[getattr(self, c) for c in self.components])

    @classmethod
    def from_cartesian(cls, other, base=None):
        return cls(*[getattr(other, c) for c in other.components])

    def transform(self, matrix, base=None, transformed_base=None):
        """Transform differentials using a 3x3 matrix in a Cartesian basis.

        This returns a new differential and does not modify the original one.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 (or stack thereof) matrix, such as a rotation matrix.
        base, transformed_base : `~astropy.coordinates.CartesianRepresentation` or None, optional
            Not used in the Cartesian transformation.
        """
        # erfa rxp: Multiply a p-vector by an r-matrix.
        p = erfa_ufunc.rxp(matrix, self.get_d_xyz(xyz_axis=-1))

        return self.__class__(p, xyz_axis=-1, copy=False)

    def get_d_xyz(self, xyz_axis=0):
        """Return a vector array of the x, y, and z coordinates.

        Parameters
        ----------
        xyz_axis : int, optional
            The axis in the final array along which the x, y, z components
            should be stored (default: 0).

        Returns
        -------
        d_xyz : `~astropy.units.Quantity`
            With dimension 3 along ``xyz_axis``.  Note that, if possible,
            this will be a view.
        """
        if self._d_xyz is not None:
            if self._xyz_axis == xyz_axis:
                return self._d_xyz
            else:
                return np.moveaxis(self._d_xyz, self._xyz_axis, xyz_axis)

        # Create combined array.  TO DO: keep it in _d_xyz for repeated use?
        # But then in-place changes have to cancel it. Likely best to
        # also update components.
        return np.stack([self._d_x, self._d_y, self._d_z], axis=xyz_axis)

    d_xyz = property(get_d_xyz)
