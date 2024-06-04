# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Affine coordinate transformations."""

from abc import abstractmethod

import numpy as np

from astropy.coordinates.transformations.base import CoordinateTransform

__all__ = [
    "BaseAffineTransform",
    "AffineTransform",
    "StaticMatrixTransform",
    "DynamicMatrixTransform",
]


class BaseAffineTransform(CoordinateTransform):
    """Base class for common functionality between the ``AffineTransform``-type
    subclasses.

    This base class is needed because `~astropy.coordinates.AffineTransform`
    and the matrix transform classes share the ``__call__()`` method, but
    differ in how they generate the affine parameters.
    `~astropy.coordinates.StaticMatrixTransform` passes in a matrix stored as a
    class attribute, and both of the matrix transforms pass in ``None`` for the
    offset. Hence, user subclasses would likely want to subclass this (rather
    than `~astropy.coordinates.AffineTransform`) if they want to provide
    alternative transformations using this machinery.

    """

    def _apply_transform(self, fromcoord, matrix, offset):
        from astropy.coordinates.representation import (
            CartesianDifferential,
            RadialDifferential,
            SphericalCosLatDifferential,
            SphericalDifferential,
            UnitSphericalRepresentation,
        )

        data = fromcoord.data
        has_velocity = "s" in data.differentials

        # Bail out if no transform is actually requested
        if matrix is None and offset is None:
            return data

        # list of unit differentials
        _unit_diffs = (
            SphericalDifferential._unit_differential,
            SphericalCosLatDifferential._unit_differential,
        )
        unit_vel_diff = has_velocity and isinstance(
            data.differentials["s"], _unit_diffs
        )
        rad_vel_diff = has_velocity and isinstance(
            data.differentials["s"], RadialDifferential
        )

        # Some initial checking to short-circuit doing any re-representation if
        # we're going to fail anyways:
        if isinstance(data, UnitSphericalRepresentation) and offset is not None:
            raise TypeError(
                "Position information stored on coordinate frame "
                "is insufficient to do a full-space position "
                "transformation (representation class: {data.__class__})"
            )

        if (
            has_velocity
            and (unit_vel_diff or rad_vel_diff)
            and offset is not None
            and "s" in offset.differentials
        ):
            # Coordinate has a velocity, but it is not a full-space velocity
            # that we need to do a velocity offset
            raise TypeError(
                "Velocity information stored on coordinate frame is insufficient to do"
                " a full-space velocity transformation (differential class:"
                f" {data.differentials['s'].__class__})"
            )

        if len(data.differentials) > 1:
            # We should never get here because the frame initializer shouldn't
            # allow more differentials, but this just adds protection for
            # subclasses that somehow skip the checks
            raise ValueError(
                "Representation passed to AffineTransform contains multiple associated"
                " differentials. Only a single differential with velocity units is"
                f" presently supported (differentials: {data.differentials})."
            )

        # If the representation is a UnitSphericalRepresentation, and this is
        # just a MatrixTransform, we have to try to turn the differential into a
        # Unit version of the differential (if no radial velocity) or a
        # sphericaldifferential with zero proper motion (if only a radial
        # velocity) so that the matrix operation works
        if (
            has_velocity
            and isinstance(data, UnitSphericalRepresentation)
            and not (unit_vel_diff or rad_vel_diff)
        ):
            # retrieve just velocity differential
            unit_diff = data.differentials["s"].represent_as(
                data.differentials["s"]._unit_differential, data
            )
            data = data.with_differentials({"s": unit_diff})  # updates key

        # If it's a RadialDifferential, we flat-out ignore the differentials
        # This is because, by this point (past the validation above), we can
        # only possibly be doing a rotation-only transformation, and that
        # won't change the radial differential. We later add it back in
        elif rad_vel_diff:
            data = data.without_differentials()

        # Convert the representation and differentials to cartesian without
        # having them attached to a frame
        rep = data.to_cartesian()
        diffs = {
            k: diff.represent_as(CartesianDifferential, data)
            for k, diff in data.differentials.items()
        }
        rep = rep.with_differentials(diffs)

        # Only do transform if matrix is specified. This is for speed in
        # transformations that only specify an offset (e.g., LSR)
        if matrix is not None:
            # Note: this applies to both representation and differentials
            rep = rep.transform(matrix)

        # TODO: if we decide to allow arithmetic between representations that
        # contain differentials, this can be tidied up
        newrep = rep.without_differentials()
        if offset is not None:
            newrep += offset.without_differentials()

        # We need a velocity (time derivative) and, for now, are strict: the
        # representation can only contain a velocity differential and no others.
        if has_velocity and not rad_vel_diff:
            veldiff = rep.differentials["s"]  # already in Cartesian form

            if offset is not None and "s" in offset.differentials:
                veldiff += offset.differentials["s"]

            newrep = newrep.with_differentials({"s": veldiff})

        if isinstance(fromcoord.data, UnitSphericalRepresentation):
            # Special-case this because otherwise the return object will think
            # it has a valid distance with the default return (a
            # CartesianRepresentation instance)

            if has_velocity and not unit_vel_diff and not rad_vel_diff:
                # We have to first represent as the Unit types we converted to,
                # then put the d_distance information back in to the
                # differentials and re-represent as their original forms
                newdiff = newrep.differentials["s"]
                _unit_cls = fromcoord.data.differentials["s"]._unit_differential
                newdiff = newdiff.represent_as(_unit_cls, newrep)

                kwargs = {comp: getattr(newdiff, comp) for comp in newdiff.components}
                kwargs["d_distance"] = fromcoord.data.differentials["s"].d_distance
                diffs = {
                    "s": type(fromcoord.data.differentials["s"])(copy=False, **kwargs)
                }

            elif has_velocity and unit_vel_diff:
                newdiff = newrep.differentials["s"].represent_as(
                    fromcoord.data.differentials["s"].__class__, newrep
                )
                diffs = {"s": newdiff}

            else:
                diffs = newrep.differentials

            newrep = newrep.represent_as(type(fromcoord.data)).with_differentials(diffs)

        elif has_velocity and unit_vel_diff:
            # Here, we're in the case where the representation is not
            # UnitSpherical, but the differential *is* one of the UnitSpherical
            # types. We have to convert back to that differential class or the
            # resulting frame will think it has a valid radial_velocity. This
            # can probably be cleaned up: we currently have to go through the
            # dimensional version of the differential before representing as the
            # unit differential so that the units work out (the distance length
            # unit shouldn't appear in the resulting proper motions)

            diff_cls = fromcoord.data.differentials["s"].__class__
            newrep = newrep.represent_as(
                type(fromcoord.data), diff_cls._dimensional_differential
            ).represent_as(type(fromcoord.data), diff_cls)

        # We pulled the radial differential off of the representation
        # earlier, so now we need to put it back. But, in order to do that, we
        # have to turn the representation into a repr that is compatible with
        # having a RadialDifferential
        if has_velocity and rad_vel_diff:
            newrep = newrep.represent_as(fromcoord.data.__class__)
            newrep = newrep.with_differentials({"s": fromcoord.data.differentials["s"]})

        return newrep

    def __call__(self, fromcoord, toframe):
        params = self._affine_params(fromcoord, toframe)
        newrep = self._apply_transform(fromcoord, *params)
        return toframe.realize_frame(newrep)

    @abstractmethod
    def _affine_params(self, fromcoord, toframe):
        pass


class AffineTransform(BaseAffineTransform):
    """
    A coordinate transformation specified as a function that yields a 3 x 3
    cartesian transformation matrix and a tuple of displacement vectors.

    See `~astropy.coordinates.Galactocentric` for
    an example.

    Parameters
    ----------
    transform_func : callable
        A callable that has the signature ``transform_func(fromcoord, toframe)``
        and returns: a (3, 3) matrix that operates on ``fromcoord`` in a
        Cartesian representation, and a ``CartesianRepresentation`` with
        (optionally) an attached velocity ``CartesianDifferential`` to represent
        a translation and offset in velocity to apply after the matrix
        operation.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `~astropy.coordinates.TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``transform_func`` is not callable

    """

    def __init__(self, transform_func, fromsys, tosys, priority=1, register_graph=None):
        if not callable(transform_func):
            raise TypeError("transform_func is not callable")
        self.transform_func = transform_func

        super().__init__(
            fromsys, tosys, priority=priority, register_graph=register_graph
        )

    def _affine_params(self, fromcoord, toframe):
        return self.transform_func(fromcoord, toframe)


class StaticMatrixTransform(BaseAffineTransform):
    """
    A coordinate transformation defined as a 3 x 3 cartesian
    transformation matrix.

    This is distinct from DynamicMatrixTransform in that this kind of matrix is
    independent of frame attributes.  That is, it depends *only* on the class of
    the frame.

    Parameters
    ----------
    matrix : array-like or callable
        A 3 x 3 matrix for transforming 3-vectors. In most cases will
        be unitary (although this is not strictly required). If a callable,
        will be called *with no arguments* to get the matrix.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `~astropy.coordinates.TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    ValueError
        If the matrix is not 3 x 3

    """

    def __init__(self, matrix, fromsys, tosys, priority=1, register_graph=None):
        if callable(matrix):
            matrix = matrix()
        self.matrix = np.array(matrix)

        if self.matrix.shape != (3, 3):
            raise ValueError("Provided matrix is not 3 x 3")

        super().__init__(
            fromsys, tosys, priority=priority, register_graph=register_graph
        )

    def _affine_params(self, fromcoord, toframe):
        return self.matrix, None


class DynamicMatrixTransform(BaseAffineTransform):
    """
    A coordinate transformation specified as a function that yields a
    3 x 3 cartesian transformation matrix.

    This is similar to, but distinct from StaticMatrixTransform, in that the
    matrix for this class might depend on frame attributes.

    Parameters
    ----------
    matrix_func : callable
        A callable that has the signature ``matrix_func(fromcoord, toframe)`` and
        returns a 3 x 3 matrix that converts ``fromcoord`` in a cartesian
        representation to the new coordinate system.
    fromsys : class
        The coordinate frame class to start from.
    tosys : class
        The coordinate frame class to transform into.
    priority : float or int
        The priority if this transform when finding the shortest
        coordinate transform path - large numbers are lower priorities.
    register_graph : `~astropy.coordinates.TransformGraph` or None
        A graph to register this transformation with on creation, or
        `None` to leave it unregistered.

    Raises
    ------
    TypeError
        If ``matrix_func`` is not callable

    """

    def __init__(self, matrix_func, fromsys, tosys, priority=1, register_graph=None):
        if not callable(matrix_func):
            raise TypeError("matrix_func is not callable")
        self.matrix_func = matrix_func

        super().__init__(
            fromsys, tosys, priority=priority, register_graph=register_graph
        )

    def _affine_params(self, fromcoord, toframe):
        return self.matrix_func(fromcoord, toframe), None
