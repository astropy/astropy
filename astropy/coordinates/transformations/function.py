# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Function-based coordinate transformations.

These are transformations that cannot be represented as an affine transformation.
"""

from contextlib import suppress
from inspect import signature
from warnings import warn

from astropy import units as u
from astropy.coordinates.transformations.base import CoordinateTransform
from astropy.utils.exceptions import AstropyWarning

__all__ = ["FunctionTransform", "FunctionTransformWithFiniteDifference"]


class FunctionTransform(CoordinateTransform):
    """
    A coordinate transformation defined by a function that accepts a
    coordinate object and returns the transformed coordinate object.

    Parameters
    ----------
    func : callable
        The transformation function. Should have a call signature
        ``func(formcoord, toframe)``. Note that, unlike
        `CoordinateTransform.__call__`, ``toframe`` is assumed to be of type
        ``tosys`` for this function.
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
        If ``func`` is not callable.
    ValueError
        If ``func`` cannot accept two arguments.


    """

    def __init__(self, func, fromsys, tosys, priority=1, register_graph=None):
        if not callable(func):
            raise TypeError("func must be callable")

        with suppress(TypeError):
            sig = signature(func)
            kinds = [x.kind for x in sig.parameters.values()]
            if (
                len(x for x in kinds if x == sig.POSITIONAL_ONLY) != 2
                and sig.VAR_POSITIONAL not in kinds
            ):
                raise ValueError("provided function does not accept two arguments")

        self.func = func

        super().__init__(
            fromsys, tosys, priority=priority, register_graph=register_graph
        )

    def __call__(self, fromcoord, toframe):
        res = self.func(fromcoord, toframe)
        if not isinstance(res, self.tosys):
            raise TypeError(
                f"the transformation function yielded {res} but "
                f"should have been of type {self.tosys}"
            )
        if fromcoord.data.differentials and not res.data.differentials:
            warn(
                "Applied a FunctionTransform to a coordinate frame with "
                "differentials, but the FunctionTransform does not handle "
                "differentials, so they have been dropped.",
                AstropyWarning,
            )
        return res


class FunctionTransformWithFiniteDifference(FunctionTransform):
    r"""Transformation based on functions using finite difference for velocities.

    A coordinate transformation that works like a
    `~astropy.coordinates.FunctionTransform`, but computes velocity shifts
    based on the finite-difference relative to one of the frame attributes.
    Note that the transform function should *not* change the differential at
    all in this case, as any differentials will be overridden.

    When a differential is in the from coordinate, the finite difference
    calculation has two components. The first part is simple the existing
    differential, but re-orientation (using finite-difference techniques) to
    point in the direction the velocity vector has in the *new* frame. The
    second component is the "induced" velocity.  That is, the velocity
    intrinsic to the frame itself, estimated by shifting the frame using the
    ``finite_difference_frameattr_name`` frame attribute a small amount
    (``finite_difference_dt``) in time and re-calculating the position.

    Parameters
    ----------
    finite_difference_frameattr_name : str or None
        The name of the frame attribute on the frames to use for the finite
        difference.  Both the to and the from frame will be checked for this
        attribute, but only one needs to have it. If None, no velocity
        component induced from the frame itself will be included - only the
        re-orientation of any existing differential.
    finite_difference_dt : `~astropy.units.Quantity` ['time'] or callable
        If a quantity, this is the size of the differential used to do the
        finite difference.  If a callable, should accept
        ``(fromcoord, toframe)`` and return the ``dt`` value.
    symmetric_finite_difference : bool
        If True, the finite difference is computed as
        :math:`\frac{x(t + \Delta t / 2) - x(t + \Delta t / 2)}{\Delta t}`, or
        if False, :math:`\frac{x(t + \Delta t) - x(t)}{\Delta t}`.  The latter
        case has slightly better performance (and more stable finite difference
        behavior).

    All other parameters are identical to the initializer for
    `~astropy.coordinates.FunctionTransform`.

    """

    def __init__(
        self,
        func,
        fromsys,
        tosys,
        priority=1,
        register_graph=None,
        finite_difference_frameattr_name="obstime",
        finite_difference_dt=1 * u.second,
        symmetric_finite_difference=True,
    ):
        super().__init__(func, fromsys, tosys, priority, register_graph)
        self.finite_difference_frameattr_name = finite_difference_frameattr_name
        self.finite_difference_dt = finite_difference_dt
        self.symmetric_finite_difference = symmetric_finite_difference

    @property
    def finite_difference_frameattr_name(self):
        return self._finite_difference_frameattr_name

    @finite_difference_frameattr_name.setter
    def finite_difference_frameattr_name(self, value):
        if value is None:
            self._diff_attr_in_fromsys = self._diff_attr_in_tosys = False
        else:
            diff_attr_in_fromsys = value in self.fromsys.frame_attributes
            diff_attr_in_tosys = value in self.tosys.frame_attributes
            if diff_attr_in_fromsys or diff_attr_in_tosys:
                self._diff_attr_in_fromsys = diff_attr_in_fromsys
                self._diff_attr_in_tosys = diff_attr_in_tosys
            else:
                raise ValueError(
                    f"Frame attribute name {value} is not a frame attribute of"
                    f" {self.fromsys} or {self.tosys}"
                )
        self._finite_difference_frameattr_name = value

    def __call__(self, fromcoord, toframe):
        from astropy.coordinates.representation import (
            CartesianDifferential,
            CartesianRepresentation,
        )

        supcall = self.func
        if not fromcoord.data.differentials:
            return supcall(fromcoord, toframe)
        # this is the finite difference case

        if callable(self.finite_difference_dt):
            dt = self.finite_difference_dt(fromcoord, toframe)
        else:
            dt = self.finite_difference_dt
        halfdt = dt / 2

        from_diffless = fromcoord.realize_frame(fromcoord.data.without_differentials())
        reprwithoutdiff = supcall(from_diffless, toframe)

        # first we use the existing differential to compute an offset due to
        # the already-existing velocity, but in the new frame
        fromcoord_cart = fromcoord.cartesian
        if self.symmetric_finite_difference:
            fwdxyz = (
                fromcoord_cart.xyz + fromcoord_cart.differentials["s"].d_xyz * halfdt
            )
            fwd = supcall(
                fromcoord.realize_frame(CartesianRepresentation(fwdxyz)), toframe
            )
            backxyz = (
                fromcoord_cart.xyz - fromcoord_cart.differentials["s"].d_xyz * halfdt
            )
            back = supcall(
                fromcoord.realize_frame(CartesianRepresentation(backxyz)), toframe
            )
        else:
            fwdxyz = fromcoord_cart.xyz + fromcoord_cart.differentials["s"].d_xyz * dt
            fwd = supcall(
                fromcoord.realize_frame(CartesianRepresentation(fwdxyz)), toframe
            )
            back = reprwithoutdiff
        diffxyz = (fwd.cartesian - back.cartesian).xyz / dt

        # now we compute the "induced" velocities due to any movement in
        # the frame itself over time
        attrname = self.finite_difference_frameattr_name
        if attrname is not None:
            if self.symmetric_finite_difference:
                if self._diff_attr_in_fromsys:
                    kws = {attrname: getattr(from_diffless, attrname) + halfdt}
                    from_diffless_fwd = from_diffless.replicate(**kws)
                else:
                    from_diffless_fwd = from_diffless
                if self._diff_attr_in_tosys:
                    kws = {attrname: getattr(toframe, attrname) + halfdt}
                    fwd_frame = toframe.replicate_without_data(**kws)
                else:
                    fwd_frame = toframe
                fwd = supcall(from_diffless_fwd, fwd_frame)

                if self._diff_attr_in_fromsys:
                    kws = {attrname: getattr(from_diffless, attrname) - halfdt}
                    from_diffless_back = from_diffless.replicate(**kws)
                else:
                    from_diffless_back = from_diffless
                if self._diff_attr_in_tosys:
                    kws = {attrname: getattr(toframe, attrname) - halfdt}
                    back_frame = toframe.replicate_without_data(**kws)
                else:
                    back_frame = toframe
                back = supcall(from_diffless_back, back_frame)
            else:
                if self._diff_attr_in_fromsys:
                    kws = {attrname: getattr(from_diffless, attrname) + dt}
                    from_diffless_fwd = from_diffless.replicate(**kws)
                else:
                    from_diffless_fwd = from_diffless
                if self._diff_attr_in_tosys:
                    kws = {attrname: getattr(toframe, attrname) + dt}
                    fwd_frame = toframe.replicate_without_data(**kws)
                else:
                    fwd_frame = toframe
                fwd = supcall(from_diffless_fwd, fwd_frame)
                back = reprwithoutdiff

            diffxyz += (fwd.cartesian - back.cartesian).xyz / dt

        newdiff = CartesianDifferential(diffxyz)
        reprwithdiff = reprwithoutdiff.data.to_cartesian().with_differentials(newdiff)
        return reprwithoutdiff.realize_frame(reprwithdiff)
