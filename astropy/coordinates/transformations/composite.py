# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Module defining the class for a composited sequence of coordinate transformations.

This module contains the class :class:`~astropy.coordinates.CompositeTransform` that
combines a sequence of transformations into a single transformation. The class has the
same API as a single-step transformation, so it can be used interchangeably with a
single-step transformation.
"""

from astropy.coordinates.transformations.affine import (
    AffineTransform,
    BaseAffineTransform,
    DynamicMatrixTransform,
    StaticMatrixTransform,
)
from astropy.coordinates.transformations.base import CoordinateTransform
from astropy.coordinates.transformations.function import (
    FunctionTransformWithFiniteDifference,
)

__all__ = ["CompositeTransform"]


class CompositeTransform(CoordinateTransform):
    """
    A transformation constructed by combining together a series of single-step
    transformations.

    Note that the intermediate frame objects are constructed using any frame
    attributes in ``toframe`` or ``fromframe`` that overlap with the intermediate
    frame (``toframe`` favored over ``fromframe`` if there's a conflict).  Any frame
    attributes that are not present use the defaults.

    Parameters
    ----------
    transforms : sequence of `~astropy.coordinates.CoordinateTransform` object
        The sequence of transformations to apply.
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
    collapse_static_mats : bool
        If `True`, consecutive `~astropy.coordinates.StaticMatrixTransform`
        will be collapsed into a single transformation to speed up the
        calculation.

    """

    def __init__(
        self,
        transforms,
        fromsys,
        tosys,
        priority=1,
        register_graph=None,
        collapse_static_mats=True,
    ):
        super().__init__(
            fromsys, tosys, priority=priority, register_graph=register_graph
        )

        if collapse_static_mats:
            transforms = self._combine_statics(transforms)

        self.transforms = tuple(transforms)

    def _combine_statics(self, transforms):
        """
        Combines together sequences of StaticMatrixTransform's into a single
        transform and returns it.
        """
        newtrans = []
        for currtrans in transforms:
            lasttrans = newtrans[-1] if len(newtrans) > 0 else None

            if isinstance(lasttrans, StaticMatrixTransform) and isinstance(
                currtrans, StaticMatrixTransform
            ):
                newtrans[-1] = StaticMatrixTransform(
                    currtrans.matrix @ lasttrans.matrix,
                    lasttrans.fromsys,
                    currtrans.tosys,
                )
            else:
                newtrans.append(currtrans)
        return newtrans

    def __call__(self, fromcoord, toframe):
        curr_coord = fromcoord
        for t in self.transforms:
            # build an intermediate frame with attributes taken from either
            # `toframe`, or if not there, `fromcoord`, or if not there, use
            # the defaults
            # TODO: caching this information when creating the transform may
            # speed things up a lot
            frattrs = {}
            for inter_frame_attr_nm in t.tosys.frame_attributes:
                if hasattr(toframe, inter_frame_attr_nm):
                    attr = getattr(toframe, inter_frame_attr_nm)
                    frattrs[inter_frame_attr_nm] = attr
                elif hasattr(fromcoord, inter_frame_attr_nm):
                    attr = getattr(fromcoord, inter_frame_attr_nm)
                    frattrs[inter_frame_attr_nm] = attr

            curr_toframe = t.tosys(**frattrs)
            curr_coord = t(curr_coord, curr_toframe)

        # this is safe even in the case where self.transforms is empty, because
        # coordinate objects are immutable, so copying is not needed
        return curr_coord

    def _as_single_transform(self):
        """
        Return an encapsulated version of the composite transform so that it appears to
        be a single transform.

        The returned transform internally calls the constituent transforms.  If all of
        the transforms are affine, the merged transform is
        `~astropy.coordinates.DynamicMatrixTransform` (if there are no
        origin shifts) or `~astropy.coordinates.AffineTransform`
        (otherwise).  If at least one of the transforms is not affine, the merged
        transform is
        `~astropy.coordinates.FunctionTransformWithFiniteDifference`.
        """
        # Create a list of the transforms including flattening any constituent CompositeTransform
        transforms = [
            t if not isinstance(t, CompositeTransform) else t._as_single_transform()
            for t in self.transforms
        ]

        if all(isinstance(t, BaseAffineTransform) for t in transforms):
            # Check if there may be an origin shift
            fixed_origin = all(
                isinstance(t, (StaticMatrixTransform, DynamicMatrixTransform))
                for t in transforms
            )

            # Dynamically define the transformation function
            def single_transform(from_coo, to_frame):
                if from_coo.is_equivalent_frame(to_frame):  # loopback to the same frame
                    return None if fixed_origin else (None, None)

                # Create a merged attribute dictionary for any intermediate frames
                # For any attributes shared by the "from"/"to" frames, the "to" frame takes
                #   precedence because this is the same choice implemented in __call__()
                merged_attr = {
                    name: getattr(from_coo, name) for name in from_coo.frame_attributes
                }
                merged_attr.update(
                    {
                        name: getattr(to_frame, name)
                        for name in to_frame.frame_attributes
                    }
                )

                affine_params = (None, None)
                # Step through each transform step (frame A -> frame B)
                for i, t in enumerate(transforms):
                    # Extract the relevant attributes for frame A
                    if i == 0:
                        # If frame A is actually the initial frame, preserve its attributes
                        a_attr = {
                            name: getattr(from_coo, name)
                            for name in from_coo.frame_attributes
                        }
                    else:
                        a_attr = {
                            k: v
                            for k, v in merged_attr.items()
                            if k in t.fromsys.frame_attributes
                        }

                    # Extract the relevant attributes for frame B
                    b_attr = {
                        k: v
                        for k, v in merged_attr.items()
                        if k in t.tosys.frame_attributes
                    }

                    # Obtain the affine parameters for the transform
                    # Note that we insert some dummy data into frame A because the transformation
                    #   machinery requires there to be data present.  Removing that limitation
                    #   is a possible TODO, but some care would need to be taken because some affine
                    #   transforms have branching code depending on the presence of differentials.
                    next_affine_params = t._affine_params(
                        t.fromsys(from_coo.data, **a_attr), t.tosys(**b_attr)
                    )

                    # Combine the affine parameters with the running set
                    affine_params = _combine_affine_params(
                        affine_params, next_affine_params
                    )

                # If there is no origin shift, return only the matrix
                return affine_params[0] if fixed_origin else affine_params

            # The return type depends on whether there is any origin shift
            transform_type = DynamicMatrixTransform if fixed_origin else AffineTransform
        else:
            # Dynamically define the transformation function
            def single_transform(from_coo, to_frame):
                if from_coo.is_equivalent_frame(to_frame):  # loopback to the same frame
                    return to_frame.realize_frame(from_coo.data)
                return self(from_coo, to_frame)

            transform_type = FunctionTransformWithFiniteDifference

        return transform_type(
            single_transform, self.fromsys, self.tosys, priority=self.priority
        )


def _combine_affine_params(params, next_params):
    """
    Combine two sets of affine parameters.

    The parameters for an affine transformation are a 3 x 3 Cartesian
    transformation matrix and a displacement vector, which can include an
    attached velocity.  Either type of parameter can be ``None``.
    """
    M, vec = params
    next_M, next_vec = next_params

    # Multiply the transformation matrices if they both exist
    if M is not None and next_M is not None:
        new_M = next_M @ M
    else:
        new_M = M if M is not None else next_M

    if vec is not None:
        # Transform the first displacement vector by the second transformation matrix
        if next_M is not None:
            vec = vec.transform(next_M)

        # Calculate the new displacement vector
        if next_vec is not None:
            if "s" in vec.differentials and "s" in next_vec.differentials:
                # Adding vectors with velocities takes more steps
                # TODO: Add support in representation.py
                new_vec_velocity = vec.differentials["s"] + next_vec.differentials["s"]
                new_vec = vec.without_differentials() + next_vec.without_differentials()
                new_vec = new_vec.with_differentials({"s": new_vec_velocity})
            else:
                new_vec = vec + next_vec
        else:
            new_vec = vec
    else:
        new_vec = next_vec

    return new_M, new_vec
