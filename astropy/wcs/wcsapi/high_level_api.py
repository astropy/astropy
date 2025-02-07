import abc
import numbers
from collections import OrderedDict, defaultdict

import numpy as np

from .utils import deserialize_class

__all__ = [
    "BaseHighLevelWCS",
    "HighLevelWCSMixin",
    "high_level_objects_to_values",
    "values_to_high_level_objects",
]


def rec_getattr(obj, att):
    for a in att.split("."):
        obj = getattr(obj, a)
    return obj


def default_order(components):
    order = []
    for key, _, _ in components:
        if key not in order:
            order.append(key)
    return order


def _toindex(value):
    """Convert value to an int or an int array.

    Input coordinates converted to integers
    corresponding to the center of the pixel.
    The convention is that the center of the pixel is
    (0, 0), while the lower left corner is (-0.5, -0.5).
    The outputs are used to index the mask.

    Examples
    --------
    >>> _toindex(np.array([-0.5, 0.49999]))
    array([0, 0])
    >>> _toindex(np.array([0.5, 1.49999]))
    array([1, 1])
    >>> _toindex(np.array([1.5, 2.49999]))
    array([2, 2])
    """
    arr = np.floor(np.asarray(value) + 0.5)

    fill_value = np.iinfo(int).min
    if np.isscalar(arr):
        if np.isnan(arr):
            arr = fill_value
    else:
        arr[np.isnan(arr)] = fill_value
    return np.asarray(arr, dtype=int)


class BaseHighLevelWCS(metaclass=abc.ABCMeta):
    """
    Abstract base class for the high-level WCS interface.

    This is described in `APE 14: A shared Python interface for World Coordinate
    Systems <https://doi.org/10.5281/zenodo.1188875>`_.
    """

    @property
    @abc.abstractmethod
    def low_level_wcs(self):
        """
        Returns a reference to the underlying low-level WCS object.
        """

    @abc.abstractmethod
    def pixel_to_world(self, *pixel_arrays):
        """
        Convert pixel coordinates to world coordinates (represented by
        high-level objects).

        If a single high-level object is used to represent the world coordinates
        (i.e., if ``len(wcs.world_axis_object_classes) == 1``), it is returned
        as-is (not in a tuple/list), otherwise a tuple of high-level objects is
        returned. See
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values` for pixel
        indexing and ordering conventions.
        """

    def array_index_to_world(self, *index_arrays):
        """
        Convert array indices to world coordinates (represented by Astropy
        objects).

        If a single high-level object is used to represent the world coordinates
        (i.e., if ``len(wcs.world_axis_object_classes) == 1``), it is returned
        as-is (not in a tuple/list), otherwise a tuple of high-level objects is
        returned. See
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.array_index_to_world_values` for
        pixel indexing and ordering conventions.
        """
        return self.pixel_to_world(*index_arrays[::-1])

    @abc.abstractmethod
    def world_to_pixel(self, *world_objects):
        """
        Convert world coordinates (represented by Astropy objects) to pixel
        coordinates.

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned. See
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_to_pixel_values` for pixel
        indexing and ordering conventions.
        """

    def world_to_array_index(self, *world_objects):
        """
        Convert world coordinates (represented by Astropy objects) to array
        indices.

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned. See
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_to_array_index_values` for
        pixel indexing and ordering conventions. The indices should be returned
        as rounded integers.
        """
        if self.low_level_wcs.pixel_n_dim == 1:
            return _toindex(self.world_to_pixel(*world_objects))
        else:
            return tuple(_toindex(self.world_to_pixel(*world_objects)[::-1]).tolist())


def high_level_objects_to_values(*world_objects, low_level_wcs):
    """
    Convert the input high level object to low level values.

    This function uses the information in ``wcs.world_axis_object_classes`` and
    ``wcs.world_axis_object_components`` to convert the high level objects
    (such as `~.SkyCoord`) to low level "values" which should be scalars or
    Numpy arrays.

    This is used in `.HighLevelWCSMixin.world_to_pixel`, but provided as a
    separate function for use in other places where needed.

    Parameters
    ----------
    *world_objects: object
        High level coordinate objects.

    low_level_wcs: `.BaseLowLevelWCS`
        The WCS object to use to interpret the coordinates.
    """
    # Cache the classes and components since this may be expensive
    serialized_classes = low_level_wcs.world_axis_object_classes
    components = low_level_wcs.world_axis_object_components

    # Deserialize world_axis_object_classes using the default order
    classes = OrderedDict()
    for key in default_order(components):
        if low_level_wcs.serialized_classes:
            classes[key] = deserialize_class(serialized_classes[key], construct=False)
        else:
            classes[key] = serialized_classes[key]

    # Check that the number of classes matches the number of inputs
    if len(world_objects) != len(classes):
        raise ValueError(
            f"Number of world inputs ({len(world_objects)}) does not match expected"
            f" ({len(classes)})"
        )

    # Determine whether the classes are uniquely matched, that is we check
    # whether there is only one of each class.
    world_by_key = {}
    unique_match = True
    for w in world_objects:
        matches = []
        for key, (klass, *_) in classes.items():
            if isinstance(w, klass):
                matches.append(key)
        if len(matches) == 1:
            world_by_key[matches[0]] = w
        else:
            unique_match = False
            break

    # If the match is not unique, the order of the classes needs to match,
    # whereas if all classes are unique, we can still intelligently match
    # them even if the order is wrong.

    objects = {}

    if unique_match:
        for key, (klass, args, kwargs, *rest) in classes.items():
            if len(rest) == 0:
                klass_gen = klass
            elif len(rest) == 1:
                klass_gen = rest[0]
            else:
                raise ValueError(
                    "Tuples in world_axis_object_classes should have length 3 or 4"
                )

            # FIXME: For now SkyCoord won't auto-convert upon initialization
            # https://github.com/astropy/astropy/issues/7689
            from astropy.coordinates import SkyCoord

            if isinstance(world_by_key[key], SkyCoord):
                if "frame" in kwargs:
                    objects[key] = world_by_key[key].transform_to(kwargs["frame"])
                else:
                    objects[key] = world_by_key[key]
            else:
                objects[key] = klass_gen(world_by_key[key], *args, **kwargs)

    else:
        for ikey, key in enumerate(classes):
            klass, args, kwargs, *rest = classes[key]

            if len(rest) == 0:
                klass_gen = klass
            elif len(rest) == 1:
                klass_gen = rest[0]
            else:
                raise ValueError(
                    "Tuples in world_axis_object_classes should have length 3 or 4"
                )

            w = world_objects[ikey]
            if not isinstance(w, klass):
                raise ValueError(
                    "Expected the following order of world arguments:"
                    f" {', '.join([k.__name__ for (k, *_) in classes.values()])}"
                )

            # FIXME: For now SkyCoord won't auto-convert upon initialization
            # https://github.com/astropy/astropy/issues/7689
            from astropy.coordinates import SkyCoord

            if isinstance(w, SkyCoord):
                if "frame" in kwargs:
                    objects[key] = w.transform_to(kwargs["frame"])
                else:
                    objects[key] = w
            else:
                objects[key] = klass_gen(w, *args, **kwargs)

    # We now extract the attributes needed for the world values
    world = []
    for key, _, attr in components:
        if callable(attr):
            world.append(attr(objects[key]))
        else:
            world.append(rec_getattr(objects[key], attr))

    # Check the type of the return values - should be scalars or plain Numpy
    # arrays, not e.g. Quantity. Note that we deliberately use type(w) because
    # we don't want to match Numpy subclasses.
    for w in world:
        if not isinstance(w, numbers.Number) and not type(w) == np.ndarray:
            raise TypeError(
                f"WCS world_axis_object_components results in "
                f"values which are not scalars or plain Numpy "
                f"arrays (got {type(w)})"
            )

    return world


def values_to_high_level_objects(*world_values, low_level_wcs):
    """
    Convert low level values into high level objects.

    This function uses the information in ``wcs.world_axis_object_classes`` and
    ``wcs.world_axis_object_components`` to convert low level "values"
    `~.Quantity` objects, to high level objects (such as `~.SkyCoord`).

    This is used in `.HighLevelWCSMixin.pixel_to_world`, but provided as a
    separate function for use in other places where needed.

    Parameters
    ----------
    *world_values: object
        Low level, "values" representations of the world coordinates.

    low_level_wcs: `.BaseLowLevelWCS`
        The WCS object to use to interpret the coordinates.
    """
    # Check the type of the input values - should be scalars or plain Numpy
    # arrays, not e.g. Quantity. Note that we deliberately use type(w) because
    # we don't want to match Numpy subclasses.
    for w in world_values:
        if not isinstance(w, numbers.Number) and not type(w) == np.ndarray:
            raise TypeError(
                f"Expected world coordinates as scalars or plain Numpy "
                f"arrays (got {type(w)})"
            )

    # Cache the classes and components since this may be expensive
    components = low_level_wcs.world_axis_object_components
    classes = low_level_wcs.world_axis_object_classes

    # Deserialize classes
    if low_level_wcs.serialized_classes:
        classes_new = {}
        for key, value in classes.items():
            classes_new[key] = deserialize_class(value, construct=False)
        classes = classes_new

    args = defaultdict(list)
    kwargs = defaultdict(dict)

    for i, (key, attr, _) in enumerate(components):
        if isinstance(attr, str):
            kwargs[key][attr] = world_values[i]
        else:
            while attr > len(args[key]) - 1:
                args[key].append(None)
            args[key][attr] = world_values[i]

    result = []

    for key in default_order(components):
        klass, ar, kw, *rest = classes[key]
        if len(rest) == 0:
            klass_gen = klass
        elif len(rest) == 1:
            klass_gen = rest[0]
        else:
            raise ValueError(
                "Tuples in world_axis_object_classes should have length 3 or 4"
            )
        result.append(klass_gen(*args[key], *ar, **kwargs[key], **kw))

    return result


class HighLevelWCSMixin(BaseHighLevelWCS):
    """
    Mix-in class that automatically provides the high-level WCS API for the
    low-level WCS object given by the `~HighLevelWCSMixin.low_level_wcs`
    property.
    """

    @property
    def low_level_wcs(self):
        return self

    def world_to_pixel(self, *world_objects):
        world_values = high_level_objects_to_values(
            *world_objects, low_level_wcs=self.low_level_wcs
        )

        # Finally we convert to pixel coordinates
        pixel_values = self.low_level_wcs.world_to_pixel_values(*world_values)

        return pixel_values

    def pixel_to_world(self, *pixel_arrays):
        # Compute the world coordinate values
        world_values = self.low_level_wcs.pixel_to_world_values(*pixel_arrays)

        if self.low_level_wcs.world_n_dim == 1:
            world_values = (world_values,)

        pixel_values = values_to_high_level_objects(
            *world_values, low_level_wcs=self.low_level_wcs
        )

        if len(pixel_values) == 1:
            return pixel_values[0]
        else:
            return pixel_values
