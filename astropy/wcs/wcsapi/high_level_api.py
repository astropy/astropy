import abc
from collections import defaultdict, OrderedDict

import numpy as np

from .utils import deserialize_class

__all__ = ['BaseHighLevelWCS', 'HighLevelWCSMixin']


def rec_getattr(obj, att):
    for a in att.split('.'):
        obj = getattr(obj, a)
    return obj


def default_order(components):
    order = []
    for key, _, _ in components:
        if key not in order:
            order.append(key)
    return order


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
        Convert pixel coordinates to world coordinates (represented by high-level
        objects).

        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values` for pixel indexing and
        ordering conventions.
        """

    @abc.abstractmethod
    def array_index_to_world(self, *index_arrays):
        """
        Convert array indices to world coordinates (represented by Astropy
        objects).

        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.array_index_to_world_values` for pixel indexing and
        ordering conventions.
        """

    @abc.abstractmethod
    def world_to_pixel(self, *world_objects):
        """
        Convert world coordinates (represented by Astropy objects) to pixel
        coordinates.

        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_to_pixel_values` for pixel indexing and
        ordering conventions.
        """

    @abc.abstractmethod
    def world_to_array_index(self, *world_objects):
        """
        Convert world coordinates (represented by Astropy objects) to array
        indices.

        See `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_to_array_index_values` for pixel indexing
        and ordering conventions. The indices should be returned as rounded
        integers.
        """


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

        # Cache the classes and components since this may be expensive
        serialized_classes = self.low_level_wcs.world_axis_object_classes
        components = self.low_level_wcs.world_axis_object_components

        # Deserialize world_axis_object_classes using the default order
        classes = OrderedDict()
        for key in default_order(components):
            if self.low_level_wcs.serialized_classes:
                classes[key] = deserialize_class(serialized_classes[key],
                                                 construct=False)
            else:
                classes[key] = serialized_classes[key]

        # Check that the number of classes matches the number of inputs
        if len(world_objects) != len(classes):
            raise ValueError("Number of world inputs ({0}) does not match "
                             "expected ({1})".format(len(world_objects), len(classes)))

        # Determine whether the classes are uniquely matched, that is we check
        # whether there is only one of each class.
        world_by_key = {}
        unique_match = True
        for w in world_objects:
            matches = []
            for key, (klass, _, _) in classes.items():
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

            for key, (klass, args, kwargs) in classes.items():

                # FIXME: For now SkyCoord won't auto-convert upon initialization
                # https://github.com/astropy/astropy/issues/7689
                from ...coordinates import SkyCoord
                if isinstance(world_by_key[key], SkyCoord):
                    if 'frame' in kwargs:
                        objects[key] = world_by_key[key].transform_to(kwargs['frame'])
                    else:
                        objects[key] = world_by_key[key]
                else:
                    objects[key] = klass(world_by_key[key], *args, **kwargs)

        else:

            for ikey, key in enumerate(classes):
                klass, args, kwargs = classes[key]
                w = world_objects[ikey]
                if not isinstance(w, klass):
                    raise ValueError("Expected the following order of world "
                                     "arguments: {0}".format(', '.join([k.__name__ for (k, _, _) in classes.values()])))

                # FIXME: For now SkyCoord won't auto-convert upon initialization
                # https://github.com/astropy/astropy/issues/7689
                from ...coordinates import SkyCoord
                if isinstance(w, SkyCoord):
                    if 'frame' in kwargs:
                        objects[key] = w.transform_to(kwargs['frame'])
                    else:
                        objects[key] = w
                else:
                    objects[key] = klass(w, *args, **kwargs)

        # We now extract the attributes needed for the world values
        world = []
        for key, _, attr in components:
            world.append(rec_getattr(objects[key], attr))

        # Finally we convert to pixel coordinates
        pixel = self.low_level_wcs.world_to_pixel_values(*world)

        return pixel

    def pixel_to_world(self, *pixel_arrays):

        # Compute the world coordinate values
        world = self.low_level_wcs.pixel_to_world_values(*pixel_arrays)

        # Cache the classes and components since this may be expensive
        components = self.low_level_wcs.world_axis_object_components
        classes = self.low_level_wcs.world_axis_object_classes

        # Deserialize classes
        if self.low_level_wcs.serialized_classes:
            classes_new = {}
            for key, value in classes.items():
                classes_new[key] = deserialize_class(value, construct=False)
            classes = classes_new

        args = defaultdict(list)
        kwargs = defaultdict(dict)

        for i, (key, attr, _) in enumerate(components):
            if isinstance(attr, str):
                kwargs[key][attr] = world[i]
            else:
                while attr > len(args[key]) - 1:
                    args[key].append(None)
                args[key][attr] = world[i]

        result = []

        for key in default_order(components):
            klass, ar, kw = classes[key]
            result.append(klass(*args[key], *ar, **kwargs[key], **kw))

        if len(result) == 1:
            return result[0]
        else:
            return result

    def array_index_to_world(self, *index_arrays):
        return self.pixel_to_world(*index_arrays[::-1])

    def world_to_array_index(self, *world_objects):
        return tuple(np.round(self.world_to_pixel(*world_objects)[::-1]).astype(int).tolist())
