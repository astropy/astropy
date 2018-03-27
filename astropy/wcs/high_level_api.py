import importlib
from collections import defaultdict, OrderedDict

__all__ = ['HighLevelWCS']


def rec_getattr(obj, att):
    for a in att.split('.'):
        obj = getattr(obj, a)
    return obj


def deserialize_class(tpl, construct=True):

    if not isinstance(tpl, tuple) or len(tpl) != 3:
        raise ValueError("Expected a tuple of three values")

    module, klass = tpl[0].rsplit('.', 1)
    module = importlib.import_module(module)
    klass = getattr(module, klass)

    args = [deserialize_class(arg) if isinstance(arg, tuple) else arg for arg in tpl[1]]

    kwargs = dict((key, deserialize_class(val)) if isinstance(val, tuple) else (key, val) for (key ,val) in tpl[2].items())

    if construct:
        return klass(*args, **kwargs)
    else:
        return klass, args, kwargs


def default_order(components):
    order = []
    for key, _, _ in components:
        if key not in order:
            order.append(key)
    return order


class HighLevelWCS(object):

    def __init__(self, wcs):
        self._wcs = wcs

    def world_to_pixel(self, *world):

        # Cache the classes and components since this may be expensive
        serialized_classes = self._wcs.world_axis_object_classes
        components = self._wcs.world_axis_object_components

        # Deserialize world_axis_object_classes using the default order
        classes = OrderedDict()
        for key in default_order(components):
            classes[key] = deserialize_class(serialized_classes[key],
                                             construct=False)

        # Check that the number of classes matches the number of inputs
        if len(world) != len(classes):
            raise ValueError("Number of world inputs ({0}) does not match "
                             "expected ({1})".format(len(world), len(classes)))

        # Determine whether the classes are uniquely matched, that is we check
        # whether there is only one of each class.
        world_by_key = {}
        unique_match = True
        for w in world:
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
                objects[key] = klass(world_by_key[key], *args, **kwargs)

        else:

            for ikey, key in enumerate(classes):
                klass, args, kwargs = classes[key]
                w = world[ikey]
                if not isinstance(w, klass):
                    raise ValueError("Expected the following order of world "
                                     "arguments: {0}".format(', '.join([k.__name__ for (k, _, _) in classes.values()])))
                objects[key] = klass(w, *args, **kwargs)

        # We now extract the attributes needed for the world values
        world = []
        for key, _, attr in components:
            world.append(rec_getattr(objects[key], attr))

        # Finally we convert to pixel coordinates
        pixel = self._wcs.world_to_pixel_values(*world)

        return pixel

    def pixel_to_world(self, *pixel):

        # Compute the world coordinate values
        world = self._wcs.pixel_to_world_values(*pixel)

        # Cache the classes and components since this may be expensive
        components = self._wcs.world_axis_object_components
        classes = self._wcs.world_axis_object_classes

        # Deserialize classes
        classes_new = {}
        for key, value in classes.items():
            classes_new[key] = deserialize_class(value, construct=False)

        args = defaultdict(list)
        kwargs = defaultdict(dict)

        for i, (key, attr, _) in enumerate(components):
            if isinstance(attr, str):
                kwargs[attr] = world[i]
            else:
                while attr > len(args[key]) - 1:
                    args[key].append(None)
                args[key][attr] = world[i]

        result = []

        for key in default_order(components):
            klass, ar, kw = classes_new[key]
            result.append(klass(*args[key], *ar, **kwargs[key], **kw))

        return result
