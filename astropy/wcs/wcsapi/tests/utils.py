import inspect

import numpy as np
from typeguard import check_type

from astropy.wcs.wcsapi import BaseLowLevelWCS


__all__ = ['validate_low_level_wcs_types']


def isproperty(obj):
    return isinstance(obj, property)


def validate_low_level_wcs_types(wcs):
    """
    Compare the return values from ``wcs`` with the type annotations on
    `astropy.wcs.wcsapi.BaseLowLevelWCS`.
    """

    if not isinstance(wcs, BaseLowLevelWCS):
        raise TypeError("wcs argument must be an instance of BaseLowLevelWCS")

    # Get a list of all the properties on BaseLowLevelWCS
    properties = [name for name, _ in inspect.getmembers(BaseLowLevelWCS, isproperty)]

    # Compare the return value of that property with the type annotations
    for prop_name in properties:
        check_type(prop_name,
                   getattr(wcs, prop_name),
                   getattr(BaseLowLevelWCS, prop_name).fget.__annotations__['return'])

    # Validate the return types of the pixel <> world functions

    # Set up valid but meaningless inputs, the idea here is just to validate
    # the *types* of the returns
    scalar_pixel_default_call = tuple([0]*wcs.pixel_n_dim)
    scalar_pixel_call = [b[0] for b in wcs.pixel_bounds] if wcs.pixel_bounds is not None else scalar_pixel_default_call
    array_pixel_call = np.broadcast_to(scalar_pixel_call, (10, len(scalar_pixel_call))).T
    scalar_world_call = wcs.pixel_to_world_values(*scalar_pixel_call)
    array_world_call = wcs.pixel_to_world_values(*array_pixel_call)
    if wcs.world_n_dim == 1:
        scalar_world_call = (scalar_world_call,)
        array_world_call = (array_world_call,)

    for pixel_call in (scalar_pixel_call, array_pixel_call):
        return_value = wcs.pixel_to_world_values(*pixel_call)
        check_type("pixel_to_world_values",
                   return_value,
                   BaseLowLevelWCS.pixel_to_world_values.__annotations__['return'])

        # The return type should not be a tuple if there is only one world dimension
        if wcs.world_n_dim == 1:
            assert not isinstance(return_value, tuple)

        return_value = wcs.array_index_to_world_values(*pixel_call[::-1])
        check_type("array_index_to_world_values",
                   return_value,
                   BaseLowLevelWCS.pixel_to_world_values.__annotations__['return'])

        # The return type should not be a tuple if there is only one world dimension
        if wcs.world_n_dim == 1:
            assert not isinstance(return_value, tuple)

    for world_call in (scalar_world_call, array_world_call):
        return_value = wcs.world_to_pixel_values(*world_call)
        check_type("world_to_pixel_values",
                   return_value,
                   BaseLowLevelWCS.world_to_pixel_values.__annotations__['return'])

        # The return type should not be a tuple if there is only one pixel dimension
        if wcs.pixel_n_dim == 1:
            assert not isinstance(return_value, tuple)

        return_value = wcs.world_to_array_index_values(*world_call[::-1])
        check_type("world_to_array_index_values",
                   return_value,
                   BaseLowLevelWCS.world_to_array_index_values.__annotations__['return'])

        # The return type should not be a tuple if there is only one pixel dimension
        if wcs.pixel_n_dim == 1:
            assert not isinstance(return_value, tuple)
