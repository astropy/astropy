# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import numpy as np

__all__ = ['deserialize_class', 'wcs_info_str']


def deserialize_class(tpl, construct=True):
    """
    Deserialize classes recursively.
    """

    if not isinstance(tpl, tuple) or len(tpl) != 3:
        raise ValueError("Expected a tuple of three values")

    module, klass = tpl[0].rsplit('.', 1)
    module = importlib.import_module(module)
    klass = getattr(module, klass)

    args = tuple([deserialize_class(arg) if isinstance(arg, tuple) else arg for arg in tpl[1]])

    kwargs = dict((key, deserialize_class(val)) if isinstance(val, tuple) else (key, val) for (key, val) in tpl[2].items())

    if construct:
        return klass(*args, **kwargs)
    else:
        return klass, args, kwargs


def wcs_info_str(wcs):

    # Overall header

    s = f'{wcs.__class__.__name__} Transformation\n\n'
    s += ('This transformation has {} pixel and {} world dimensions\n\n'
            .format(wcs.pixel_n_dim, wcs.world_n_dim))
    s += f'Array shape (Numpy order): {wcs.array_shape}\n\n'

    # Pixel dimensions table

    array_shape = wcs.array_shape or (0,)
    pixel_shape = wcs.pixel_shape or (None,) * wcs.pixel_n_dim

    # Find largest between header size and value length
    pixel_dim_width = max(9, len(str(wcs.pixel_n_dim)))
    pixel_nam_width = max(9, max(len(x) for x in wcs.pixel_axis_names))
    pixel_siz_width = max(9, len(str(max(array_shape))))

    s += (('{0:' + str(pixel_dim_width) + 's}').format('Pixel Dim') + '  ' +
            ('{0:' + str(pixel_nam_width) + 's}').format('Axis Name') + '  ' +
            ('{0:' + str(pixel_siz_width) + 's}').format('Data size') + '  ' +
            'Bounds\n')

    for ipix in range(wcs.pixel_n_dim):
        s += (('{0:' + str(pixel_dim_width) + 'g}').format(ipix) + '  ' +
                ('{0:' + str(pixel_nam_width) + 's}').format(wcs.pixel_axis_names[ipix] or 'None') + '  ' +
                (" " * 5 + str(None) if pixel_shape[ipix] is None else
                ('{0:' + str(pixel_siz_width) + 'g}').format(pixel_shape[ipix])) + '  ' +
                '{:s}'.format(str(None if wcs.pixel_bounds is None else wcs.pixel_bounds[ipix]) + '\n'))

    s += '\n'

    # World dimensions table

    # Find largest between header size and value length
    world_dim_width = max(9, len(str(wcs.world_n_dim)))
    world_nam_width = max(9, max(len(x) if x is not None else 0 for x in wcs.world_axis_names))
    world_typ_width = max(13, max(len(x) if x is not None else 0 for x in wcs.world_axis_physical_types))

    s += (('{0:' + str(world_dim_width) + 's}').format('World Dim') + '  ' +
            ('{0:' + str(world_nam_width) + 's}').format('Axis Name') + '  ' +
            ('{0:' + str(world_typ_width) + 's}').format('Physical Type') + '  ' +
            'Units\n')

    for iwrl in range(wcs.world_n_dim):

        name = wcs.world_axis_names[iwrl] or 'None'
        typ = wcs.world_axis_physical_types[iwrl] or 'None'
        unit = wcs.world_axis_units[iwrl] or 'unknown'

        s += (('{0:' + str(world_dim_width) + 'd}').format(iwrl) + '  ' +
                ('{0:' + str(world_nam_width) + 's}').format(name) + '  ' +
                ('{0:' + str(world_typ_width) + 's}').format(typ) + '  ' +
                '{:s}'.format(unit + '\n'))
    s += '\n'

    # Axis correlation matrix

    pixel_dim_width = max(3, len(str(wcs.world_n_dim)))

    s += 'Correlation between pixel and world axes:\n\n'

    s += (' ' * world_dim_width + '  ' +
            ('{0:^' + str(wcs.pixel_n_dim * 5 - 2) + 's}').format('Pixel Dim') +
            '\n')

    s += (('{0:' + str(world_dim_width) + 's}').format('World Dim') +
            ''.join(['  ' + ('{0:' + str(pixel_dim_width) + 'd}').format(ipix)
                    for ipix in range(wcs.pixel_n_dim)]) +
            '\n')

    matrix = wcs.axis_correlation_matrix
    matrix_str = np.empty(matrix.shape, dtype='U3')
    matrix_str[matrix] = 'yes'
    matrix_str[~matrix] = 'no'

    for iwrl in range(wcs.world_n_dim):
        s += (('{0:' + str(world_dim_width) + 'd}').format(iwrl) +
                ''.join(['  ' + ('{0:>' + str(pixel_dim_width) + 's}').format(matrix_str[iwrl, ipix])
                        for ipix in range(wcs.pixel_n_dim)]) +
                '\n')

    # Make sure we get rid of the extra whitespace at the end of some lines
    return '\n'.join([l.rstrip() for l in s.splitlines()])


def convert_between_array_and_pixel_axes(axis, naxes):
    """
    Converts between array axis indices and WCS pixel axis indices.

    The translation is symmetric so works in either direction.

    Parameters
    ----------
    axis: `numpy.ndarray` of `int`
        The array or WCS pixel axis number(s) to convert.

    naxes: `int`
        The number of array axes.
        The number of array axes should be the same as the number of WCS pixel axes.

    Returns
    -------
    new_axis: `numpy.ndarray` of `int`
        The axis number(s) after reflection.
    """
    # Check type of input.
    if not isinstance(axis, np.ndarray):
        raise TypeError("input must be of array type. Got type: {type(axis)}")
    if axis.dtype.char not in np.typecodes['AllInteger']:
        raise TypeError("input dtype must be of int type.  Got dtype: {axis.dtype})")
    # Convert negative indices to positive equivalents.
    axis[axis < 0] += naxes
    if any(axis > naxes - 1):
        raise IndexError("Axis out of range.  "
                         f"Number of axes = {naxes}; Axis numbers requested = {axes}")
    # Reflect axis about center of number of axes.
    new_axis = naxes - 1 - axis
    return new_axis


def pixel_axis_to_world_axes(pixel_axis, axis_correlation_matrix):
    """
    Retrieves the indices of the world axis index corresponding to a WCS pixel axis.

    A world axis index is equivalent to the index of the relevant physical type
    in `~astropy.wcs.WCS.world_axis_physical_types`.

    Parameters
    ----------
    pixel_axis: `int`
        The WCS pixel axis index/indices for which the world axes are desired.

    axis_correlation_matrix: `numpy.ndarray` of `bool`
        2D boolean correlation matrix defining the dependence between the pixel and world axes.
        Format same as `astropy.wcs.BaseLowLevelWCS.axis_correlation_matrix`.

    Returns
    -------
    world_axes: `numpy.ndarray`
        The world axis indices corresponding to the pixel axis.
    """
    return np.arange(axis_correlation_matrix.shape[0])[axis_correlation_matrix[:, pixel_axis]]


def world_axis_to_pixel_axes(world_axis, axis_correlation_matrix):
    """
    Gets the WCS pixel axis indices corresponding to a given WCS world axis.

    A world axis index is equivalent to the index of the relevant physical type
    in `~astropy.wcs.WCS.world_axis_physical_types`.

    Parameters
    ----------
    world_axis: `int`
        The index of the world axis for which the pixel axes are desired.

    axis_correlation_matrix: `numpy.ndarray` of `bool`
        2D boolean correlation matrix defining the dependence between the pixel and world axes.
        Format same as `astropy.wcs.BaseLowLevelWCS.axis_correlation_matrix`.

    Returns
    -------
    pixel_axes: `numpy.ndarray`
        The pixel axis indices corresponding to the world axis.
    """
    return np.arange(axis_correlation_matrix.shape[1])[axis_correlation_matrix[world_axis]]


def pixel_axis_to_physical_types(pixel_axis, wcs):
    """
    Gets the world axis physical types corresponding to a WCS pixel axis.

    Parameters
    ----------
    pixel_axis: `int`
        The pixel axis number(s) for which the world axis numbers are desired.

    wcs: `astropy.wcs.BaseLowLevelWCS`
        The WCS object defining the relationship between pixel and world axes.

    Returns
    -------
    physical_types: `numpy.ndarray` of `str`
        The physical types corresponding to the pixel axis.
    """
    return np.array(wcs.world_axis_physical_types)[wcs.axis_correlation_matrix[:, pixel_axis]]


def physical_type_to_pixel_axes(physical_type, wcs):
    """
    Gets the WCS pixel axis indices corresponding to a world axis physical type.

    Parameters
    ----------
    physical_type: `int`
        The pixel axis number(s) for which the world axis numbers are desired.

    wcs: `astropy.wcs.BaseLowLevelWCS`
        The WCS object defining the relationship between pixel and world axes.

    Returns
    -------
    pixel_axes: `numpy.ndarray`
        The pixel axis indices corresponding to the physical type.
    """
    world_axis = physical_type_to_world_axis(physical_type, wcs.world_axis_physical_types)
    return world_axis_to_pixel_axes(world_axis, wcs.axis_correlation_matrix)


def physical_type_to_world_axis(physical_type, world_axis_physical_types):
    """
    Returns world axis index of a physical type.

    The world axis index is equivalent to the index of the relevant physical type
    in `~astropy.wcs.WCS.world_axis_physical_types`.
    Input can be a substring of a physical type, but
    substrings must be unique to a single physical type.

    Parameters
    ----------
    physical_type: `str`
        The physical type or a substring unique to a single physical type.

    world_axis_physical_types: sequence of `str`
        All available physical types.  Ordering must be same as
        `astropy.wcs.BaseLowLevelWCS.world_axis_physical_types`

    Returns
    -------
    world_axis: `numbers.Integral`
        The world axis index of the physical type.
    """
    # Find world axis index described by physical type.
    widx = np.where(world_axis_physical_types == physical_type)[0]
    # If physical type does not correspond to entry in world_axis_physical_types,
    # check if it is a substring of any physical types.
    if len(widx) == 0:
        widx = [physical_type in world_axis_physical_type
                for world_axis_physical_type in world_axis_physical_types]
        widx = np.arange(len(world_axis_physical_types))[widx]
    if len(widx) != 1:
        raise ValueError(
                "Input does not uniquely correspond to a physical type."
                f" Expected unique substring of one of {world_axis_physical_types}."
                f"  Got: {physical_type}")
    # Return axes with duplicates removed.
    return widx[0]


def get_dependent_pixel_axes(pixel_axis, axis_correlation_matrix):
    """
    Find indices of all pixel axes associated with the world axes linked to the input pixel axis.

    For example, say the input pixel axis is 0 and it is associated with two world axes
    corresponding to longitude and latitude. Let's say that pixel axis 1 is also
    associated with longitude and latitude. Thus, this function would return pixel axes 0 and 1.
    On the other hand let's say pixel axis 2 is associated with only one world axis,
    e.g. wavelength, which does not depend on any other pixel axis (i.e. it is independent).
    In that case this function would only return pixel axis 2.
    Both input and output pixel axis indices are in the WCS ordering convention
    (reverse of numpy ordering convention).
    The returned axis indices include the input axis.

    Parameters
    ----------
    wcs_axis: `int`
        Index of axis (in WCS ordering convention) for which dependent axes are desired.

    axis_correlation_matrix: `numpy.ndarray` of `bool`
        2D boolean correlation matrix defining the dependence between the pixel and world axes.
        Format same as `astropy.wcs.BaseLowLevelWCS.axis_correlation_matrix`.

    Returns
    -------
    dependent_pixel_axes: `np.ndarray` of `int`
        Sorted indices of pixel axes dependent on input axis in WCS ordering convention.
    """
    # The axis_correlation_matrix is (n_world, n_pixel) but we want to know
    # which pixel coordinates are linked to which other pixel coordinates.
    # To do this we take a column from the matrix and find if there are
    # any entries in common with all other columns in the matrix.
    world_dep = axis_correlation_matrix[:, pixel_axis:pixel_axis + 1]
    dependent_pixel_axes = np.sort(np.nonzero((world_dep & axis_correlation_matrix).any(axis=0))[0])
    return dependent_pixel_axes


def get_dependent_array_axes(array_axis, axis_correlation_matrix):
    """
    Find indices of all array axes associated with the world axes linked to the input array axis.

    For example, say the input array axis is 0 and it is associated with two world axes
    corresponding to longitude and latitude. Let's say that array axis 1 is also
    associated with longitude and latitude. Thus, this function would return array axes 0 and 1.
    On the other hand let's say array axis 2 is associated with only one world axis,
    e.g. wavelength, which does not depend on any other array axis (i.e. it is independent).
    In that case this function would only return array axis 2.
    Both input and output array axis indices are in the numpy array ordering convention
    (reverse of WCS ordering convention).
    The returned axis indices include the input axis.

    Parameters
    ----------
    array_axis: `int`
        Index of array axis (in numpy ordering convention) for which dependent axes are desired.

    axis_correlation_matrix: `numpy.ndarray` of `bool`
        2D boolean correlation matrix defining the dependence between the pixel and world axes.
        Format same as `astropy.wcs.BaseLowLevelWCS.axis_correlation_matrix`.

    Returns
    -------
    dependent_array_axes: `np.ndarray` of `int`
        Sorted indices of array axes dependent on input axis in numpy ordering convention.
    """
    naxes = axis_correlation_matrix.shape[1]
    pixel_axis = convert_between_array_and_pixel_axes(np.array([array_axis], dtype=int), naxes)[0]
    dependent_pixel_axes = get_dependent_pixel_axes(pixel_axis, axis_correlation_matrix)
    dependent_array_axes = convert_between_array_and_pixel_axes(dependent_pixel_axes, naxes)
    return np.sort(dependent_array_axes)


def get_dependent_world_axes(world_axis, axis_correlation_matrix):
    """
    Given a WCS world axis index, return indices of dependent WCS world axes.

    Both input and output axis indices are in the WCS ordering convention
    (reverse of numpy ordering convention).
    The returned axis indices include the input axis.

    Parameters
    ----------
    world_axis: `int`
        Index of axis (in WCS ordering convention) for which dependent axes are desired.

    axis_correlation_matrix: `numpy.ndarray` of `bool`
        2D boolean correlation matrix defining the dependence between the pixel and world axes.
        Format same as `astropy.wcs.BaseLowLevelWCS.axis_correlation_matrix`.

    Returns
    -------
    dependent_world_axes: `np.ndarray` of `int`
        Sorted indices of pixel axes dependent on input axis in WCS ordering convention.
    """
    # The axis_correlation_matrix is (n_world, n_pixel) but we want to know
    # which world coordinates are linked to which other world coordinates.
    # To do this we take a row from the matrix and find if there are
    # any entries in common with all other rows in the matrix.
    pixel_dep = axis_correlation_matrix[world_axis:world_axis + 1].T
    dependent_world_axes = np.sort(np.nonzero((pixel_dep & axis_correlation_matrix).any(axis=1))[0])
    return dependent_world_axes


def get_dependent_physical_types(physical_type, wcs):
    """
    Given a world axis physical type, return the dependent physical types including the input type.

    Parameters
    ----------
    physical_type: `str`
        The world axis physical types whose dependent physical types are desired.

    wcs: `astropy.wcs.BaseLowLevelWCS`
        The WCS object defining the relationship between pixel and world axes.

    Returns
    -------
    dependent_physical_types: `np.ndarray` of `str`
        Physical types dependent on the input physical type.
    """
    world_axis_physical_types = wcs.world_axis_physical_types
    world_axis = physical_type_to_world_axis(physical_type, world_axis_physical_types)
    dependent_world_axes = get_dependent_world_axes(world_axis, wcs.axis_correlation_matrix)
    dependent_physical_types = np.array(world_axis_physical_types)[dependent_world_axes]
    return dependent_physical_types
