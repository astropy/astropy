# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import numpy as np


__all__ = ['deserialize_class', 'wcs_info_str']

import numpy as np

from astropy.utils.misc import unbroadcast

__all__ = ['deserialize_class', 'efficient_pixel_to_pixel']


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

    print(pixel_nam_width)

    s += (('{0:' + str(pixel_dim_width) + 's}').format('Pixel Dim') + '  ' +
            ('{0:' + str(pixel_nam_width) + 's}').format('Axis Name') + '  ' +
            ('{0:' + str(pixel_siz_width) + 's}').format('Data size') + '  ' +
            'Bounds\n')

    for ipix in range(wcs.pixel_n_dim):
        s += (('{0:' + str(pixel_dim_width) + 'd}').format(ipix) + '  ' +
                ('{0:' + str(pixel_nam_width) + 's}').format(wcs.pixel_axis_names[ipix] or 'None') + '  ' +
                (" " * 5 + str(None) if pixel_shape[ipix] is None else
                ('{0:' + str(pixel_siz_width) + 'd}').format(pixel_shape[ipix])) + '  ' +
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

def unique_with_order_preserved(items):
    """
    Return a list of unique items in the list provided, preserving the order
    in which they are found.
    """
    new_items = []
    for item in items:
        if item not in new_items:
            new_items.append(item)
    return new_items


def pixel_to_world_correlation_matrix(wcs):
    """
    Return a correlation matrix between the pixel coordinates and the
    high level world coordinates, along with the list of high level world
    coordinate classes.

    The shape of the matrix is ``(n_world, n_pix)``, where ``n_world`` is the
    number of high level world coordinates.
    """

    # We basically want to collapse the world dimensions together that are
    # combined into the same high-level objects.

    # Get the following in advance as getting these properties can be expensive
    all_components = wcs.low_level_wcs.world_axis_object_components
    all_classes = wcs.low_level_wcs.world_axis_object_classes
    axis_correlation_matrix = wcs.low_level_wcs.axis_correlation_matrix

    components = unique_with_order_preserved([c[0] for c in all_components])

    matrix = np.zeros((len(components), wcs.pixel_n_dim), dtype=bool)

    for iworld in range(wcs.world_n_dim):
        iworld_unique = components.index(all_components[iworld][0])
        matrix[iworld_unique] |= axis_correlation_matrix[iworld]

    classes = [all_classes[component][0] for component in components]

    return matrix, classes


def pixel_to_pixel_correlation_matrix(wcs_in, wcs_out):
    """
    Correlation matrix between the input and output pixel coordinates for a
    pixel -> world -> pixel transformation specified by two WCS instances.

    The first WCS specified is the one used for the pixel -> world
    transformation and the second WCS specified is the one used for the world ->
    pixel transformation. The shape of the matrix is
    ``(n_pixel_out, n_pixel_in)``.
    """

    matrix1, classes1 = pixel_to_world_correlation_matrix(wcs_in)
    matrix2, classes2 = pixel_to_world_correlation_matrix(wcs_out)

    if len(classes1) != len(classes2):
        raise ValueError("The two WCS return a different number of world coordinates")

    # Check if classes match uniquely
    unique_match = True
    mapping = []
    for class1 in classes1:
        matches = classes2.count(class1)
        if matches == 0:
            raise ValueError("The world coordinate types of the two WCS do not match")
        elif matches > 1:
            unique_match = False
            break
        else:
            mapping.append(classes2.index(class1))

    if unique_match:

        # Classes are unique, so we need to re-order matrix2 along the world
        # axis using the mapping we found above.
        matrix2 = matrix2[mapping]

    elif classes1 != classes2:

        raise ValueError("World coordinate order doesn't match and automatic matching is ambiguous")

    matrix = np.matmul(matrix2.T, matrix1)

    return matrix


def split_matrix(matrix):
    """
    Given an axis correlation matrix from a WCS object, return information about
    the individual WCS that can be split out.

    The output is a list of tuples, where each tuple contains a list of
    pixel dimensions and a list of world dimensions that can be extracted to
    form a new WCS. For example, in the case of a spectral cube with the first
    two world coordinates being the celestial coordinates and the third
    coordinate being an uncorrelated spectral axis, the matrix would look like::

        array([[ True,  True, False],
               [ True,  True, False],
               [False, False,  True]])

    and this function will return ``[([0, 1], [0, 1]), ([2], [2])]``.
    """

    pixel_used = []

    split_info = []

    for ipix in range(matrix.shape[1]):
        if ipix in pixel_used:
            continue
        pixel_include = np.zeros(matrix.shape[1], dtype=bool)
        pixel_include[ipix] = True
        n_pix_prev, n_pix = 0, 1
        while n_pix > n_pix_prev:
            world_include = matrix[:, pixel_include].any(axis=1)
            pixel_include = matrix[world_include, :].any(axis=0)
            n_pix_prev, n_pix = n_pix, np.sum(pixel_include)
        pixel_indices = list(np.nonzero(pixel_include)[0])
        world_indices = list(np.nonzero(world_include)[0])
        pixel_used.extend(pixel_indices)
        split_info.append((pixel_indices, world_indices))

    return split_info


def efficient_pixel_to_pixel(wcs_in, wcs_out, *inputs):
    """
    Transform pixel coordinates in a dataset with a WCS to pixel coordinates
    in another dataset with a different WCS. This function is designed to
    efficiently deal with input pixel arrays that are broadcasted views of
    smaller arrays.

    Parameters
    ----------
    wcs_in : `~astropy.wcs.wcsapi.BaseHighLevelWCS`
        A WCS object for the original dataset which complies with the
        high-level shared APE 14 WCS API.
    wcs_out : `~astropy.wcs.wcsapi.BaseHighLevelWCS`
        A WCS object for the target dataset which complies with the
        high-level shared APE 14 WCS API.
    *inputs :
        Scalars or arrays giving the pixel coordinates to transform.
    """

    # Shortcut for scalars
    if np.isscalar(inputs[0]):
        world_outputs = wcs_in.pixel_to_world(*inputs)
        if not isinstance(world_outputs, (tuple, list)):
            world_outputs = (world_outputs,)
        return wcs_out.world_to_pixel(*world_outputs)

    # Remember original shape
    original_shape = inputs[0].shape

    matrix = pixel_to_pixel_correlation_matrix(wcs_in, wcs_out)
    split_info = split_matrix(matrix)

    outputs = [None] * wcs_out.pixel_n_dim

    for (pixel_in_indices, pixel_out_indices) in split_info:

        pixel_inputs = []
        for ipix in range(wcs_in.pixel_n_dim):
            if ipix in pixel_in_indices:
                pixel_inputs.append(unbroadcast(inputs[ipix]))
            else:
                pixel_inputs.append(inputs[ipix].flat[0])

        pixel_inputs = np.broadcast_arrays(*pixel_inputs)

        world_outputs = wcs_in.pixel_to_world(*pixel_inputs)
        if not isinstance(world_outputs, (tuple, list)):
            world_outputs = (world_outputs,)
        pixel_outputs = wcs_out.world_to_pixel(*world_outputs)

        for ipix in range(wcs_out.pixel_n_dim):
            if ipix in pixel_out_indices:
                outputs[ipix] = np.broadcast_to(pixel_outputs[ipix], original_shape)

    return outputs
