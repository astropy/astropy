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
