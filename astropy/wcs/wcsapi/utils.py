# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib

import numpy as np

__all__ = ["deserialize_class", "wcs_info_str"]


def deserialize_class(tpl, construct=True):
    """
    Deserialize classes recursively.
    """
    if not isinstance(tpl, tuple) or len(tpl) != 3:
        raise ValueError("Expected a tuple of three values")

    module, klass = tpl[0].rsplit(".", 1)
    module = importlib.import_module(module)
    klass = getattr(module, klass)

    args = tuple(
        deserialize_class(arg) if isinstance(arg, tuple) else arg for arg in tpl[1]
    )

    kwargs = dict(
        (key, deserialize_class(val)) if isinstance(val, tuple) else (key, val)
        for (key, val) in tpl[2].items()
    )

    if construct:
        return klass(*args, **kwargs)
    else:
        return klass, args, kwargs


def wcs_info_str(wcs):
    # Overall header

    if wcs.array_shape is None:
        array_shape = None
    else:
        array_shape = tuple(int(n) for n in wcs.array_shape)

    s = (
        f"{type(wcs).__name__} Transformation\n\n"
        f"This transformation has {wcs.pixel_n_dim} pixel and {wcs.world_n_dim} "
        "world dimensions\n\n"
        f"Array shape (Numpy order): {array_shape}\n\n"
    )

    # Pixel dimensions table

    array_shape = array_shape or (0,)
    pixel_shape = wcs.pixel_shape or (None,) * wcs.pixel_n_dim

    # Find largest between header size and value length
    pixel_dim_width = max(9, len(str(wcs.pixel_n_dim)))
    pixel_nam_width = max(9, *map(len, wcs.pixel_axis_names))
    pixel_siz_width = max(9, len(str(max(array_shape))))

    # fmt: off
    s += (('{0:' + str(pixel_dim_width) + 's}').format('Pixel Dim') + '  ' +
            ('{0:' + str(pixel_nam_width) + 's}').format('Axis Name') + '  ' +
            ('{0:' + str(pixel_siz_width) + 's}').format('Data size') + '  ' +
            'Bounds\n')
    # fmt: on

    if wcs.pixel_bounds is None:
        pixel_bounds = [None for _ in range(wcs.pixel_n_dim)]
    else:
        # converting to scalar arrays and back to Python with np.array(val).item()
        # guarantees that we end up with Python scalars (int or float) with
        # simple reprs, while not making any unnecessary type promotion
        # (e.g. int to float)
        pixel_bounds = [
            tuple(np.array(b).item() for b in bounds) for bounds in wcs.pixel_bounds
        ]

    for ipix in range(wcs.pixel_n_dim):
        # fmt: off
        s += (('{0:' + str(pixel_dim_width) + 'g}').format(ipix) + '  ' +
                ('{0:' + str(pixel_nam_width) + 's}').format(wcs.pixel_axis_names[ipix] or 'None') + '  ' +
                (" " * 5 + str(None) if pixel_shape[ipix] is None else
                ('{0:' + str(pixel_siz_width) + 'g}').format(pixel_shape[ipix])) + '  ' +
                f"{pixel_bounds[ipix]}\n"
              )
        # fmt: on

    s += "\n"

    # World dimensions table

    # Find largest between header size and value length
    world_dim_width = max(9, len(str(wcs.world_n_dim)))
    world_nam_width = max(9, *(len(x) for x in wcs.world_axis_names if x is not None))
    world_typ_width = max(
        [13] + [len(x) for x in wcs.world_axis_physical_types if x is not None]
    )

    # fmt: off
    s += (('{0:' + str(world_dim_width) + 's}').format('World Dim') + '  ' +
            ('{0:' + str(world_nam_width) + 's}').format('Axis Name') + '  ' +
            ('{0:' + str(world_typ_width) + 's}').format('Physical Type') + '  ' +
            'Units\n')
    # fmt: on

    for iwrl in range(wcs.world_n_dim):
        name = wcs.world_axis_names[iwrl] or "None"
        typ = wcs.world_axis_physical_types[iwrl] or "None"
        unit = wcs.world_axis_units[iwrl] or "unknown"

        # fmt: off
        s += (('{0:' + str(world_dim_width) + 'd}').format(iwrl) + '  ' +
                ('{0:' + str(world_nam_width) + 's}').format(name) + '  ' +
                ('{0:' + str(world_typ_width) + 's}').format(typ) + '  ' +
                '{:s}'.format(unit + '\n'))
        # fmt: on

    s += "\n"

    # Axis correlation matrices

    s += "Dependence of world axes on pixel axes (pixel to world):\n\n"
    s += _matrix_str(wcs.axis_correlation_matrix, "World Dim", "Pixel Dim")

    s += "\nDependence of pixel axes on world axes (world to pixel):\n\n"
    s += _matrix_str(wcs.reverse_axis_correlation_matrix, "Pixel Dim", "World Dim")

    # Make sure we get rid of the extra whitespace at the end of some lines
    return "\n".join([l.rstrip() for l in s.splitlines()])


def _matrix_str(matrix, row_label, col_label):
    """
    Format a boolean correlation matrix as a table with yes/no entries.
    """
    n_row, n_col = matrix.shape
    row_width = max(len(row_label), len(str(n_row)))
    col_width = max(3, len(str(n_col)))  # wide enough for "yes"
    lines = [
        " " * row_width + "  " + f"{col_label:^{n_col * (col_width + 2) - 2}s}",
        f"{row_label:{row_width}s}"
        + "".join(f"  {i:{col_width}d}" for i in range(n_col)),
    ]
    for irow, row in enumerate(np.where(matrix, "yes", "no")):
        lines.append(
            f"{irow:{row_width}d}" + "".join(f"  {v:>{col_width}s}" for v in row)
        )
    return "\n".join(lines) + "\n"


def _split_matrix(matrix):
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
