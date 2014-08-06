# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Image processing utilities for Astropy.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .scale_image import *
    from .array_utils import *
    from .sampling import *

__all__ = ['find_cutlevels', 'normalize_image', 'scale_image',
           'sigmaclip_stats', 'downsample', 'upsample', 'extract_array_2d',
           'add_array_2d', 'subpixel_indices', 'mask_to_mirrored_num']
