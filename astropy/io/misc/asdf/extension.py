# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from asdf.extension import AsdfExtension, BuiltinExtension
from asdf.util import filepath_to_url

# Make sure that all tag implementations are imported by the time we create
# the extension class so that _astropy_asdf_types is populated correctly. We
# could do this using __init__ files, except it causes pytest import errors in
# the case that asdf is not installed.
from .tags.coordinates.angle import *  # noqa: F401, F403
from .tags.coordinates.earthlocation import *  # noqa: F401, F403
from .tags.coordinates.frames import *  # noqa: F401, F403
from .tags.coordinates.representation import *  # noqa: F401, F403
from .tags.coordinates.skycoord import *  # noqa: F401, F403
from .tags.coordinates.spectralcoord import *  # noqa: F401, F403
from .tags.fits.fits import *  # noqa: F401, F403
from .tags.table.table import *  # noqa: F401, F403
from .tags.time.time import *  # noqa: F401, F403
from .tags.time.timedelta import *  # noqa: F401, F403
from .tags.transform.basic import *  # noqa: F401, F403
from .tags.transform.compound import *  # noqa: F401, F403
from .tags.transform.functional_models import *  # noqa: F401, F403
from .tags.transform.math import *  # noqa: F401, F403
from .tags.transform.physical_models import *  # noqa: F401, F403
from .tags.transform.polynomial import *  # noqa: F401, F403
from .tags.transform.powerlaws import *  # noqa: F401, F403
from .tags.transform.projections import *  # noqa: F401, F403
from .tags.transform.spline import *  # noqa: F401, F403
from .tags.transform.tabular import *  # noqa: F401, F403
from .tags.unit.equivalency import *  # noqa: F401, F403
from .tags.unit.quantity import *  # noqa: F401, F403
from .tags.unit.unit import *  # noqa: F401, F403
from .types import _astropy_asdf_types, _astropy_types

__all__ = ["AstropyExtension", "AstropyAsdfExtension"]


ASTROPY_SCHEMA_URI_BASE = "http://astropy.org/schemas/"
SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "data", "schemas")
)
ASTROPY_URL_MAPPING = [
    (
        ASTROPY_SCHEMA_URI_BASE,
        filepath_to_url(os.path.join(SCHEMA_PATH, "astropy.org"))
        + "/{url_suffix}.yaml",
    )
]


# This extension is used to register custom types that have both tags and
# schemas defined by Astropy.
class AstropyExtension(AsdfExtension):
    @property
    def types(self):
        return _astropy_types

    @property
    def tag_mapping(self):
        return [
            ("tag:astropy.org:astropy", ASTROPY_SCHEMA_URI_BASE + "astropy{tag_suffix}")
        ]

    @property
    def url_mapping(self):
        return ASTROPY_URL_MAPPING


# This extension is used to register custom tag types that have schemas defined
# by ASDF, but have tag implementations defined in astropy.
class AstropyAsdfExtension(BuiltinExtension):
    @property
    def types(self):
        return _astropy_asdf_types
