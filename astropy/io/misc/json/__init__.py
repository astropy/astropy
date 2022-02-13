# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Extended Encoder and Decoder for serializing objects via the JSON protocol.

First we import the necessary libraries: :mod:`json` and custom en/decoders.

>>> import json
>>> from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder

Without using the custom en/decoders even some python builtins, like
`complex` numbers cannot be serialized with JSON.

>>> v = 1 + 2j
>>> try:
...     json.dumps(v)
... except TypeError:
...     print(":(")
:(

Using Astropy's custom en/decoders the python builtins and many :mod:`numpy`
and Astropy objects can be dumped to and from JSON.

Returning to the complex number example:

>>> serialized = json.dumps(v, cls=JSONExtendedEncoder)
>>> serialized
'{"!": "builtins.complex", "value": [1.0, 2.0]}'

>>> json.loads(serialized, cls=JSONExtendedDecoder)
(1+2j)

`~numpy.ndarray` are saved in a similar manner, with the class specified
and the value converted to a JSON-compatible format. However, :mod:`numpy`
objects are generally more complicated, with `~numpy.dtype` and other
properties. The en/decoders recursively iterate through known properties,
serializing to / loading from JSON.

>>> import numpy as np
>>> arr = np.array([3, 4], dtype=float)
>>> serialized = json.dumps(arr, cls=JSONExtendedEncoder)
>>> serialized
'{"!": "numpy.ndarray", "value": ["3.0", "4.0"], "dtype": "float64"}'

>>> json.loads(serialized, cls=JSONExtendedDecoder)
array([3., 4.])

Precision is kept; for example a `numpy.int64` number will remain an integer.

>>> arr = np.int64(10.001)
>>> serialized = json.dumps(arr, cls=JSONExtendedEncoder)
>>> serialized
'{"!": "numpy.int64", "value": "10"}'

>>> json.loads(serialized, cls=JSONExtendedDecoder)
10

Astropy objects can also be serialized to / loaded from JSON.
We start with `astropy.units.Quantity`, which builds upon the
:mod:`numpy` extensions, adding a ``"units"`` key.

>>> import astropy.units as u
>>> q = np.array([3], dtype=float) * u.km
>>> serialized = json.dumps(q, cls=JSONExtendedEncoder)
>>> serialized
'{"!": "astropy.units.Quantity",
  "value": {"!": "numpy.ndarray", "value": ["3.0"], "dtype": "float64"},
  "unit": "km"}'

>>> json.loads(serialized, cls=JSONExtendedDecoder)
<Quantity [3.] km>

`astropy.cosmology.Cosmology` can also be saved. These are complex objects
composed of many |Quantity| and arbitrary metadata.

>>> from astropy.cosmology import Planck18, units as cu
>>> serialized = json.dumps(Planck18, cls=JSONExtendedEncoder)
>>> serialized
'{"!": "astropy.cosmology.FlatLambdaCDM",
  "value": {"name": "Planck18",
            "H0": {"!": "astropy.units.Quantity",
                   "value": 67.66, "unit": "km / (Mpc s)"},
            ...
  "meta": {"Oc0": 0.2607, "n": 0.9665, ...

>>> with u.add_enabled_units(cu):
...     json.loads(serialized, cls=JSONExtendedDecoder)
FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
              Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.   0.   0.06] eV, Ob0=0.04897)
"""


from astropy.io.utils import load_all_entry_points

from . import builtins as _  # imported before entry points are loaded
from . import numpy as _
from .core import *
from .core import _json_base_encode  # imported before entry points are loaded

# After importing, load all entry points
load_all_entry_points('astropy_io_json_extensions')

del load_all_entry_points
