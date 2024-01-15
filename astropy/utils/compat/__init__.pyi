# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .misc import (
    PYTHON_LT_3_10 as PYTHON_LT_3_10,
    PYTHON_LT_3_11 as PYTHON_LT_3_11,
    override__dir__ as override__dir__,
)
from .numpycompat import (
    NUMPY_LT_1_22_1 as NUMPY_LT_1_22_1,
    NUMPY_LT_1_23 as NUMPY_LT_1_23,
    NUMPY_LT_1_24 as NUMPY_LT_1_24,
    NUMPY_LT_1_25 as NUMPY_LT_1_25,
    NUMPY_LT_1_26 as NUMPY_LT_1_26,
    NUMPY_LT_2_0 as NUMPY_LT_2_0,
)
from . import (
    misc as misc,
    numpycompat as numpycompat,
    optional_deps as optional_deps,
)
