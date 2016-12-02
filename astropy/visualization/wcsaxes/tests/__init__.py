import os
from distutils.version import LooseVersion

import matplotlib

MPL_VERSION = LooseVersion(matplotlib.__version__)

ROOT = "https://astropy.stsci.edu/data/wcsaxes/2016-05-04T10:26:13.545916"

if MPL_VERSION >= LooseVersion('1.5.0'):
    baseline_dir = ROOT + '/1.5.x/'
else:
    baseline_dir = ROOT + '/1.4.x/'
