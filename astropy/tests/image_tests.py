from distutils.version import LooseVersion

import matplotlib

MPL_VERSION = LooseVersion(matplotlib.__version__)

ROOT = "http://data.astropy.org/testing/astropy/2016-05-04T10:26:13.545916"

if MPL_VERSION >= LooseVersion('1.5.0'):
    IMAGE_REFERENCE_DIR = ROOT + '/1.5.x/'
else:
    IMAGE_REFERENCE_DIR = ROOT + '/1.4.x/'
