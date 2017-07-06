from distutils.version import LooseVersion

import matplotlib

MPL_VERSION = LooseVersion(matplotlib.__version__)

ROOT = "http://data.astropy.org/testing/astropy/2017-07-06T17:59:08.793939"

if MPL_VERSION >= LooseVersion('1.5.0'):
    IMAGE_REFERENCE_DIR = ROOT + '/1.5.x/'
else:
    IMAGE_REFERENCE_DIR = ROOT + '/1.4.x/'
