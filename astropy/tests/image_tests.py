from distutils.version import LooseVersion

import matplotlib

MPL_VERSION = LooseVersion(matplotlib.__version__)

ROOT = "http://{server}/testing/astropy/2017-07-12T14:12:26.217559/{mpl_version}/"

IMAGE_REFERENCE_DIR = (ROOT.format(server='data.astropy.org', mpl_version='1.5.x') + ',' +
                       ROOT.format(server='www.astropy.org/astropy-data', mpl_version='1.5.x'))
