from distutils.version import LooseVersion

import matplotlib
from matplotlib import pyplot as plt

from ..utils.decorators import wraps

MPL_VERSION = LooseVersion(matplotlib.__version__)

ROOT = "http://{server}/testing/astropy/2017-07-12T14:12:26.217559/{mpl_version}/"

IMAGE_REFERENCE_DIR = (ROOT.format(server='data.astropy.org', mpl_version='1.5.x') + ',' +
                       ROOT.format(server='www.astropy.org/astropy-data', mpl_version='1.5.x'))


def ignore_matplotlibrc(func):
    # This is a decorator for tests that use matplotlib but not pytest-mpl
    # (which already handles rcParams)
    @wraps(func)
    def wrapper(*args, **kwargs):
        with plt.style.context({}, after_reset=True):
            return func(*args, **kwargs)
    return wrapper
