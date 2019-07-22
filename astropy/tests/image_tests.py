import matplotlib
from matplotlib import pyplot as plt

from astropy.utils.decorators import wraps

MPL_VERSION = matplotlib.__version__

# The developer versions of the form 3.1.x+... contain changes that will only
# be included in the 3.2.x release, so we update this here.
if MPL_VERSION[:3] == '3.1' and '+' in MPL_VERSION:
    MPL_VERSION = '3.2'

ROOT = "http://{server}/testing/astropy/2018-10-24T12:38:34.134556/{mpl_version}/"
IMAGE_REFERENCE_DIR = (ROOT.format(server='data.astropy.org', mpl_version=MPL_VERSION[:3] + '.x') + ',' +
                       ROOT.format(server='www.astropy.org/astropy-data', mpl_version=MPL_VERSION[:3] + '.x'))


def ignore_matplotlibrc(func):
    # This is a decorator for tests that use matplotlib but not pytest-mpl
    # (which already handles rcParams)
    @wraps(func)
    def wrapper(*args, **kwargs):
        with plt.style.context({}, after_reset=True):
            return func(*args, **kwargs)
    return wrapper
