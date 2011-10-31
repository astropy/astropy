# When setup.py code loads this while trying to get setup_package.py,
# these in general wouldn't be built. But when user imports this after
# installation, these will be available.
try:
    from . import skyc
    from . import astrom
except ImportError:
    pass
