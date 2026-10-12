import configobj as _configobj

__all__ = [name for name in dir(_configobj) if not name.startswith('_')]

# We need to also add the following __getattr__ to support importing objects
# not in __all__ (such as Section)

def __getattr__(name):
    return getattr(_configobj, name)
