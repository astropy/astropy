"""
A simple class to manage a piece of global science state.
"""
from contextlib import contextmanager
import copy
import warnings

from ..config import ConfigItem
from .exceptions import AstropyDeprecationWarning


class ScienceState(object):
    def __init__(self):
        raise RuntimeError(
            "This class is a singleton.  Do not instantiate.")

    @classmethod
    def get(cls):
        return cls.validate(cls._value)

    @classmethod
    def set(cls, value):
        class _Context(object):
            def __init__(self, parent, value):
                self._value = value
                self._parent = parent

            def __enter__(self):
                pass

            def __exit__(self, type, value, tb):
                self._parent._value = self._value

        ctx = _Context(cls, copy.copy(cls._value))
        value = cls.validate(value)
        cls._value = value
        return ctx

    @classmethod
    def validate(cls, value):
        return value


class ScienceStateAlias(ConfigItem):
    """
    This is a backward compatibility layer for configuration items
    that moved to ScienceState classes in astropy 0.4.
    """
    # REMOVE in astropy 0.5

    def __init__(self, old_name, science_state):
        self._is_initing = True
        self._old_name = old_name
        self._science_state = science_state
        super(ScienceStateAlias, self).__init__(
            science_state._value,
            science_state.__doc__,
            module=science_state.__module__)
        self.name = science_state.__name__
        # Set the default value of the science state from the config
        # file, if defined, otherwise this should just effectively be
        # a no-op
        science_state._value = super(ScienceStateAlias, self).__call__()
        del self._is_initing

    def _deprecation_warning(self, extra=None):
        if hasattr(self, '_is_initing'):
            return

        message = ("'{0}.{1}' is deprecated, and is no longer defined "
                   "as a configuration item.".format(self.module, self._old_name))
        if extra is not None:
            message += " Use '{0}.{1}{2}' instead.".format(
                self.module, self.name, extra)
        warnings.warn(message, AstropyDeprecationWarning)

    def set(self, value):
        self._deprecation_warning('.set()')
        self._science_state.set(value)

    @contextmanager
    def set_temp(self, value):
        self._deprecation_warning('.set_temp()')
        with self._science_state.set(value):
            yield

    def reload(self):
        self._deprecation_warning()
        super(ScienceStateAlias, self).reload()
        result = super(ScienceStateAlias, self).reload()
        if result is not None:
            self._science_state.set(result)

    def __call__(self):
        self._deprecation_warning('.get()')
        return self._science_state.get()
