"""
A simple class to manage a piece of global science state.  See
:ref:`config-developer` for more details.
"""
from contextlib import contextmanager
import warnings

from ..config import ConfigItem
from .exceptions import AstropyDeprecationWarning
from ..utils import find_current_module


__all__ = ['ScienceState', 'ScienceStateAlias']


class ScienceState(object):
    """
    Science state subclasses are used to manage global items that can
    affect science results.  Subclasses will generally override
    `validate` to convert from any of the acceptable inputs (such as
    strings) to the appropriate internal objects, and set an initial
    value to the ``_value`` member so it has a default.

    Examples
    --------

    ::

        class MyState(ScienceState):
            @classmethod
            def validate(cls, value):
                if value not in ('A', 'B', 'C'):
                    raise ValueError("Must be one of A, B, C")
                return value
    """

    def __init__(self):
        raise RuntimeError(
            "This class is a singleton.  Do not instantiate.")

    @classmethod
    def get(cls):
        """
        Get the current science state value.
        """
        return cls.validate(cls._value)

    @classmethod
    def set(cls, value):
        """
        Set the current science state value.
        """
        class _Context(object):
            def __init__(self, parent, value):
                self._value = value
                self._parent = parent

            def __enter__(self):
                pass

            def __exit__(self, type, value, tb):
                self._parent._value = self._value

        ctx = _Context(cls, cls._value)
        value = cls.validate(value)
        cls._value = value
        return ctx

    @classmethod
    def validate(cls, value):
        """
        Validate the value and convert it to its native type, if
        necessary.
        """
        return value


class ScienceStateAlias(ConfigItem):
    """
    This is a backward compatibility layer for configuration items
    that moved to ScienceState classes in astropy 0.4.

    Parameters
    ----------
    since : str
        The version in which the configuration item was converted into
        science state.

    python_name : str
        The old name of the Python variable for the configuration item.

    config_name : str
        The old name of the configuration item in the configuration file.

    science_state : ScienceState subclass
        The science state class that now manages this information.

    cfgtype : str or `None`, optional
        A type specifier like those used as the *values* of a particular key
        in a ``configspec`` file of ``configobj``. If `None`, the type will
        be inferred from the default value.

    module : str, optional
        The module containing the old configuration item.
    """
    # REMOVE in astropy 0.5

    def __init__(self, since, python_name, config_name, science_state,
                 cfgtype=None, module=None):
        # We have to do the automatic module determination here, not
        # just in ConfigItem, otherwise the extra stack frame will
        # make it come up with the wrong answer.
        if module is None:
            module = find_current_module(2)
            if module is None:
                msg1 = 'Cannot automatically determine get_config module, '
                msg2 = 'because it is not called from inside a valid module'
                raise RuntimeError(msg1 + msg2)
            else:
                module = module.__name__

        self._dont_warn = True
        self._since = since
        self._python_name = python_name
        self._config_name = config_name
        self._science_state = science_state
        super(ScienceStateAlias, self).__init__(
            science_state._value,
            science_state.__doc__,
            cfgtype=cfgtype,
            module=module)
        self.name = config_name

        # Set the default value of the science state from the config
        # file, if defined.  This is what pulls any old values in the
        # config file and applies them to the science state.
        value = super(ScienceStateAlias, self).__call__()

        # We got a value in the config file
        if science_state.validate(value) != science_state._value:
            warnings.warn(
                "Config parameter '{0}' in section [{1}] is deprecated. "
                "Use science state {2}.{3} instead.".format(
                    self.name, self.module, self._science_state.__module__,
                    self._science_state.__name__),
                AstropyDeprecationWarning)

        del self._dont_warn

    def _deprecation_warning(self, extra=None):
        if hasattr(self, '_dont_warn'):
            return

        message = ("'{0}.{1}' is deprecated, and is no longer defined "
                   "as a configuration item.".format(
                       self.module, self._python_name))
        if extra is not None:
            message += " Use '{0}.{1}{2}' instead.".format(
                self._science_state.__module__,
                self._science_state.__name__,
                extra)
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
        self._dont_warn = True
        try:
            result = super(ScienceStateAlias, self).reload()
        finally:
            del self._dont_warn
        if result is not None:
            self._science_state.set(result)

    def __call__(self):
        self._deprecation_warning('.get()')
        return self._science_state.get()
