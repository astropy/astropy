"""
A simple class to manage a piece of global science state.  See
:ref:`astropy:config-developer` for more details.
"""


__all__ = ['ScienceState']


class _ScienceStateContext:
    def __init__(self, parent, value):
        self._value = value
        self._parent = parent

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        self._parent._value = self._value

    def __repr__(self):
        # Ensure we have a single-line repr, just in case our
        # value is not something simple like a string.
        value_repr, lb, _ = repr(self._parent._value).partition("\n")
        if lb:
            value_repr += "..."
        return f"<ScienceState {self._parent.__name__}: {value_repr}>"


class ScienceState:
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
        raise RuntimeError("This class is a singleton.  Do not instantiate.")

    @classmethod
    def get(cls):
        """
        Get the current science state value.
        """
        return cls.validate(cls._value)

    @classmethod
    def set(cls, value):
        """Set the current science state value."""
        # Create context with current value
        ctx = _ScienceStateContext(cls, cls._value)

        # Set new value
        value = cls.validate(value)
        cls._value = value

        # Return context manager
        return ctx

    @classmethod
    def validate(cls, value):
        """
        Validate the value and convert it to its native type, if
        necessary.
        """
        return value
