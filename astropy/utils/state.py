"""
A simple class to manage a piece of global science state.  See
:ref:`config-developer` for more details.
"""
import copy
from collections.abc import MappingView
from types import MappingProxyType

__all__ = ['ScienceState']


class _StateProxy(MappingView):
    """
    `~collections.abc.MappingView` with mix-in read-only methods from
    `~types.MappingProxyType`. Adds ``getitem``, ``keys``, ``values``,
    and ``items`` methods. Deepcopying also now returns decoupled dict.

    Parameters
    ----------
    mapping : Mapping
        The mapping on which to apply the `~collections.abc.MappingView`

    """

    def __init__(self, mapping):
        super().__init__(mapping)
        self._mappingproxy = MappingProxyType(self._mapping)  # read-only

    def __getitem__(self, key):
        """Read-only ``getitem``."""
        return self._mappingproxy[key]

    def __deepcopy__(self, memo):
        return copy.deepcopy(self._mapping, memo=memo)

    def keys(self):
        return self._mappingproxy.keys()

    def values(self):
        return self._mappingproxy.values()

    def items(self):
        return self._mappingproxy.items()


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
        class _Context:
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
                value_repr, lb, _ = repr(self._parent._value).partition('\n')
                if lb:
                    value_repr += '...'
                return ('<ScienceState {}: {}>'
                        .format(self._parent.__name__, value_repr))

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
