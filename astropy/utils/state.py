"""
A simple class to manage a piece of global science state.  See
:ref:`config-developer` for more details.
"""

import copy
from collections.abc import MappingView
from types import MappingProxyType

__all__ = ["ScienceState"]


class _IndexedMappingView(MappingView):
    """
    `~collections.abc.MappingView` with a read-only ``getitem`` through
    `~types.MappingProxyType`.

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


class ScienceStateContext:
    """
    Science state context. Returned by :func:`~ScienceState.set`.

    Parameters
    ----------
    parent : ClassType
    state : dict

    Examples
    --------
    ::

        class MyState(ScienceState):
            _value = None

        old_context = MyState.set("test")

        assert old_context._value is None
        assert MyState._value == "test"

    """

    def __init__(self, parent, state):
        super().__setattr__("_parent", parent)
        super().__setattr__("_old_state", state)

    @property
    def _value(cls):
        """Get value from parent ScienceState."""
        return cls._parent._value

    def __getattr__(self, key):
        """
        Get attribute from the class. Items in ``_state`` can also be accessed
        as attributes of the class.

        Raises
        ------
        AttributeError
            If key is not in class or ``_state``

        """
        try:
            return super().__getattribute__(key)
        except AttributeError:
            try:
                return super().__getattribute__("_parent")._state[key]
            except KeyError:
                raise AttributeError

    def __setattr__(cls, key, value):
        """Cannot set. Read-only."""
        raise NotImplementedError("Contexts are read-only.")

    def __enter__(self):
        """Enter ``with`` statement."""
        return self

    def __exit__(self, type, value, tb):
        """Exit ``with`` statement, resetting value and clearing context."""
        self._parent._state = self._old_state._mapping  # reset
        self.__dict__.clear()  # clears on exit

    def __repr__(self):
        # Ensure we have a single-line repr, just in case our
        # value is not something simple like a string.
        value_repr, lb, _ = repr(self._parent._value).partition("\n")
        if lb:
            value_repr += "..."
        return "<ScienceState {}: {}>".format(
            self._parent.__name__, value_repr
        )

    def get(self):
        """Get Value."""
        return _IndexedMappingView(self._parent._state)

    def set(self):
        raise NotImplementedError("Contexts are read-only.")


class _StateMixinMeta(type):
    """
    Metaclass for ScienceStates. Mixes ``_state`` into the class attributes.
    Ensures existence of ``_state`` attribute. ScienceStates generally define
    a ``_value`` attribute, which this metaclass moves to ``_state``.
    """

    @classmethod
    def __prepare__(cls, *args):
        """
        Prepare class namespace. Called before class attributes.
        the namespace is initialized with a ``_state`` field.
        """
        namespace = dict(_state=dict())
        return namespace

    def __new__(metacls, name, bases, namespace, **kwargs):
        """
        Creates a new instance of the class. Called after class attributes.
        ScienceStates generally define a ``_value`` attribute, which must
        be moved into the ``_state`` attribute.
        """
        # move _value to _state. This is for backward compatibility.
        if "_value" in namespace:
            key = namespace.get("_state_value_name_", "value")
            namespace["_state"][key] = namespace.pop("_value")

        return super().__new__(metacls, name, bases, namespace, **kwargs)

    def __getattr__(cls, key):
        """
        Get attribute from the class. Items in ``_state`` can also be accessed
        as attributes of the class.

        Raises
        ------
        AttributeError
            If key is not in class or ``_state``

        """
        try:  # first check if in class
            return super().__getattribute__(key)
        except AttributeError:  # it is not, check the _state
            try:
                return cls._state[key]
            except KeyError:  # match Exception to method.
                raise AttributeError
            # All other Exceptions are raised normally.

    def __setattr__(cls, key, value):
        """
        Set attribute to the class, if attribute already in class.
        Since ``_state`` holds most things of interest, this
        method should not change the state.

        To forcibly set, use ``super().__setattr__(key, value)`` or
        ``object.__setattr__(cls, key, value)``

        Raises
        ------
        AttributeError
            If key is not in class or ``_state``

        """
        try:
            super().__getattribute__(key)
        except Exception:
            raise NotImplementedError(
                "Can only set values to existing attributes."
            )
        else:  # Can only set if has
            super().__setattr__(key, value)

    @property
    def _value(cls):
        """Get value from ``_state``."""
        return cls._state[cls._state_value_name_]

    @_value.setter
    def _value(cls, value):
        """Set value in ``_state``."""
        cls._state[cls._state_value_name_] = value


class ScienceState(metaclass=_StateMixinMeta):
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

    _context_ = ScienceStateContext  #
    _state_value_name_ = "value"  # what cls._value gets from `_state`

    def __init__(self):
        raise RuntimeError("This class is a singleton. Do not instantiate.")

    def __repr__(self):
        # Ensure we have a single-line repr, just in case our
        # value is not something simple like a string.
        value_repr, lb, _ = repr(self._value).partition("\n")
        if lb:
            value_repr += "..."
        return "<ScienceState {}: {}>".format(self.__name__, value_repr)

    @classmethod
    def get(cls):
        """
        Get the current science state value.
        """
        return cls.validate(cls._value)

    @classmethod
    def get_state(cls):
        """
        Get view of current science state.
        """
        return _IndexedMappingView(cls._state)

    @classmethod
    def set(cls, value):
        """
        Set the current science state value.
        """
        # create context, storing old state
        ctx = cls._context_(cls, _IndexedMappingView(copy.deepcopy(cls._state)))  # FIXME

        # clear current state
        # cls._state = {}  # FIXME

        try:  # get new state
            value = cls.validate(value)
        except Exception as e:  # If problems, restore previous state
            cls._state = ctx._old_state._mapping
            raise e
        else:
            cls._state[
                cls._state_value_name_
            ] = value  # TODO replace by cls._value = value

        return ctx

    @classmethod
    def validate(cls, value):
        """
        Validate the value and convert it to its native type, if
        necessary.
        """
        return value
