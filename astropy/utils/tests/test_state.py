# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.state import ScienceState
import pytest


def test_ScienceState_and_Context():
    """
    Tests a ScienceState and spawned contexts.
    """
    class MyState(ScienceState):
        _value = "A"
        _state = dict(foo="bar")

        @classmethod
        def validate(cls, value):
            if value not in ('A', 'B', 'C'):
                raise ValueError("Must be one of A, B, C")
            return value

    assert MyState.get() == "A"
    assert MyState.validate("B") == "B"
    assert MyState._state == dict(foo="bar", value="A")

    old_context = MyState.set("B")
    assert old_context._parent is MyState
    assert old_context._value == "B"

    with old_context:
        assert MyState._value == "B"
        assert MyState._state == dict(foo="bar", value="B")

    # test context descopes
    with pytest.raises(AttributeError):
        old_context._value
