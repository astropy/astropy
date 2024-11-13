# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.state import ScienceState


def test_ScienceState_and_Context():
    """
    Tests a ScienceState and spawned contexts.
    """

    class MyState(ScienceState):
        _value = "A"
        _state = dict(foo="bar")

    state = {"foo": "bar"}

    # test created ScienceState
    assert MyState.get() == "A"
    assert MyState.validate("B") == "B"
    assert MyState._state == state

    # test setting
    with MyState.set("B"):
        assert MyState.get() == "B"
    assert MyState.get() == "A"  # test returning
