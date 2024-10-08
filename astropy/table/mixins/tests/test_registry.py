from copy import copy

import pytest
from numpy.testing import assert_equal

from astropy.table import Column, Table
from astropy.table.mixins.registry import (
    MixinRegistryError,
    _handlers,
    get_mixin_handler,
    register_mixin_handler,
)
from astropy.table.table_helpers import ArrayWrapper

ORIGINAL = {}


def setup_function(function):
    ORIGINAL["handlers"] = copy(_handlers)
    _handlers.clear()


def teardown_function(function):
    _handlers.clear()
    _handlers.update(ORIGINAL["handlers"])


class SpamData:
    pass


class SpamWrapper(ArrayWrapper):
    def __init__(self):
        super().__init__([0, 1, 3, 4, 5])


FULL_QUALNAME = "astropy.table.mixins.tests.test_registry.SpamData"


def handle_spam(obj):
    return SpamWrapper()


def handle_spam_alt(obj):
    return SpamWrapper()


def test_no_handler():
    data = SpamData()
    assert get_mixin_handler(data) is None


def test_register_handler():
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    assert get_mixin_handler(SpamData()) is handle_spam


def test_register_handler_override():
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    with pytest.raises(MixinRegistryError) as exc:
        register_mixin_handler(FULL_QUALNAME, handle_spam_alt)
    assert (
        exc.value.args[0]
        == "Handler for class astropy.table.mixins.tests.test_registry.SpamData is"
        " already defined"
    )
    register_mixin_handler(FULL_QUALNAME, handle_spam_alt, force=True)
    assert get_mixin_handler(SpamData()) is handle_spam_alt


def test_get_mixin_handler_str():
    # Check that we can also pass a fully qualified name to get_mixin_handler
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    assert get_mixin_handler(FULL_QUALNAME) is handle_spam


def test_add_column_to_empty_table():
    t = Table()
    t["a"] = SpamData()
    # By default, we get an object column.
    assert isinstance(t["a"], Column)
    assert t["a"].dtype == object
    # But after registration, we can get the intended mixin.
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    t = Table()
    t["a"] = SpamData()

    assert len(t) == 5
    assert isinstance(t["a"], SpamWrapper)
    assert_equal(t["a"].data, [0, 1, 3, 4, 5])


def test_add_column_to_existing_table():
    # As above, but for a table that already has a column
    # (addition used to depend on whether or a table was empty; gh-17102).
    t = Table([[5, 6, 7, 8, 9]], names=["x"])
    t["a"] = SpamData()
    assert isinstance(t["a"], Column)
    assert t["a"].dtype == object

    register_mixin_handler(FULL_QUALNAME, handle_spam)

    t["a"] = SpamData()

    assert len(t) == 5
    assert isinstance(t["a"], SpamWrapper)
    assert_equal(t["a"].data, [0, 1, 3, 4, 5])


def invalid_handler(obj):
    return "invalid"


def test_invalid_handler():
    t = Table()

    register_mixin_handler(FULL_QUALNAME, invalid_handler)

    with pytest.raises(TypeError) as exc:
        t["a"] = SpamData()
    assert (
        exc.value.args[0] == f"Mixin handler for object of type {FULL_QUALNAME} "
        "did not return a valid mixin column"
    )
