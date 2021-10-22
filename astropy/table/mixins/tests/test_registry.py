from copy import copy
from unittest.mock import MagicMock

import pytest
import numpy as np
from numpy.testing import assert_equal

from astropy.table import Table
from astropy.table.table_helpers import ArrayWrapper
from astropy.table.mixins.registry import (_handlers, register_mixin_handler,
                                           MixinRegistryError, get_mixin_handler)

ORIGINAL = {}


def setup_function(function):
    ORIGINAL['handlers'] = copy(_handlers)
    _handlers.clear()


def teardown_function(function):
    _handlers.clear()
    _handlers.update(ORIGINAL['handlers'])


class SpamData:
    pass


class SpamWrapper(ArrayWrapper):
    def __init__(self):
        super().__init__([0, 1, 3, 4, 5])


FULL_QUALNAME = 'astropy.table.mixins.tests.test_registry.SpamData'


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
    assert exc.value.args[0] == 'Handler for class astropy.table.mixins.tests.test_registry.SpamData is already defined'
    register_mixin_handler(FULL_QUALNAME, handle_spam_alt, force=True)
    assert get_mixin_handler(SpamData()) is handle_spam_alt


def test_get_mixin_handler_str():
    # Check that we can also pass a fully qualified name to get_mixin_handler
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    assert get_mixin_handler(FULL_QUALNAME) is handle_spam


def test_add_column():

    t = Table()
    with pytest.raises(TypeError):
        t['a'] = SpamData()

    register_mixin_handler(FULL_QUALNAME, handle_spam)

    t['a'] = SpamData()

    assert len(t) == 5
    assert isinstance(t['a'], SpamWrapper)
    assert_equal(t['a'].data, [0, 1, 3, 4, 5])


def invalid_handler(obj):
    return 'invalid'


def test_invalid_handler():

    t = Table()

    register_mixin_handler(FULL_QUALNAME, invalid_handler)

    with pytest.raises(TypeError) as exc:
        t['a'] = SpamData()
    assert exc.value.args[0] == (f'Mixin handler for object of type {FULL_QUALNAME} '
                                 f'did not return a valid mixin column')
