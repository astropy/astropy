from copy import copy
from unittest.mock import MagicMock

import numpy as np

from astropy.table.mixins.registry import _handlers

ORIGINAL = {}


def setup_function(function):
    ORIGINAL['handlers'] = copy(_handlers)
    _handlers.clear()


def teardown_function(function):
    _handlers.clear()
    _handlers.update(ORIGINAL['handlers'])


class SpamData:
    def __array__(self):
        return np.arange(5)


FULL_QUALNAME = 'astropy.table.mixins.tests.test_registry.SpamData'


def handle_spam(obj):
    return ArrayWrapper(obj)


def test_no_handler():
    data = SpamData()
    assert mixin_handler(data) is None


def test_register_handler():
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    data = SpamData()
    assert mixin_handler(array) is handle_spam


def test_register_handler_override():
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    register_mixin_handler(FULL_QUALNAME, handle_spam)
    register_mixin_handler(FULL_QUALNAME, handle_spam)


def test_table_operations():

    t = Table()
    t['a'] = SpamData()

    register_mixin_handler(FULL_QUALNAME, handle_spam)

    t['a'] = SpamData()

    assert len(t) == 5
    assert isinstance(t['a'], ArrayWrapper)
    assert_equal(t['a'], [0, 1, 2, 3, 4])

    subt = t[1:4]
    assert_equal(t['a'], [0, 1, 2, 3, 4])

