# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
import os
from datetime import datetime
import locale

import numpy as np

from .. import data, misc
from ...tests.helper import remote_data
from ...extern import six
from ...tests.helper import pytest

try:
    locale.setlocale(locale.LC_ALL, 'en_US')
    locale.setlocale(locale.LC_ALL, 'de_DE')
except:
    HAS_LOCALES = False
else:
    HAS_LOCALES = True


def test_isiterable():
    assert misc.isiterable(2) is False
    assert misc.isiterable([2]) is True
    assert misc.isiterable([1, 2, 3]) is True
    assert misc.isiterable(np.array(2)) is False
    assert misc.isiterable(np.array([1, 2, 3])) is True


@remote_data
def test_api_lookup():
    strurl = misc.find_api_page('astropy.utils.misc', 'dev', False, timeout=3)
    objurl = misc.find_api_page(misc, 'dev', False, timeout=3)

    assert strurl == objurl
    assert strurl == 'http://devdocs.astropy.org/utils/index.html#module-astropy.utils.misc'


def test_skip_hidden():
    path = data._find_pkg_data_path('data')
    for root, dirs, files in os.walk(path):
        assert '.hidden_file.txt' in files
        assert 'local.dat' in files
        # break after the first level since the data dir contains some other
        # subdirectories that don't have these files
        break

    for root, dirs, files in misc.walk_skip_hidden(path):
        assert '.hidden_file.txt' not in files
        assert 'local.dat' in files
        break


def test_JsonCustomEncoder():
    assert json.dumps(np.arange(3), cls=misc.JsonCustomEncoder) == '[0, 1, 2]'
    assert json.dumps(1+2j, cls=misc.JsonCustomEncoder) == '[1.0, 2.0]'
    assert json.dumps(set([1, 2, 1]), cls=misc.JsonCustomEncoder) == '[1, 2]'
    assert json.dumps(b'hello world \xc3\x85',
                      cls=misc.JsonCustomEncoder) == '"hello world \\u00c5"'
    assert json.dumps({1: 2},
                      cls=misc.JsonCustomEncoder) == '{"1": 2}'  # default


def test_inherit_docstrings():
    @six.add_metaclass(misc.InheritDocstrings)
    class Base(object):
        def __call__(self, *args):
            "FOO"
            pass

    class Subclass(Base):
        def __call__(self, *args):
            pass

    if Base.__call__.__doc__ is not None:
        # TODO: Maybe if __doc__ is None this test should be skipped instead?
        assert Subclass.__call__.__doc__ == "FOO"


@pytest.mark.skipif('not HAS_LOCALES')
def test_set_locale():
    date = datetime(2000, 10, 1, 0, 0, 0)
    day_mon = date.strftime('%a, %b')

    with misc.set_locale('en_US'):
        assert date.strftime('%a, %b') == 'Sun, Oct'

    with misc.set_locale('de_DE'):
        assert date.strftime('%a, %b') == 'So, Okt'

    # Back to original
    assert date.strftime('%a, %b') == day_mon

    current = locale.setlocale(locale.LC_ALL)
    with misc.set_locale(current):
        assert date.strftime('%a, %b') == day_mon
