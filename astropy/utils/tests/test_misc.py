# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
import os

import numpy as np

from .. import data, misc
from ...tests.helper import remote_data
from ...extern import six


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

    assert Subclass.__call__.__doc__ == "FOO"


def test_gen_ref_comments():
    comments = misc._generate_referee_comments(seed=1234, num=5)
    np.testing.assert_array_equal(
        comments,
        ['I cannot get past the abstract.',
         'Change your font to MS Comic Sans.',
         'No.',
         'Replace Figure X with a cat pic.',
         'Cite my paper.'])


def test_gen_ref_response():
    response = misc._generate_referee_response(seed=1234, num=5)
    np.testing.assert_array_equal(
        response, ['Maybe.', 'No.', 'No.', 'No.', 'No.'])
