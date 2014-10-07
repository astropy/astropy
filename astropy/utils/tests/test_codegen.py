# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import traceback

from ..codegen import make_function_with_signature
from ...tests.helper import pytest


def test_make_function_with_signature_lineno():
    """
    Tests that a function made with ``make_function_with_signature`` is give
    the correct line number into the module it was created from (i.e. the line
    ``make_function_with_signature`` was called from).
    """

    def crashy_function(*args, **kwargs):
        1 / 0

    # Make a wrapper around this function with the signature:
    # crashy_function(a, b)
    # Note: the signature is not really relevant to this test
    wrapped = make_function_with_signature(crashy_function, ('a', 'b'))
    line = """
    wrapped = make_function_with_signature(crashy_function, ('a', 'b'))
    """.strip()

    try:
        wrapped(1, 2)
    except:
        exc_cls, exc, tb = sys.exc_info()
        assert exc_cls is ZeroDivisionError
        # The *last* line in the traceback should be the 1 / 0 line in
        # crashy_function; the next line up should be the line that the
        # make_function_with_signature call was one
        tb_lines = traceback.format_tb(tb)
        assert '1 / 0' in tb_lines[-1]
        assert line in tb_lines[-2] and 'line =' not in tb_lines[-2]
    else:
        pytest.fail('This should have caused an exception')
