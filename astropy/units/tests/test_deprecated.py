# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

"""Regression tests for deprecated units."""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from .. import deprecated
from ... import units as u
from ...tests.helper import pytest  # TODO: Stop using bundled pytest


def test_emu():
    with pytest.raises(AttributeError):
        u.emu

    assert u.Bi.to(deprecated.emu, 1) == 1

    with deprecated.enable():
        assert u.Bi.compose()[0] == deprecated.emu

    assert u.Bi.compose()[0] == u.Bi
