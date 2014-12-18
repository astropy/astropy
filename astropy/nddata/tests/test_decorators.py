# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect

import numpy as np

from ...tests.helper import catch_warnings, pytest
from ... import units as u

from ..nddata import NDData
from ..decorators import support_nddata


@support_nddata
def wrapped_function_1(data, wcs=None, unit=None):
    return data, wcs, unit


def test_pass_numpy():

    data_in = np.array([1, 2, 3])
    data_out, wcs_out, unit_out = wrapped_function_1(data=data_in)

    assert data_out is data_in
    assert wcs_out is None
    assert unit_out is None


def test_pass_all_separate():

    data_in = np.array([1, 2, 3])
    wcs_in = "the wcs"
    unit_in = u.Jy

    data_out, wcs_out, unit_out = wrapped_function_1(data=data_in, wcs=wcs_in, unit=unit_in)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


def test_pass_nddata():

    data_in = np.array([1, 2, 3])
    wcs_in = "the wcs"
    unit_in = u.Jy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in)

    data_out, wcs_out, unit_out = wrapped_function_1(nddata_in)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


def test_pass_nddata_and_explicit():

    data_in = np.array([1, 2, 3])
    wcs_in = "the wcs"
    unit_in = u.Jy
    unit_in_alt = u.mJy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in)

    with catch_warnings() as w:
        data_out, wcs_out, unit_out = wrapped_function_1(nddata_in, unit=unit_in_alt)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in_alt

    assert len(w) == 1
    assert str(w[0].message) == ("Property unit has been passed explicitly and as "
                                 "an NDData property, using explicitly specified value")


def test_pass_nddata_ignored():

    data_in = np.array([1, 2, 3])
    wcs_in = "the wcs"
    unit_in = u.Jy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in, mask=[0, 1, 0])

    with catch_warnings() as w:
        data_out, wcs_out, unit_out = wrapped_function_1(nddata_in)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in

    assert len(w) == 1
    assert str(w[0].message) == ("The following attributes were set on the data "
                                 "object, but will be ignored by the function: mask")


def test_incorrect_first_argument():

    with pytest.raises(ValueError) as exc:
        @support_nddata
        def wrapped_function_2(something, wcs=None, unit=None):
            pass
    assert exc.value.args[0] == "Can only wrap functions whose first positional argument is `data`"

    with pytest.raises(ValueError) as exc:
        @support_nddata
        def wrapped_function_3(something, data, wcs=None, unit=None):
            pass
    assert exc.value.args[0] == "Can only wrap functions whose first positional argument is `data`"

    with pytest.raises(ValueError) as exc:
        @support_nddata
        def wrapped_function_4(wcs=None, unit=None):
            pass
    assert exc.value.args[0] == "Can only wrap functions whose first positional argument is `data`"


def test_wrap_function_no_kwargs():

    @support_nddata
    def wrapped_function_5(data, other_data):
        return data

    data_in = np.array([1, 2, 3])
    nddata_in = NDData(data_in)

    assert wrapped_function_5(nddata_in, [1, 2, 3]) is data_in


def test_wrap_function_repack_valid():

    @support_nddata(repack=True, returns=['data'])
    def wrapped_function_5(data, other_data):
        return data

    data_in = np.array([1, 2, 3])
    nddata_in = NDData(data_in)

    nddata_out = wrapped_function_5(nddata_in, [1, 2, 3])

    assert isinstance(nddata_out, NDData)
    assert nddata_out.data is data_in


def test_wrap_preserve_signature_docstring():

    @support_nddata
    def wrapped_function_6(data, wcs=None, unit=None):
        """
        An awesome function
        """
        pass

    assert wrapped_function_6.__doc__.strip() == "An awesome function"

    signature = inspect.formatargspec(*inspect.getargspec(wrapped_function_6))

    assert signature == "(data, wcs=None, unit=None)"
