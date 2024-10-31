# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect

import numpy as np
import pytest

from astropy import units as u
from astropy.nddata.decorators import support_nddata
from astropy.nddata.nddata import NDData
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import WCS


class CCDData(NDData):
    pass


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
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy

    data_out, wcs_out, unit_out = wrapped_function_1(
        data=data_in, wcs=wcs_in, unit=unit_in
    )

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


def test_pass_nddata():
    data_in = np.array([1, 2, 3])
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in)

    data_out, wcs_out, unit_out = wrapped_function_1(nddata_in)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


@pytest.mark.parametrize(
    "func",
    (
        lambda *, data, wcs, unit=None: (data, wcs, unit),
        lambda *, wcs=None, data, unit=None: (data, wcs, unit),
    ),
)
def test_pass_nddata_kwarg_only(func):
    wrapped_function = support_nddata(func)

    data_in = np.array([1, 2, 3])
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in)

    data_out, wcs_out, unit_out = wrapped_function(data=nddata_in)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


@pytest.mark.parametrize(
    "func, call_type",
    (
        pytest.param(
            lambda data, wcs=None, mask=None, /, *, unit=None: (data, wcs, unit, mask),
            "data_as_pos",
            id="data_wcs_mask_pos-only",
        ),
        pytest.param(
            lambda data, wcs=None, /, mask=None, *, unit: (data, wcs, unit, mask),
            "data_as_pos",
            id="data_wcs_pos-only",
        ),
        pytest.param(
            lambda wcs=None, /, data=None, mask=None, *, unit: (data, wcs, unit, mask),
            "data_as_kw",
            id="data_pos-or-kw",
        ),
        pytest.param(
            lambda wcs=None, /, mask=None, *, data, unit: (data, wcs, unit, mask),
            "data_as_kw",
            id="data_kw-only",
        ),
    ),
)
def test_pass_nddata_constrained_signature(func, call_type):
    wrapped_function = support_nddata(func)

    data_in = np.array([1, 2, 3])
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy
    mask_in = np.array([True, False, False])

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in, mask=mask_in)

    if call_type == "data_as_pos":
        args = (nddata_in,)
        kwargs = {}
    elif call_type == "data_as_kw":
        args = ()
        kwargs = {"data": nddata_in}

    data_out, wcs_out, unit_out, mask_out = wrapped_function(*args, **kwargs)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in
    assert mask_out is mask_in

    nddata2 = NDData(data_in, unit=unit_in, mask=mask_in)

    if call_type == "data_as_pos":
        args = (nddata2,)
        kwargs = {}
    elif call_type == "data_as_kw":
        args = ()
        kwargs = {"data": nddata2}

    data_out, wcs_out, unit_out, mask_out = wrapped_function(*args, **kwargs)

    assert data_out is data_in
    assert wcs_out is None
    assert unit_out is unit_in
    assert mask_out is mask_in

    if call_type == "data_as_pos":
        args = (nddata2, wcs_in)
        kwargs = {}
    elif call_type == "data_as_kw":
        args = (wcs_in,)
        kwargs = {"data": nddata2}

    data_out, wcs_out, unit_out, mask_out = wrapped_function(*args, **kwargs)

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in
    assert mask_out is mask_in

    if call_type == "data_as_pos":
        args = (nddata_in, wcs_in)
        kwargs = {}
    elif call_type == "data_as_kw":
        args = (wcs_in,)
        kwargs = {"data": nddata_in}

    with pytest.warns(
        AstropyUserWarning,
        match=(
            "Property wcs has been passed explicitly and as "
            "an NDData property, using explicitly specified value"
        ),
    ):
        data_out, wcs_out, unit_out, mask_out = wrapped_function(*args, **kwargs)
    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in
    assert mask_out is mask_in


def test_pass_nddata_and_explicit():
    data_in = np.array([1, 2, 3])
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy
    unit_in_alt = u.mJy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in)

    with pytest.warns(
        AstropyUserWarning,
        match=(
            "Property unit has been passed explicitly and as "
            "an NDData property, using explicitly specified value"
        ),
    ) as w:
        data_out, wcs_out, unit_out = wrapped_function_1(nddata_in, unit=unit_in_alt)
    assert len(w) == 1

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in_alt


def test_pass_nddata_ignored():
    data_in = np.array([1, 2, 3])
    wcs_in = WCS(naxis=1)
    unit_in = u.Jy

    nddata_in = NDData(data_in, wcs=wcs_in, unit=unit_in, mask=[0, 1, 0])

    with pytest.warns(
        AstropyUserWarning,
        match=(
            "The following attributes were set on the data "
            "object, but will be ignored by the function: mask"
        ),
    ) as w:
        data_out, wcs_out, unit_out = wrapped_function_1(nddata_in)
    assert len(w) == 1

    assert data_out is data_in
    assert wcs_out is wcs_in
    assert unit_out is unit_in


@pytest.mark.parametrize(
    "func",
    (
        lambda something, wcs=None, unit=None: None,
        lambda wcs=None, unit=None: None,
    ),
)
def test_incorrect_first_argument(func):
    with pytest.raises(
        ValueError, match="Can only wrap a function with a data argument"
    ):
        support_nddata(func)


def test_wrap_function_no_kwargs():
    @support_nddata
    def wrapped_function_5(data, other_data):
        return data

    data_in = np.array([1, 2, 3])
    nddata_in = NDData(data_in)

    assert wrapped_function_5(nddata_in, [1, 2, 3]) is data_in


def test_wrap_function_repack_valid():
    @support_nddata(repack=True, returns=["data"])
    def wrapped_function_5(data, other_data):
        return data

    data_in = np.array([1, 2, 3])
    nddata_in = NDData(data_in)

    nddata_out = wrapped_function_5(nddata_in, [1, 2, 3])

    assert isinstance(nddata_out, NDData)
    assert nddata_out.data is data_in


def test_wrap_function_accepts():
    class MyData(NDData):
        pass

    @support_nddata(accepts=MyData)
    def wrapped_function_5(data, other_data):
        return data

    data_in = np.array([1, 2, 3])
    nddata_in = NDData(data_in)
    mydata_in = MyData(data_in)

    assert wrapped_function_5(mydata_in, [1, 2, 3]) is data_in

    with pytest.raises(
        TypeError,
        match=(
            "Only NDData sub-classes that inherit "
            "from MyData can be used by this function"
        ),
    ):
        wrapped_function_5(nddata_in, [1, 2, 3])


def test_wrap_preserve_signature_docstring():
    @support_nddata
    def wrapped_function_6(data, wcs=None, unit=None):
        """
        An awesome function
        """

    if wrapped_function_6.__doc__ is not None:
        assert wrapped_function_6.__doc__.strip() == "An awesome function"

    signature = inspect.signature(wrapped_function_6)

    assert str(signature) == "(data, wcs=None, unit=None)"


def test_setup_failures1():
    # repack but no returns
    with pytest.raises(ValueError):
        support_nddata(repack=True)


def test_setup_failures2():
    # returns but no repack
    with pytest.raises(ValueError):
        support_nddata(returns=["data"])


def test_setup_failures9():
    # keeps but no repack
    with pytest.raises(ValueError):
        support_nddata(keeps=["unit"])


def test_setup_failures3():
    # same attribute in keeps and returns
    with pytest.raises(ValueError):
        support_nddata(repack=True, keeps=["mask"], returns=["data", "mask"])


def test_setup_failures4():
    # function accepts *args
    with pytest.raises(ValueError):

        @support_nddata
        def test(data, *args):
            pass


def test_setup_failures10():
    # function accepts **kwargs
    with pytest.raises(ValueError):

        @support_nddata
        def test(data, **kwargs):
            pass


def test_setup_failures5():
    # function accepts *args (or **kwargs)
    with pytest.raises(ValueError):

        @support_nddata
        def test(data, *args):
            pass


def test_setup_failures6():
    # First argument is not data
    with pytest.raises(ValueError):

        @support_nddata
        def test(img):
            pass


def test_setup_failures7():
    # accepts CCDData but was given just an NDData
    with pytest.raises(TypeError):

        @support_nddata(accepts=CCDData)
        def test(data):
            pass

        test(NDData(np.ones((3, 3))))


def test_setup_failures8():
    # function returns a different amount of arguments than specified. Using
    # NDData here so we don't get into troubles when creating a CCDData without
    # unit!
    with pytest.raises(ValueError):

        @support_nddata(repack=True, returns=["data", "mask"])
        def test(data):
            return 10

        test(NDData(np.ones((3, 3))))  # do NOT use CCDData here.


def test_setup_failures11():
    # function accepts no arguments
    with pytest.raises(ValueError):

        @support_nddata
        def test():
            pass


def test_setup_numpyarray_default():
    # It should be possible (even if it's not advisable to use mutable
    # defaults) to have a numpy array as default value.
    @support_nddata
    def func(data, wcs=np.array([1, 2, 3])):
        return wcs


def test_still_accepts_other_input():
    @support_nddata(repack=True, returns=["data"])
    def test(data):
        return data

    assert isinstance(test(NDData(np.ones((3, 3)))), NDData)
    assert isinstance(test(10), int)
    assert isinstance(test([1, 2, 3]), list)


def test_accepting_property_normal():
    # Accepts a mask attribute and takes it from the input
    @support_nddata
    def test(data, mask=None):
        return mask

    ndd = NDData(np.ones((3, 3)))
    assert test(ndd) is None
    ndd._mask = np.zeros((3, 3))
    assert np.all(test(ndd) == 0)
    # Use the explicitly given one (raises a Warning)
    with pytest.warns(AstropyUserWarning) as w:
        assert test(ndd, mask=10) == 10
    assert len(w) == 1


def test_parameter_default_identical_to_explicit_passed_argument():
    # If the default is identical to the explicitly passed argument this
    # should still raise a Warning and use the explicit one.
    @support_nddata
    def func(data, meta={"a": 1}):
        return meta

    with pytest.warns(AstropyUserWarning) as w:
        assert func(NDData(1, meta={"b": 2}), {"a": 1}) == {"a": 1}
    assert len(w) == 1

    assert func(NDData(1, meta={"b": 2})) == {"b": 2}


def test_accepting_property_notexist():
    # Accepts flags attribute but NDData doesn't have one
    @support_nddata
    def test(data, flags=10):
        return flags

    ndd = NDData(np.ones((3, 3)))
    test(ndd)


def test_accepting_property_translated():
    # Accepts a error attribute and we want to pass in uncertainty!
    @support_nddata(mask="masked")
    def test(data, masked=None):
        return masked

    ndd = NDData(np.ones((3, 3)))
    assert test(ndd) is None
    ndd._mask = np.zeros((3, 3))
    assert np.all(test(ndd) == 0)
    # Use the explicitly given one (raises a Warning)
    with pytest.warns(AstropyUserWarning) as w:
        assert test(ndd, masked=10) == 10
    assert len(w) == 1


def test_accepting_property_meta_empty():
    # Meta is always set (dict) so it has a special case that it's
    # ignored if it's empty but not None
    @support_nddata
    def test(data, meta=None):
        return meta

    ndd = NDData(np.ones((3, 3)))
    assert test(ndd) is None
    ndd._meta = {"a": 10}
    assert test(ndd) == {"a": 10}
