# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re

import numpy as np
import pytest

from astropy.time import Time, conf, TimeYearDayTime

iso_times = ['2000-02-29', '1981-12-31 12:13', '1981-12-31 12:13:14', '2020-12-31 12:13:14.56']
isot_times = [re.sub(' ', 'T', tm) for tm in iso_times]
yday_times = ['2000:060', '1981:365:12:13:14', '1981:365:12:13', '2020:366:12:13:14.56']


def test_fast_conf():
    # Default is to try C parser and then Python parser. Both fail so we get the
    # Python message.
    assert conf.use_fast_parser == 'True'  # default
    with pytest.raises(ValueError, match='Time 2000:0601 does not match yday format'):
        Time('2000:0601', format='yday')

    # This is one case where Python parser is different from C parser because the
    # Python parser has a bug and fails with a trailing ".", but C parser works.
    Time('2020:150:12:13:14.', format='yday')
    with conf.set_temp('use_fast_parser', 'force'):
        Time('2020:150:12:13:14.', format='yday')
    with conf.set_temp('use_fast_parser', 'False'):
        with pytest.raises(ValueError, match='could not convert string to float'):
            Time('2020:150:12:13:14.', format='yday')

    with conf.set_temp('use_fast_parser', 'False'):
        assert conf.use_fast_parser == 'False'
        # Make sure that this is really giving the Python parser
        with pytest.raises(ValueError, match='Time 2000:0601 does not match yday format'):
            Time('2000:0601', format='yday')

    with conf.set_temp('use_fast_parser', 'force'):
        assert conf.use_fast_parser == 'force'
        # Make sure that this is really giving the Python parser
        err = 'fast C time string parser failed: time string ends in middle of component'
        with pytest.raises(ValueError, match=err):
            Time('2000:0601', format='yday')


@pytest.mark.parametrize('times,format', [(iso_times, 'iso'),
                                          (isot_times, 'isot'),
                                          (yday_times, 'yday')])
@pytest.mark.parametrize('variant', [0, 1, 2])
def test_fast_matches_python(times, format, variant):
    if variant == 0:
        pass  # list/array of different values (null terminated strings included)
    elif variant == 1:
        times = times[-1]  # scalar
    elif variant == 2:
        times = [times[-1]] * 2  # list/array of identical values (no null terminations)

    with conf.set_temp('use_fast_parser', 'False'):
        tms_py = Time(times, format=format)

    with conf.set_temp('use_fast_parser', 'force'):
        tms_c = Time(times, format=format)

    # Times are binary identical
    assert np.all(tms_py == tms_c)


def test_fast_yday_exceptions():
    # msgs = {1: 'time string ends at beginning of component where break is not allowed',
    #         2: 'time string ends in middle of component',
    #         3: 'required delimiter character not found',
    #         4: 'non-digit found where digit (0-9) required',
    #         5: 'bad day of year (1 <= doy <= 365 or 366 for leap year'}

    with conf.set_temp('use_fast_parser', 'force'):
        for times, err in [('2020:150:12', 'time string ends at beginning of component'),
                           ('2020:150:1', 'time string ends in middle of component'),
                           ('2020:150*12:13:14', 'required delimiter character'),
                           ('2020:15*:12:13:14', 'non-digit found where digit'),
                           ('2020:999:12:13:14', 'bad day of year')]:
            with pytest.raises(ValueError, match=err):
                Time(times, format='yday')


def test_fast_iso_exceptions():
    with conf.set_temp('use_fast_parser', 'force'):
        for times, err in [('2020-10-10 12', 'time string ends at beginning of component'),
                           ('2020-10-10 1', 'time string ends in middle of component'),
                           ('2020*10-10 12:13:14', 'required delimiter character'),
                           ('2020-10-10 *2:13:14', 'non-digit found where digit')]:
            with pytest.raises(ValueError, match=err):
                Time(times, format='iso')


def test_fast_non_ascii():
    with pytest.raises(ValueError, match='input is not pure ASCII'):
        with conf.set_temp('use_fast_parser', 'force'):
            Time('2020-01-01 1ᛦ:13:14.4324')


def test_fast_subclass():
    """Test subclass where use_fast_parser class attribute is not in __dict__"""
    class TimeYearDayTimeSubClass(TimeYearDayTime):
        name = 'yday_subclass'

    # Inheritance works
    assert hasattr(TimeYearDayTimeSubClass, 'fast_parser_pars')
    assert 'fast_parser_pars' not in TimeYearDayTimeSubClass.__dict__

    try:
        # For YearDayTime, forcing the fast parser with a bad date will give
        # "fast C time string parser failed: time string ends in middle of component".
        # But since YearDayTimeSubClass does not have fast_parser_pars it will
        # use the Python parser.
        with pytest.raises(ValueError, match='Time 2000:0601 does not match yday_subclass format'):
            with conf.set_temp('use_fast_parser', 'force'):
                Time('2000:0601', format='yday_subclass')
    finally:
        del TimeYearDayTimeSubClass._registry['yday_subclass']
