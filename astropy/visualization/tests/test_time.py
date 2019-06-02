# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True

from astropy.time import Time
from astropy.visualization.time import time_support


def get_ticklabels(axis):
    axis.figure.canvas.draw()
    return [x.get_text() for x in axis.get_ticklabels()]


# We first check that we get the expected labels for different time intervals
# for standard ISO formatting. This is a way to check both the locator and
# formatter code.

RANGE_CASES = [

    # Interval of many years
    (('2014-03-22T12:30:30.9', '2077-03-22T12:30:32.1'),
     ['2020-01-01',
      '2040-01-01',
      '2060-01-01']),

    # Interval of a few years
    (('2014-03-22T12:30:30.9', '2017-03-22T12:30:32.1'),
     ['2015-01-01',
      '2016-01-01',
      '2017-01-01']),

    # Interval of just under a year
    (('2014-03-22T12:30:30.9', '2015-01-22T12:30:32.1'),
     ['2014-05-01',
      '2014-10-01']),

    # Interval of just over a month
    (('2014-03-22T12:30:30.9', '2014-04-23T12:30:32.1'),
     ['2014-04-01']),

    # Interval of just under a month
    (('2014-03-22T12:30:30.9', '2014-04-21T12:30:32.1'),
     ['2014-03-14',  # FIXME: incorrect!
      '2014-03-24',
      '2014-04-03',
      '2014-04-13',
      '2014-04-23']),

    # Interval of just over an hour
    (('2014-03-22T12:30:30.9', '2014-03-22T13:31:30.9'),
     ['2014-03-22T12:40:00.000',
      '2014-03-22T13:00:00.000',
      '2014-03-22T13:20:00.000']),

    # Interval of just under an hour
    (('2014-03-22T12:30:30.9', '2014-03-22T13:28:30.9'),
     ['2014-03-22T12:40:00.000',
      '2014-03-22T13:00:00.000',
      '2014-03-22T13:20:00.000']),

    # Interval of a few minutes
    (('2014-03-22T12:30:30.9', '2014-03-22T12:38:30.9'),
     ['2014-03-22T12:33:00.000',
      '2014-03-22T12:36:00.000']),

    # Interval of a few seconds
    (('2014-03-22T12:30:30.9', '2014-03-22T12:30:40.9'),
     ['2014-03-22T12:30:33.000',
      '2014-03-22T12:30:36.000',
      '2014-03-22T12:30:39.000']),

    # Interval of a couple of seconds
    (('2014-03-22T12:30:30.9', '2014-03-22T12:30:32.1'),
     ['2014-03-22T12:30:31.000',
      '2014-03-22T12:30:31.500',
      '2014-03-22T12:30:32.000']),

    # Interval of under a second
    (('2014-03-22T12:30:30.89', '2014-03-22T12:30:31.19'),
     ['2014-03-22T12:30:30.900',
      '2014-03-22T12:30:31.000',
      '2014-03-22T12:30:31.100']),

]


@pytest.mark.skipif('not HAS_PLT')
@pytest.mark.parametrize(('interval', 'expected'), RANGE_CASES)
def test_formatter_locator(interval, expected):

    # Check that the ticks and labels returned for the above cases are correct.

    with time_support():
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(Time(interval[0]), Time(interval[1]))
        assert get_ticklabels(ax.xaxis) == expected


FORMAT_CASES = [
  ('byear', ['2000', '2020', '2040', '2060', '2080']),
  ('byear_str', ['B2000.000', 'B2020.000', 'B2040.000', 'B2060.000', 'B2080.000']),
  ('cxcsec', ['500000000', '1000000000', '1500000000', '2000000000', '2500000000', '3000000000']),
  ('decimalyear', ['2000', '2020', '2040', '2060', '2080']),
  ('fits', ['2020-01-01T00:00:00.000', '2040-01-01T00:00:00.000', '2060-01-01T00:00:00.000']),
  ('gps', ['1000000000', '1500000000', '2000000000', '2500000000', '3000000000', '3500000000']),
  ('iso', ['2020-01-01 00:00:00.000', '2040-01-01 00:00:00.000', '2060-01-01 00:00:00.000']),
  ('isot', ['2020-01-01T00:00:00.000', '2040-01-01T00:00:00.000', '2060-01-01T00:00:00.000']),
  ('jd', ['2452000', '2458000', '2464000', '2470000', '2476000', '2482000']),
  ('jyear', ['2000', '2020', '2040', '2060', '2080']),
  ('jyear_str', ['J2000.000', 'J2020.000', 'J2040.000', 'J2060.000', 'J2080.000']),
  ('mjd', ['54000', '60000', '66000', '72000', '78000', '84000']),
  ('plot_date', ['732000', '738000', '744000', '750000', '756000', '762000']),
  ('unix', ['1000000000', '1500000000', '2000000000', '2500000000', '3000000000', '3500000000']),
  ('yday', ['2020:001:00:00:00.000', '2040:001:00:00:00.000', '2060:001:00:00:00.000']),
]


@pytest.mark.parametrize(('format', 'expected'), FORMAT_CASES)
def test_formats(format, expected):
    # Check that the locators/formatters work fine for all time formats
    with time_support(format=format, simplify=False):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(Time('2014-03-22T12:30:30.9'), Time('2077-03-22T12:30:32.1'))
        assert get_ticklabels(ax.xaxis) == expected
        ax.get_xlabel() == f'Time ({format})'


@pytest.mark.parametrize(('format', 'expected'), FORMAT_CASES)
def test_auto_formats(format, expected):
    # Check that the format/scale is taken from the first time used.
    with time_support(simplify=False):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(Time(Time('2014-03-22T12:30:30.9'), format=format),
                    Time('2077-03-22T12:30:32.1'))
        assert get_ticklabels(ax.xaxis) == expected


FORMAT_CASES_SIMPLIFY = [
  ('fits', ['2020-01-01', '2040-01-01', '2060-01-01']),
  ('iso', ['2020-01-01', '2040-01-01', '2060-01-01']),
  ('isot', ['2020-01-01', '2040-01-01', '2060-01-01']),
  ('yday', ['2020', '2040', '2060']),
]


@pytest.mark.parametrize(('format', 'expected'), FORMAT_CASES_SIMPLIFY)
def test_formats_simplify(format, expected):
    # Check the use of the simplify= option
    with time_support(format=format, simplify=True):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(Time('2014-03-22T12:30:30.9'), Time('2077-03-22T12:30:32.1'))
        assert get_ticklabels(ax.xaxis) == expected
