# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

import pytest

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True

from astropy.time import Time
from astropy.visualization.time import time_support


#
# converter = AstropyTimeConverter(format='mjd')
# units.registry[Time] = converter
#
# times = Time(['2016-03-22T12:30:31', '2016-03-22T12:30:38', '2016-03-22T12:34:40'])
#
# converter.format = 'iso'
#
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# ax.plot(times, [2, 3, 5])
# fig.savefig('test_times.pdf')
# print([x.get_text() for x in ax.xaxis.get_ticklabels()])


def get_ticklabels(axis):
    axis.figure.canvas.draw()
    return [x.get_text() for x in axis.get_ticklabels()]


# We first check that we get the expected labels for different time intervals
# for standard ISO formatting. This is a way to check both the locator and
# formatter code.

CASES = [

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
@pytest.mark.parametrize(('interval', 'expected'), CASES)
def test_formatter_locator(interval, expected):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    with time_support():
        ax.set_xlim(Time(interval[0]), Time(interval[1]))
        assert get_ticklabels(ax.xaxis) == expected
