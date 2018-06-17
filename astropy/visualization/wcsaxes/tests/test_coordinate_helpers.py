# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from unittest.mock import patch

import pytest
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits

from ..core import WCSAxes
from .... import units as u
from ....tests.image_tests import ignore_matplotlibrc

ROOT = os.path.join(os.path.dirname(__file__))
MSX_HEADER = fits.Header.fromtextfile(os.path.join(ROOT, 'data', 'msx_header'))


@ignore_matplotlibrc
def test_getaxislabel():

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal')

    ax.coords[0].set_axislabel("X")
    ax.coords[1].set_axislabel("Y")
    assert ax.coords[0].get_axislabel() == "X"
    assert ax.coords[1].get_axislabel() == "Y"


@pytest.fixture
def ax():
    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal')
    fig.add_axes(ax)

    return ax


def assert_label_draw(ax, x_label, y_label):
    ax.coords[0].set_axislabel("Label 1")
    ax.coords[1].set_axislabel("Label 2")

    with patch.object(ax.coords[0].axislabels, 'set_position') as pos1:
        with patch.object(ax.coords[1].axislabels, 'set_position') as pos2:
            ax.figure.canvas.draw()

    assert pos1.call_count == x_label
    assert pos2.call_count == y_label


@ignore_matplotlibrc
def test_label_visibility_rules_default(ax):
    assert_label_draw(ax, True, True)


@ignore_matplotlibrc
def test_label_visibility_rules_label(ax):

    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks(values=[-9999]*u.one)

    assert_label_draw(ax, False, False)


@ignore_matplotlibrc
def test_label_visibility_rules_ticks(ax):

    ax.coords[0].set_axislabel_visibility_rule('ticks')
    ax.coords[1].set_axislabel_visibility_rule('ticks')

    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks(values=[-9999]*u.one)

    assert_label_draw(ax, True, False)


@ignore_matplotlibrc
def test_label_visibility_rules_always(ax):

    ax.coords[0].set_axislabel_visibility_rule('always')
    ax.coords[1].set_axislabel_visibility_rule('always')

    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks(values=[-9999]*u.one)

    assert_label_draw(ax, True, True)


def test_set_separator(tmpdir):

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS(MSX_HEADER))
    fig.add_axes(ax)

    # Force a draw which is required for format_coord to work
    ax.figure.canvas.draw()

    ax.coords[1].set_format_unit('deg')
    assert ax.coords[1].format_coord(4) == '4\xb000\'00\"'
    ax.coords[1].set_separator((':', ':', ''))
    assert ax.coords[1].format_coord(4) == '4:00:00'
    ax.coords[1].set_separator('abc')
    assert ax.coords[1].format_coord(4) == '4a00b00c'
    ax.coords[1].set_separator(None)
    assert ax.coords[1].format_coord(4) == '4\xb000\'00\"'
