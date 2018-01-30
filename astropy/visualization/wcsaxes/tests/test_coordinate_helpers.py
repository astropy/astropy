# Licensed under a 3-clause BSD style license - see LICENSE.rst
import matplotlib.pyplot as plt

from ..core import WCSAxes
from ....tests.image_tests import ignore_matplotlibrc


@ignore_matplotlibrc
def test_getaxislabel():

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], aspect='equal')

    ax.coords[0].set_axislabel("X")
    ax.coords[1].set_axislabel("Y")
    assert ax.coords[0].get_axislabel() == "X"
    assert ax.coords[1].get_axislabel() == "Y"
