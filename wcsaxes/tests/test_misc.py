import matplotlib.pyplot as plt
from ..core import WCSAxes


def test_grid_regression():
    # Regression test for a bug that meant that if the rc parameter
    # axes.grid was set to True, WCSAxes would crash upon initalization.
    plt.rc('axes', grid=True)
    fig = plt.figure(figsize=(3, 3))
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8])
