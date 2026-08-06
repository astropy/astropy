import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Path

from astropy.visualization.wcsaxes import conf
from astropy.visualization.wcsaxes.grid_paths import get_lon_lat_path
from astropy.wcs import WCS


def test_gridline_two_samples_no_error():
    # Regression test for a bug where drawing a gridline sampled with only
    # two points raised an IndexError while checking for discontinuities.
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [0, 0]
    wcs.wcs.crpix = [1, 1]
    wcs.wcs.cdelt = [-1, 1]

    fig = Figure()
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    original_grid_samples = conf.grid_samples
    conf.grid_samples = 2
    try:
        ax.coords[0].grid()
        canvas.draw()
    finally:
        conf.grid_samples = original_grid_samples


@pytest.mark.parametrize("step_in_degrees", [10, 1, 0.01])
def test_round_trip_visibility(step_in_degrees):
    zero = np.zeros(100)

    # The pixel values are irrelevant for this test
    pixel = np.stack([zero, zero]).T

    # Create a grid line of constant latitude with a point every step
    line = np.stack([np.arange(100), zero]).T * step_in_degrees

    # Create a modified grid line where the point spacing is larger by 5%
    # Starting with point 20, the discrepancy between `line` and `line_round` is greater than `step`
    line_round = line * 1.05

    # Perform the round-trip check
    path = get_lon_lat_path(line, pixel, line_round)

    # The grid line should be visible for only the initial part line (19 points)
    codes_check = np.full(100, Path.MOVETO)
    codes_check[line_round[:, 0] - line[:, 0] < step_in_degrees] = Path.LINETO
    assert np.all(path.codes[1:] == codes_check[1:])
