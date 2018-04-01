# Various tests of models not related to evaluation, fitting, or parameters

from ...tests.helper import assert_quantity_allclose
from ... import units as u
from ..functional_models import Gaussian1D


def test_gaussian1d_bounding_box():
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm, amplitude=3 * u.Jy)
    bbox = g.bounding_box
    assert_quantity_allclose(bbox[0], 2.835 * u.m)
    assert_quantity_allclose(bbox[1], 3.165 * u.m)


def test_gaussian1d_n_models():
    g = Gaussian1D(amplitude=[1 * u.J, 2. * u.J],
                   mean=[1 * u.m, 5000 * u.AA],
                   stddev=[0.1 * u.m, 100 * u.AA],
                   n_models=2)
    assert_quantity_allclose(g(1.01 * u.m), [0.99501248, 0.] * u.J)
    assert_quantity_allclose(g(u.Quantity([1.01 * u.m, 5010 * u.AA])), [0.99501248, 1.990025] * u.J)
    # FIXME: The following doesn't work as np.asanyarray doesn't work with a
    # list of quantity objects.
    # assert_quantity_allclose(g([1.01 * u.m, 5010 * u.AA]),
    #                            [ 0.99501248, 1.990025] * u.J)
