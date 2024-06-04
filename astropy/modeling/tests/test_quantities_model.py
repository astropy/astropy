# Various tests of models not related to evaluation, fitting, or parameters
# pylint: disable=invalid-name, no-member
import pytest

from astropy import units as u
from astropy.modeling import models
from astropy.modeling.core import _ModelMeta
from astropy.modeling.models import Gaussian1D, Mapping, Pix2Sky_TAN
from astropy.tests.helper import assert_quantity_allclose


def test_gaussian1d_bounding_box():
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm, amplitude=3 * u.Jy)
    bbox = g.bounding_box.bounding_box()
    assert_quantity_allclose(bbox[0], 2.835 * u.m)
    assert_quantity_allclose(bbox[1], 3.165 * u.m)


def test_gaussian1d_n_models():
    g = Gaussian1D(
        amplitude=[1 * u.J, 2.0 * u.J],
        mean=[1 * u.m, 5000 * u.AA],
        stddev=[0.1 * u.m, 100 * u.AA],
        n_models=2,
    )
    assert_quantity_allclose(g(1.01 * u.m), [0.99501248, 0.0] * u.J)
    assert_quantity_allclose(
        g(u.Quantity([1.01 * u.m, 5010 * u.AA])), [0.99501248, 1.990025] * u.J
    )
    # FIXME: The following doesn't work as np.asanyarray doesn't work with a
    # list of quantity objects.
    # assert_quantity_allclose(g([1.01 * u.m, 5010 * u.AA]),
    #                            [ 0.99501248, 1.990025] * u.J)


"""
Test the "rules" of model units.
"""


def test_quantity_call():
    """
    Test that if constructed with Quanties models must be called with quantities.
    """
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm, amplitude=3 * u.Jy)

    g(10 * u.m)

    MESSAGE = (
        r".* Units of input 'x', .* could not be converted to required input units of"
        r" m .*"
    )
    with pytest.raises(u.UnitsError, match=MESSAGE):
        g(10)


def test_no_quantity_call():
    """
    Test that if not constructed with Quantities they can be called without quantities.
    """
    g = Gaussian1D(mean=3, stddev=3, amplitude=3)
    assert isinstance(g, Gaussian1D)
    g(10)


def test_default_parameters():
    # Test that calling with a quantity works when one of the parameters
    # defaults to dimensionless
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm)
    assert isinstance(g, Gaussian1D)
    g(10 * u.m)


def test_uses_quantity():
    """
    Test Quantity
    """
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm, amplitude=3 * u.Jy)

    assert g.uses_quantity

    g = Gaussian1D(mean=3, stddev=3, amplitude=3)

    assert not g.uses_quantity

    g.mean = 3 * u.m

    assert g.uses_quantity


def test_uses_quantity_compound():
    """
    Test Quantity
    """
    g = Gaussian1D(mean=3 * u.m, stddev=3 * u.cm, amplitude=3 * u.Jy)
    g2 = Gaussian1D(mean=5 * u.m, stddev=5 * u.cm, amplitude=5 * u.Jy)

    assert (g | g2).uses_quantity

    g = Gaussian1D(mean=3, stddev=3, amplitude=3)
    g2 = Gaussian1D(mean=5, stddev=5, amplitude=5)

    comp = g | g2

    assert not (comp).uses_quantity


def test_uses_quantity_no_param():
    comp = Mapping((0, 1)) | Pix2Sky_TAN()

    assert comp.uses_quantity


def _allmodels():
    allmodels = []
    for name in dir(models):
        model = getattr(models, name)
        if type(model) is _ModelMeta:
            try:
                m = model()
            except Exception:
                pass
            allmodels.append(m)
    return allmodels


@pytest.mark.parametrize("m", _allmodels())
def test_read_only(m):
    """
    input_units
    return_units
    input_units_allow_dimensionless
    input_units_strict
    """
    with pytest.raises(AttributeError):
        m.input_units = {}
    with pytest.raises(AttributeError):
        m.return_units = {}
    with pytest.raises(AttributeError):
        m.input_units_allow_dimensionless = {}
    with pytest.raises(AttributeError):
        m.input_units_strict = {}
