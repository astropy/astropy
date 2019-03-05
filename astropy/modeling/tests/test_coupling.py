import numpy as np
import pytest
from astropy.modeling.models import (Shift, Identity, Pix2Sky_AZP,
                                     Pix2Sky_CylindricalPerspective)
from astropy.modeling.coupled import coupling_matrix


@pytest.mark.parametrize("model", ((Pix2Sky_CylindricalPerspective()),
                                   (Pix2Sky_AZP())))
def test_invaritant_separable_model(model):
    sep = coupling_matrix(model)

    assert sep.all()


def test_compound_model():
    comp = coupling_matrix(Pix2Sky_CylindricalPerspective() & (Identity(1) | Shift(10)))

    answer = np.array([[True, True, False],
                       [True, True, False],
                       [False, False, True]])

    np.testing.assert_array_equal(comp, answer)
