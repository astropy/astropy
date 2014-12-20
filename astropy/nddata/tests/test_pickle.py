from ...extern.six.moves import cPickle as pickle
from ..nddata import NDData
from ...utils import NumpyRNGContext
from ... import coordinates as coord
from ...tests.helper import pytest, pickle_protocol, check_pickling_recovery
import numpy as np

originals = [NDData(np.random.random((10, 10)))]
xfails = [True]

@pytest.mark.parametrize("original,xfail",
                         zip(originals, xfails))
def test_simple_object(pickle_protocol, original, xfail):
    if xfail:
        pytest.xfail()
    check_pickling_recovery(original, pickle_protocol)
