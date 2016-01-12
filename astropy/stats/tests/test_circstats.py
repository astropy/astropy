from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from numpy.testing import assert_equal
from numpy.testing.utils import assert_allclose

try:
    import scipy.stats
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

from ...tests.helper import pytest
from ..circstats import _length, circmean, circvar, circmoment, circcorrcoef
from ..circstats import rayleightest, vtest, vonmisesmle


def test__length():
    # testing against R CircStats package
    # Ref. [1] pages 6 and 125
    data = np.deg2rad(np.array([0, 3.6, 36, 72, 108, 169.2, 324]))
    w = np.array([12, 1, 6, 1, 2, 1, 1])
    answer = 0.766282

    assert_allclose(answer, _length(data, w), atol=1e-4)

    # testing against scipy.stats.circmean function
    data = np.random.uniform(-np.pi, np.pi, size=(100, 5))
    answer = scipy.stats.circmean(data, np.pi, -np.pi)
    assert_allclose(answer, circmean(data), atol=1e-6)


def test_circmean():
    # testing against R CircStats package
    # Ref[1], page 23
    data = np.deg2rad(np.array([51, 67, 40, 109, 31, 358]))
    answer = 0.8487044
    assert_allclose(answer, circmean(data), atol=1e-4)


def test_circvar():
    # testing against R CircStats package
    # Ref[1], page 23
    data = np.deg2rad(np.array([51, 67, 40, 109, 31, 358]))
    answer = 0.1635635
    assert_allclose(answer, circvar(data), atol=1e-4)


def test_circmoment():
    # testing against R CircStats package
    # Ref[1], page 23
    data = np.deg2rad(np.array([51, 67, 40, 109, 31, 358]))
    # 2nd, 3rd, and 4th moments
    answer = np.array([[1.588121, 0.4800428], [1.963919, 0.236541],
                       [2.685556, 0.2255761]])
    assert_allclose(answer, (circmoment(data, p=2),
                             circmoment(data, p=3),
                             circmoment(data, p=4)), atol=1e-4)


def test_circcorrcoef():
    # testing against R CircStats package
    # Ref[1], page 180
    alpha = np.deg2rad(np.array([356, 97, 211, 232, 343, 292, 157, 302, 335,
                                 302, 324, 85, 324, 340, 157, 238, 254, 146,
                                 232, 122, 329]))
    beta = np.deg2rad(np.array([119, 162, 221, 259, 270, 29, 97, 292, 40, 313,
                                94, 45, 47, 108, 221, 270, 119, 248, 270, 45,
                                23]))
    answer = 0.2704648
    assert_allclose(answer, circcorrcoef(alpha, beta), atol=1e-4)


def test_rayleightest():
    # testing against R CircStats package
    data = np.array([3.3192729, 3.0628233, 2.7219562, 3.8019564, 2.7290561,
                     0.8840816, 3.7646373, 2.9782072, 3.3060406, 2.9012014])
    # answer was obtained through R CircStats function r.test(x)
    answer = (0.0009917748, 0.774055)
    result = (rayleightest(data), _length(data))
    assert_allclose(answer[0], result[0], atol=1e-9)
    assert_allclose(answer[1], result[1], atol=1e-4)

    data = np.array([3.8955614, 3.1700932, 3.1804325, 2.7887885, 2.9513829,
                     3.1140079, 4.5052974, 3.7984484, 3.7510773, 2.9764646,
                     2.2100108, 2.9529469, 2.9722527, 3.5099678, 2.6844710,
                     4.5020762, 2.9285673, 2.2871037, 3.5805262, 3.4247286,
                     2.4355773, 3.7549614, 4.1346585, 2.8426607, 3.3291801,
                     2.5847791, 2.7217270, 3.4323084, 3.3058256, 3.2147845,
                     2.4918005, 3.4849414, 0.8249985, 3.4881397, 3.2070389,
                     3.0859854, 4.0566486, 2.3463984, 3.7370984, 5.6853310,
                     5.8108009, 3.3921987, 3.2166975, 3.6827617, 4.5055291,
                     3.5295258, 3.0162183, 3.2317242, 3.2778354, 3.1713455])
    # answer was obtained through R CircStats function r.test(x)
    answer = (2.726928e-13, 0.7606633)
    result = (rayleightest(data), _length(data))
    assert_allclose(answer[0], result[0], atol=1e-18)
    assert_allclose(answer[1], result[1], atol=1e-5)


@pytest.mark.skipif('not HAS_SCIPY')
def test_vtest():
    # testing against R CircStats package
    data = np.array([1.316075, 4.439193, 3.096231, 4.807068, 2.986021,
                     1.756324, 3.046718, 3.299150, 3.360557, 4.842499])
    # answer was obtained through R CircStats function v0.test(x)
    answer = 0.987142
    assert_allclose(answer, vtest(data), atol=1e-5)


def test_vonmissesmle():
    # testing against R CircStats package
    data = np.array([3.3699057, 4.0411630, 0.5014477, 2.6223103, 3.7336524,
                     1.8136389, 4.1566039, 2.7806317, 2.4672173, 2.8493644])
    # answer was obtained through R CircStats function vm.ml(x)
    answer = (3.006514, 1.474132)
    result = vonmisesmle(data)
    assert_allclose(answer, result, atol=1e-5)
