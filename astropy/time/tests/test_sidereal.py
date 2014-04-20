# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import itertools

import numpy as np

from ...tests.helper import pytest
from .. import Time
from ..core import SIDEREAL_TIME_MODELS

allclose_hours = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)
# 1 nanosec atol
within_1_second = functools.partial(np.allclose, rtol=1., atol=1./3600.)
within_2_seconds = functools.partial(np.allclose, rtol=1., atol=2./3600.)


def test_doc_string_contains_models():
    """The doc string is formatted; this ensures this remains working."""
    for kind in ('mean', 'apparent'):
        for model in SIDEREAL_TIME_MODELS[kind]:
            assert model in Time.sidereal_time.__doc__


class TestERFATestCases():
    """Test that we reproduce the test cases given in erfa/src/t_erfa_c.c"""
    # all tests use the following JD inputs
    time_ut1 = Time(2400000.5, 53736.0, scale='ut1', format='jd')
    time_tt = Time(2400000.5, 53736.0, scale='tt', format='jd')
    # but tt!=ut1 at these dates, unlike what is assumed, so we cannot
    # reproduce this exactly. Now it does not really matter,
    # but may as well fake this (and avoid IERS table lookup here)
    time_ut1.delta_ut1_utc = 0.
    time_ut1.delta_ut1_utc = 24*3600*(time_ut1.tt.jd2-time_tt.jd2)
    assert np.allclose(time_ut1.tt.jd2 - time_tt.jd2, 0., atol=1.e-14)

    @pytest.mark.parametrize('erfa_test_input',
                             ((1.754174972210740592, 1e-12, "eraGmst00"),
                              (1.754174971870091203, 1e-12, "eraGmst06"),
                              (1.754174981860675096, 1e-12, "eraGmst82"),
                              (1.754166138018281369, 1e-12, "eraGst00a"),
                              (1.754166136510680589, 1e-12, "eraGst00b"),
                              (1.754166137675019159, 1e-12, "eraGst06a"),
                              (1.754166136020645203, 1e-12, "eraGst94")))
    def test_iau_models(self, erfa_test_input):
        result, precision, name = erfa_test_input
        if name[4] == 'm':
            kind = 'mean'
            model_name = 'IAU{0:2d}{1:s}'.format(20 if name[7] == '0' else 19,
                                                 name[7:])
        else:
            kind = 'apparent'
            model_name = 'IAU{0:2d}{1:s}'.format(20 if name[6] == '0' else 19,
                                                 name[6:].upper())

        assert kind in SIDEREAL_TIME_MODELS.keys()
        assert model_name in SIDEREAL_TIME_MODELS[kind]

        model = SIDEREAL_TIME_MODELS[kind][model_name]
        gst_erfa = self.time_ut1._erfa_sidereal_time(model)
        assert np.allclose(gst_erfa.to('radian').value, result,
                           rtol=1., atol=precision)

        gst = self.time_ut1.sidereal_time(kind, 'greenwich', model_name)
        assert np.allclose(gst.to('radian').value, result,
                           rtol=1., atol=precision)


class TestST():
    """Test Greenwich Sidereal Time.  Unlike above, this is relative to
    what was found earlier, so checks changes in implementation, including
    leap seconds, rather than correctness"""

    t1 = Time(['2012-06-30 12:00:00', '2012-06-30 23:59:59',
               '2012-06-30 23:59:60', '2012-07-01 00:00:00',
               '2012-07-01 12:00:00'], scale='utc')
    t2 = Time(t1, location=('120d', '10d'))

    def test_gmst(self):
        """Compare Greenwich Mean Sidereal Time with what was found earlier
        """
        gmst_compare = np.array([6.5968497894730564, 18.629426164144697,
                                 18.629704702452862, 18.629983240761003,
                                 6.6628381828899643])
        gmst = self.t1.sidereal_time('mean', 'greenwich')
        assert allclose_hours(gmst.value, gmst_compare)

    def test_gst(self):
        """Compare Greenwich Apparent Sidereal Time with what was found earlier
        """
        gst_compare = np.array([6.5971168570494854, 18.629694220878296,
                                18.62997275921186, 18.630251297545389,
                                6.6631074284018244])
        gst = self.t1.sidereal_time('apparent', 'greenwich')
        assert allclose_hours(gst.value, gst_compare)

    def test_gmst_gst_close(self):
        """Check that Mean and Apparent are within a few seconds."""
        gmst = self.t1.sidereal_time('mean', 'greenwich')
        gst = self.t1.sidereal_time('apparent', 'greenwich')
        assert within_2_seconds(gst.value, gmst.value)

    def test_gmst_independent_of_self_location(self):
        """Check that Greenwich time does not depend on self.location"""
        gmst1 = self.t1.sidereal_time('mean', 'greenwich')
        gmst2 = self.t2.sidereal_time('mean', 'greenwich')
        assert allclose_hours(gmst1.value, gmst2.value)

    @pytest.mark.parametrize('kind', ('mean', 'apparent'))
    def test_lst(self, kind):
        """Compare Local Sidereal Time with what was found earlier,
        as well as with what is expected from GMST
        """
        lst_compare = {
            'mean': np.array([14.596849789473058, 2.629426164144693,
                              2.6297047024528588, 2.6299832407610033,
                              14.662838182889967]),
            'apparent': np.array([14.597116857049487, 2.6296942208782959,
                                  2.6299727592118565, 2.6302512975453887,
                                  14.663107428401826])}

        gmst2 = self.t2.sidereal_time(kind, 'greenwich')
        lmst2 = self.t2.sidereal_time(kind)
        assert allclose_hours(lmst2.value, lst_compare[kind])
        assert allclose_hours((lmst2-gmst2).wrap_at('12h').value,
                              self.t2.location.longitude.to('hourangle').value)
        # check it also works when one gives longitude explicitly
        lmst1 = self.t1.sidereal_time(kind, self.t2.location.longitude)
        assert allclose_hours(lmst1.value, lst_compare[kind])

    def test_lst_needs_location(self):
        with pytest.raises(ValueError):
            self.t1.sidereal_time('mean')
        with pytest.raises(ValueError):
            self.t1.sidereal_time('mean', None)


class TestModelInterpretation():
    """Check that models are different, and that wrong models are recognized"""
    t = Time(['2012-06-30 12:00:00'], scale='utc', location=('120d', '10d'))

    @pytest.mark.parametrize('kind', ('mean', 'apparent'))
    def test_model_uniqueness(self, kind):
        """Check models give different answers, yet are close."""
        for model1, model2 in itertools.combinations(
                SIDEREAL_TIME_MODELS[kind].keys(), 2):
            gst1 = self.t.sidereal_time(kind, 'greenwich', model1)
            gst2 = self.t.sidereal_time(kind, 'greenwich', model2)
            assert np.all(gst1.value != gst2.value)
            assert within_1_second(gst1.value, gst2.value)
            lst1 = self.t.sidereal_time(kind, None, model1)
            lst2 = self.t.sidereal_time(kind, None, model2)
            assert np.all(lst1.value != lst2.value)
            assert within_1_second(lst1.value, lst2.value)

    @pytest.mark.parametrize(('kind', 'other'), (('mean', 'apparent'),
                                                 ('apparent', 'mean')))
    def test_wrong_models_raise_exceptions(self, kind, other):

        with pytest.raises(ValueError):
            self.t.sidereal_time(kind, 'greenwich', 'nonsense')

        for model in (set(SIDEREAL_TIME_MODELS[other].keys()) -
                      set(SIDEREAL_TIME_MODELS[kind].keys())):
            with pytest.raises(ValueError):
                self.t.sidereal_time(kind, 'greenwich', model)
            with pytest.raises(ValueError):
                self.t.sidereal_time(kind, None, model)
