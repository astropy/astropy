# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from .. import Time


class TestGuess():
    """Test guessing the input value format"""

    def test_guess1(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
        t = Time(times, scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='iso' "
                "value=['1999-01-01 00:00:00.123' '2010-01-01 00:00:00.000']>")

    def test_guess2(self):
        times = ['1999-01-01 00:00:00.123456789', '2010-01 00:00:00']
        with pytest.raises(ValueError):
            Time(times, scale='utc')

    def test_guess3(self):
        times = ['1999:001:00:00:00.123456789', '2010:001']
        t = Time(times, scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='yday' "
                "value=['1999:001:00:00:00.123' '2010:001:00:00:00.000']>")

    def test_guess4(self):
        times = [10, 20]
        with pytest.raises(ValueError):
            Time(times, scale='utc')

    def test_guess5(self):
        # Fractional days ISO
        times = ['1999-01-01.123', '2010-01-01.']
        t = Time(times, scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='iso' "
                "value=['1999-01-01 02:57:07.200' '2010-01-01 00:00:00.000']>")

    def test_guess6(self):
        # Fractional days yday
        times = ['1999:123.123', '2010:1.0']
        t = Time(times, scale='utc')
        assert (repr(t) == "<Time object: scale='utc' format='yday' "
                "value=['1999:123:02:57:07.200' '2010:001:00:00:00.000']>")
