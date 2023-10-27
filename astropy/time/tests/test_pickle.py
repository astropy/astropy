# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pickle

import numpy as np

from astropy.time import Time


class TestPickle:
    """Basic pickle test of time"""

    def test_pickle(self):
        times = ["1999-01-01 00:00:00.123456789", "2010-01-01 00:00:00"]
        t1 = Time(times, scale="utc")

        for prot in range(pickle.HIGHEST_PROTOCOL):
            t1d = pickle.dumps(t1, prot)
            t1l = pickle.loads(t1d)
            assert np.all(t1l == t1)

        t2 = Time("2012-06-30 12:00:00", scale="utc")

        for prot in range(pickle.HIGHEST_PROTOCOL):
            t2d = pickle.dumps(t2, prot)
            t2l = pickle.loads(t2d)
            assert t2l == t2

    def test_cache_not_shared(self):
        t = Time(["2001:020", "2001:040", "2001:060", "2001:080"], out_subfmt="date")
        # Ensure something is in the cache.
        t.value
        assert "format" in t.cache
        td = pickle.dumps(t)
        assert "format" in t.cache
        tl = pickle.loads(td)
        assert "format" in t.cache
        assert "format" not in tl.cache
        t[0] = "1999:099"
        assert t.value[0] == "1999:099"
        assert tl.value[0] == "2001:020"
