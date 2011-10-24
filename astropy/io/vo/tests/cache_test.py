from __future__ import absolute_import

# STDLIB
import BaseHTTPServer
import datetime
import multiprocessing
import os
import shutil
import SimpleHTTPServer
import socket
import tempfile
import time
import urllib2

# THIRD PARTY
from numpy.testing import assert_raises
from numpy.testing.decorators import knownfailureif
import numpy as np

# LOCAL
from astropy.io.vo import cache, conesearch, image, ssa, vos_catalog
from astropy.io.vo.util import IS_PY3K

"""
To test this stuff, a local webserver is started up in another process
containing dummy content that can be queried.  A temporary directory
is used in place of the local cache so that running tests will not
interfere with real cached content.

The port is selected until an available one is found.
"""

# TODO: Some tests are disabled because they rely upon "real" external
# servers which go down from time to time and probably don't
# appreciate the spamming either.

class SuperSimpleHTTPRequestHandler(SimpleHTTPServer.SimpleHTTPRequestHandler):
    def do_GET(self, *args, **kwargs):
        question_mark = self.path.find('?')
        if question_mark >= 0:
            self.path = self.path[:question_mark]
        return SimpleHTTPServer.SimpleHTTPRequestHandler.do_GET(self, *args, **kwargs)


def run_server(port):
    # SimpleHTTPServer serves a directory of files in the current
    # directory, so we chdir to the data directory here.
    data_dir = os.path.join(os.path.dirname(__file__) or '.', 'data')
    os.chdir(data_dir)
    httpd = BaseHTTPServer.HTTPServer(
        ('', port), SuperSimpleHTTPRequestHandler)
    httpd.serve_forever()


class TestCache:
    def setup_class(self):
        port = None
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        for i in range(8888, 10000, 107):
            try:
                s.bind(('', i))
                port = i
            except:
                pass
            else:
                break
        s.close()

        if port is None:
            raise RuntimeError(
                "Could not find open port for dummy HTTP server.")

        self.port = port
        self.server_url = 'http://127.0.0.1:%d/' % port
        self.server = multiprocessing.Process(target=run_server, args=(port,))
        self.server.start()

        # Give some time for the server to start before hammering
        # requests at it
        time.sleep(1)

        # Use a tmp directory for the cache dir
        self.cache_dir = tempfile.mkdtemp()
        self.cache = cache.Cache(
            self.cache_dir, self.server_url,
            cache_delta=datetime.timedelta(seconds=1))

    def teardown_class(self):
        # Kill the server process
        self.server.terminate()
        time.sleep(1)

        # This is actually an additional test -- now that we've killed
        # the server, check that getting a cached file still works,
        # even though the server itself would throw a 404.
        x = vos_catalog.get_remote_catalog_db('basic', self.cache)
        assert x._tree['content'] == ['A', 'B', 'C']

        # Mark the file as really old, so an HTTP attempt will be
        # requested and fail.
        cache_file = self.cache.get_cache_file_path('basic.json')
        os.utime(cache_file, (0, 0))
        x = vos_catalog.get_remote_catalog_db('basic', self.cache)
        assert x._tree['content'] == ['A', 'B', 'C']

        # Delete the temporary cache dir
        shutil.rmtree(self.cache_dir)

    def test_load(self):
        x = vos_catalog.get_remote_catalog_db('basic', self.cache)
        assert x._tree['content'] == ['A', 'B', 'C']
        assert self.cache.has_cached_copy('basic.json')
        x = vos_catalog.get_remote_catalog_db('basic', self.cache)
        assert len(self.cache.get_cache_backups('basic.json')) == 0

        # Let some time pass to invalidate the cache
        time.sleep(2)
        x = vos_catalog.get_remote_catalog_db('basic', self.cache)
        assert x._tree['content'] == ['A', 'B', 'C']
        assert self.cache.has_cached_copy('basic.json')

    def test_404_failure(self):
        def raises():
            x = vos_catalog.get_remote_catalog_db('foo', self.cache)
        assert_raises(urllib2.HTTPError, raises)

    # def test_conesearch(self):
    #     tab = conesearch.conesearch(
    #         catalog_db=vos_catalog.get_remote_catalog_db('conesearch', self.cache),
    #         ra=1.0, dec=-10.0, sr=0.1, verb=1)
    #     assert len(tab.array)

    def test_dummy_conesearch_query(self):
        tab = conesearch.conesearch(
            catalog_db=self.server_url + "vo_cone.xml?",
            ra=1.0, dec=-10.0, sr=0.1, verb=1, pedantic=False)
        assert len(tab.array)

    # def test_bad_conesearch(self):
    #     def raises():
    #         tab = conesearch.conesearch(
    #             catalog_db=vos_catalog.get_remote_catalog_db('conesearch', self.cache),
    #             ra=1.0, dec=91.0, sr=0.1, verb=1)
    #     assert_raises(IOError, raises)

    def test_conesearch_list_catalogs(self):
        catalog_db = vos_catalog.get_remote_catalog_db('conesearch', self.cache)
        catalogs = catalog_db.list_catalogs()
        assert set(catalogs) == set(['BROKEN', 'USNO NOMAD', 'USNO-B1', 'USNO-A2', 'USNO ACT'])

    # def test_ssa(self):
    #     tab = ssa.query_data(
    #         catalog_db=vos_catalog.get_remote_catalog_db('ssa', self.cache),
    #         pos=(0.0, 0.0), size=5.0)
    #     assert len(tab.array)

    def test_dummy_ssa_query(self):
        tab = ssa.query_data(
            catalog_db=self.server_url + "hst_ssa.xml?",
            pos=(0.0, 0.0), size=5.0)
        assert len(tab.array)

    @knownfailureif(IS_PY3K, "Requires pyfits which is not available on Python 3.x")
    def test_ssa_get_data(self):
        import pyfits
        a = np.array([(self.server_url + "tb.fits", 0)],
                     dtype=[(('AcRef', '_AcRef'), 'O'),
                            (('foo', 'bar'), 'i')])
        x = ssa.get_data(a[0], mimetype='fits')
        if hasattr(pyfits, 'NP_pyfits'):
            assert isinstance(x, pyfits.NP_pyfits.HDUList)
        else:
            assert isinstance(x, pyfits.core.HDUList)

    # def test_image_query(self):
    #     tab = image.query(
    #         "http://archive.noao.edu/nvo/sim/voquery.php?",
    #         pos=(217,35), size=0.08)
    #     assert len(tab.array)

    def test_dummy_image_query(self):
        tab = image.query(
            catalog_db=self.server_url + "image.xml?",
            pos=(217, 35), size=0.08)
        assert len(tab.array)

if __name__ == '__main__':
    TestCache.setup_class()
    db = TestCache()
    time.sleep(1)
    db.test_load()
    TestCache.teardown_class()
