import sys

import pytest

from ...utils.data import get_pkg_data_filename

fd = None
@pytest.mark.skipif("sys.version_info >= (3, 0)")
def test_open_file_detection():
    global fd
    fd = open(get_pkg_data_filename('data/open_file_detection.txt'))


def teardown():
    if fd is not None:
        fd.close()
