from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


from ...utils.data import get_pkg_data_filename

fd = None


def test_open_file_detection():
    global fd
    fd = open(get_pkg_data_filename('data/open_file_detection.txt'))


def teardown():
    if fd is not None:
        fd.close()
