from astropy.utils.data import get_pkg_data_filename

def get_data_filename():
    return get_pkg_data_filename('data/foo.txt')
