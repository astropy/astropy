from distutils.core import Extension
from os.path import join

from astropy import setup_helpers


def get_extensions(build_type='release'):
    VO_DIR = 'astropy/io/vo/src'

    return [Extension(
        "astropy.io.vo.tablewriter",
        [join(VO_DIR, "tablewriter.c")],
        include_dirs=[VO_DIR])]


def get_package_data():
    return {
        'astropy.io.vo': [
            'data/ucd1p-words.txt', 'data/*.xsd', 'data/*.dtd'],
        'astropy.io.vo.tests': [
            'data/*.xml', 'data/*.gz', 'data/*.json', 'data/*.fits',
            'data/*.txt'],
        'astropy.io.vo.validator': [
            'urls/*.dat.gz']}


def get_legacy_alias():
    return setup_helpers.add_legacy_alias(
        'vo', 'astropy.io.vo', '0.8')
