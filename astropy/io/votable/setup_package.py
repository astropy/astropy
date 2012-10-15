from distutils.core import Extension
from os.path import join

from astropy import setup_helpers


def get_extensions(build_type='release'):
    VO_DIR = 'astropy/io/votable/src'

    return [Extension(
        "astropy.io.votable.tablewriter",
        [join(VO_DIR, "tablewriter.c")],
        include_dirs=[VO_DIR])]


def get_package_data():
    return {
        'astropy.io.votable': [
            'data/ucd1p-words.txt', 'data/*.xsd', 'data/*.dtd'],
        'astropy.io.votable.tests': [
            'data/*.xml', 'data/*.gz', 'data/*.json', 'data/*.fits',
            'data/*.txt'],
        'astropy.io.votable.validator': [
            'urls/*.dat.gz']}


def get_legacy_alias():
    return setup_helpers.add_legacy_alias(
        'vo', 'astropy.io.votable', '0.8')
