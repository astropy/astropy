from distutils.core import setup, Extension
from os.path import join
import sys

from astropy import setup_helpers


def get_extensions(build_type='release'):
    EXPAT_DIR = 'cextern/expat/lib'
    VO_DIR = 'astropy/io/vo/src'

    defines = [("HAVE_EXPAT_CONFIG_H", 1)]
    if sys.byteorder == 'big':
        defines.append(('BYTEORDER', '4321'))
    else:
        defines.append(('BYTEORDER', '1234'))
    if sys.platform != 'win32':
        defines.append(('HAVE_UNISTD_H', None))

    return [Extension(
        "astropy.io.vo.iterparser",
        [join(VO_DIR, "iterparse.c"),
         join(EXPAT_DIR, "xmlparse.c"),
         join(EXPAT_DIR, "xmlrole.c"),
         join(EXPAT_DIR, "xmltok.c"),
         join(EXPAT_DIR, "xmltok_impl.c")],
        define_macros=defines,
        include_dirs=[VO_DIR, EXPAT_DIR])]


def get_package_data():
    return {
        'astropy.io.vo': [
            'data/ucd1p-words.txt', 'data/*.xsd', 'data/*.dtd'],
        'astropy.io.vo.tests': [
            'data/*.xml', 'data/*.gz', 'data/*.json', 'data/*.fits']}


def get_legacy_alias():
    return setup_helpers.add_legacy_alias('vo', 'astropy.io.vo')
