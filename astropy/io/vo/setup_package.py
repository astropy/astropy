from distutils.core import setup, Extension
import sys


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
        "vo.iterparser",
        [VO_DIR + "/iterparse.c",
         EXPAT_DIR + "/xmlparse.c",
         EXPAT_DIR + "/xmlrole.c",
         EXPAT_DIR + "/xmltok.c",
         EXPAT_DIR + "/xmltok_impl.c"],
        define_macros=defines,
        include_dirs=[VO_DIR, EXPAT_DIR])]


def get_package_data():
    return {
        'astropy.io.vo': [
            'data/ucd1p-words.txt', 'data/*.xsd', 'data/*.dtd'],
        'astropy.io.vo.tests': [
            'data/*.xml', 'data/*.gz', 'data/*.json', 'data/*.fits']}
