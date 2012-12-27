from distutils.core import Extension
from os.path import dirname, join, relpath, exists

ASTROPY_UTILS_ROOT = dirname(__file__)

def get_extensions():
    _generate_cython_pyx()

    return [
        Extension('astropy.utils._compiler',
                  [relpath(join(ASTROPY_UTILS_ROOT, 'src', 'compiler.c'))]),
        Extension('astropy.utils._cython',
                  [relpath(join(ASTROPY_UTILS_ROOT, 'src', 'cython.pyx'))])
        ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.utils.tests': [
            'data/*.dat',
            'data/*.dat.gz',
            'data/*.dat.bz2',
            'data/*.txt']
        }


def _generate_cython_pyx():
    # Check whether cython.pyx needs to be regenerated
    cython_pyx = relpath(join(ASTROPY_UTILS_ROOT, 'src', 'cython.pyx'))
    regenerate = True

    try:
        import Cython
        current_cython_version = Cython.__version__
    except ImportError:
        current_cython_version = 'unknown'

    if exists(cython_pyx):
        try:
            execfile(cython_pyx, globals())
            if cython_version == current_cython_version:
                regenerate = False
        except:
            # If there are any errors in the cython_pyx or the cython_version
            # variable is undefined then ignore and regenerate the file
            pass

    if regenerate:
        with open(cython_pyx, 'w') as f:
            f.write('# Generated file; do not modify\n')
            f.write('cython_version = {0!r}\n'.format(current_cython_version))
