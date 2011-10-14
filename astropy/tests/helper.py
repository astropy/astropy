import sys
import base64
import zlib
import imp

try:
    import pytest
    HAVE_PYTEST = True
except ImportError:
    HAVE_PYTEST = False

from .. import __path__ as astropy_path
from ..extern import pytest as extern_pytest

def pytest_main(args=None,plugins=None):
    """
    Implements pytest.main() after importing from a stand-alone module produced
    with py.test --genscript. Method adapted from file created by
    py.test --genscript.

    See Also
    --------
    pytest.main : This takes the same arguments.

    """
    if sys.version_info >= (3, 0):
        exec("def do_exec_def(co, loc): exec(co, loc)\n")
        extern_pytest.do_exec = do_exec_def

        import pickle
        unpacked_sources = extern_pytest.sources.encode("ascii") # ensure bytes
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodebytes(unpacked_sources)))
    else:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        extern_pytest.do_exec = do_exec_def

        import cPickle as pickle
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodestring(extern_pytest.sources)))

    importer = extern_pytest.DictImporter(unpacked_sources)
    sys.meta_path.append(importer)

    pytest = importer.load_module('pytest')
    pytest.main(args=args,plugins=plugins)


def run_tests(module=None, args=None, plugins=None, verbose=False, pastebin=None):

    import os.path

    if HAVE_PYTEST:
        main = pytest.main
    else:
        main = pytest_main

    if module is None:
        module_path = astropy_path[0]
    else:
        module_path = os.path.join(astropy_path[0],module.replace('.',os.path.sep))

        if not os.path.isdir(module_path):
            raise ValueError('Module not found: {0}'.format(module))

    all_args = module_path
    if args is not None:
        all_args += " {}".format(args)
    if verbose:
        all_args += " -v"
    if pastebin is not None:
        if pastebin in ['failed', 'all']:
            all_args += " --pastebin={}".format(pastebin)
        else:
            raise Exception("pastebin should be 'failed' or 'all'")

    main(args=all_args, plugins=plugins)
