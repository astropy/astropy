import sys
import base64
import zlib
import functools

from .. import __path__ as astropy_path


try:
    import pytest

except ImportError:
    from ..extern import pytest as extern_pytest

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
    

# pytest marker to mark tests which get data from the web
big_data = pytest.mark.big_data

# these pytest hooks allow us to mark tests and run the marked tests with
# specific command line options.
def pytest_addoption(parser):
    parser.addoption("--runbigdata", action="store_true",
        help="run tests with online data")

def pytest_runtest_setup(item):
    if 'big_data' in item.keywords and not item.config.getvalue("runbigdata"):
        pytest.skip("need --runbigdata option to run")
        
        
def run_tests(module=None, args=None, plugins=None, verbose=False, pastebin=None):
    """
    Run Astropy tests using py.test. A proper set of arguments is constructed
    and passed to `pytest.main`.

    Parameters
    ----------
    module : str, optional
        The name of a specific module to test, e.g. 'io.fits' or 'utils'.
        If nothing is specified all default Astropy tests are run.

    args : str, optional
        Additional arguments to be passed to `pytest.main` in the `args`
        keyword argument.

    plugins : str, optional
        Arguments to passed to `pytest.main` in the `plugins` keyword argument.

    verbose : bool, optional
        Convenience option to turn on verbose output from py.test. Passing True
        is the same as specifying `-v` in `args`.

    pastebin : {'failed','all',True,None}, optional
>>>>>>> Updated run_tests docstring with better explanation of the pastebin keyword.

    See Also
    --------
    pytest.main : py.test function wrapped by `run_tests`.

    """
    import os.path

    if module is None:
        module_path = astropy_path[0]
    else:
        module_path = os.path.join(astropy_path[0],module.replace('.',os.path.sep))

        if not os.path.isdir(module_path):
            raise ValueError('Module not found: {0}'.format(module))

    all_args = module_path
    if args is not None:
        all_args += " {0}".format(args)
    if verbose:
        all_args += " -v"
    if pastebin is not None:
        if pastebin in ['failed', 'all']:
            all_args += " --pastebin={0}".format(pastebin)
        elif pastebin is True:
            all_args += " --pastebin=failed"
        else:
            raise ValueError("pastebin should be 'failed' or 'all'")

    pytest.main(args=all_args, plugins=plugins)


class raises:
    """
    A decorator to mark that a test should raise a given exception.
    Use as follows::

        @raises(ZeroDivisionError)
        def test_foo():
            x = 1/0
    """
    # pep-8 naming exception -- this is a decorator class
    def __init__(self, exc):
        self._exc = exc

    def __call__(self, func):
        @functools.wraps(func)
        def run_raises_test():
            pytest.raises(self._exc, func)
        return run_raises_test
