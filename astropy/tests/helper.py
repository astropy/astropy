import sys
import base64
import zlib

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
    

def run_tests(module=None, args=None, plugins=None, verbose=False, pastebin=None):
    """
    Run AstroPy tests using py.test. A proper set of arguments is constructed
    and passed to `pytest.main`.
    
    Parameters
    ----------
    module : str, optional
        The name of a specific module to test, e.g. 'io.fits' or 'utils'.
        If nothing is specified all default AstroPy tests are run.
        
    args : str, optional
        Additional arguments to be passed to `pytest.main` in the `args`
        keyword argument. 
        
    plugins : str, optional
        Arguments to passed to `pytest.main` in the `plugins` keyword argument.
        
    verbose : bool, optional
        Convenience option to turn on verbose output from py.test. Passing True
        is the same as specifying `-v` in `args`.
        
    pastebin : {'failed','all',True,None}, optional
        Convenience option for turning on py.test pastebin output.
        
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
        elif pastebin:
            all_args += " --pastebin=failed"
        else:
            raise ValueError("pastebin should be 'failed' or 'all'")

    pytest.main(args=all_args, plugins=plugins)
