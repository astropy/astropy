import sys
import base64
import zlib
import imp

try:
    import pytest
    HAVE_PYTEST = True
except ImportError:
    HAVE_PYTEST = False
    
import astropy
import astropy.extern.pytest

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
        astropy.extern.pytest.do_exec = do_exec_def
        
        import pickle
        unpacked_sources = astropy.extern.pytest.sources.encode("ascii") # ensure bytes
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodebytes(unpacked_sources)))
    else:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        astropy.extern.pytest.do_exec = do_exec_def

        import cPickle as pickle
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodestring(astropy.extern.pytest.sources)))

    importer = astropy.extern.pytest.DictImporter(unpacked_sources)
    sys.meta_path.append(importer)

    pytest = importer.load_module('pytest')
    pytest.main(args=args,plugins=plugins)
    

def runtests(module=None):
    import os.path
    
    if HAVE_PYTEST:
        main = pytest.main
    else:
        main = pytest_main
    
    if module is None:
        main(astropy.__path__[0])
    else:
        module_path = os.path.join(astropy.__path__[0],module.replace('.',os.path.sep))
        
        if not os.path.isdir(module_path):
            raise ValueError('Module not found: {0}'.format(module))
            
        main(args=module)
