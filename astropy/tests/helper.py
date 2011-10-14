"""
Implements pytest.main() after importing from a stand-alone module produced with
py.test --genscript. Method adapted from file created by py.test --genscript.

"""

import sys
import base64
import zlib
import imp

import astropy.extern.pytest as pytest

def main(args=None,plugins=None):
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
        pytest.do_exec = do_exec_def
        
        import pickle
        unpacked_sources = pytest.sources.encode("ascii") # ensure bytes
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodebytes(unpacked_sources)))
    else:
        exec("def do_exec_def(co, loc): exec co in loc\n")
        pytest.do_exec = do_exec_def

        import cPickle as pickle
        unpacked_sources = pickle.loads(zlib.decompress(base64.decodestring(pytest.sources)))

    importer = pytest.DictImporter(unpacked_sources)
    sys.meta_path.append(importer)

    py = importer.load_module('pytest')
    py.main(args=args,plugins=plugins)
