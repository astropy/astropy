# -*- coding: utf-8 *-*
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

def do_import_test(pkg=None, skiptests=True, queue=None):
    """
    Recursively imports packages.

    While this is used here as part of the `test_all_imports` test, affiliated
    packages can call this function, and it will run for their packages.

    Parameters
    ----------
    pkg : str, module, or None
        The root name of the package to recursively import. If None, the
        root package from whatever function calls this will be used, and if a
        module, it will be taken as the root package
    skiptests : bool
        If True, won't import modules that begin with 'test'
    queue : multiprocessing.Queue or None
        The return value will be put in this queue if it is not None

    Returns
    -------
    failedmods : list of 2-tuples
        A list of modules that raise exceptions on import. Each entry is a
        tuple (module_name, exception_at_import, traceback_str).

    """
    from pkgutil import walk_packages
    from ...utils import find_current_module
    from sys import exc_info
    from traceback import format_exc
    from inspect import ismodule

    exc_list = []

    if pkg is None:
        pkg = find_current_module(2)
        pkgname = pkg.__name__
    elif isinstance(pkg, basestring):
        try:
            pkg = __import__(pkg)
        except Exception, e:
            return [(pkgname, e)]
    elif ismodule(pkg):
        pass
    else:
        raise TypeError('pkg is not None, module, or string')

    rootpkg = __import__(pkg.__name__.split('.')[0])
    rootpath = rootpkg.__path__
    rootprefix = rootpkg.__name__ + '.'

    #capture any errors during import of the packages
    pkgerr = lambda n: exc_list.append((n, exc_info()[1], format_exc()))

    for importer, modname, ispkg in walk_packages(rootpath, rootprefix, pkgerr):
        if not ispkg:  # packages are already imported by walk_packages
            if skiptests:
                unqual_name = modname.split('.')[-1]
                if unqual_name != 'test' and unqual_name.startswith('test'):
                    continue
            try:
                importer.find_module(modname).load_module(modname)
            except:
                exc_list.append((modname, exc_info()[1], format_exc()))

    if queue:
        queue.put(exc_list)
    return exc_list


#this is the list of modules that fail for some known good reason
known_import_fails = [
'astropy.utils.compat._gzip_32',  # fails if not py 3.x
#These all fail in py 3.x because of the compiler moule, which is only in 2.x
'astropy.sphinx.ext.comment_eater',
'astropy.sphinx.ext.compiler_unparse',
'astropy.sphinx.ext.phantom_import',
'astropy.sphinx.ext.traitsdoc',
#These fail if sphinx is not present
'astropy.sphinx.ext.docscrape_sphinx',
'astropy.sphinx.ext.numpydoc',
'astropy.sphinx.ext.automodsumm'
]


@pytest.mark.importtest
def test_all_imports():
    """
    Tries to import all astropy packages (that aren't tests) to make sure they
    will all import without errors.

    Note that this uses multiprocessing to manage independent processes so as
    not to screw up other tests that might be doing something odd.
    """
    from multiprocessing import Process, Queue
    from ...utils import find_current_module

    q = Queue()
    p1 = Process(target=do_import_test, args=('astropy', True, q))
    p1.start()
    p1.join()
    p2 = Process(target=do_import_test, args=(find_current_module(1), True, q))
    p2.start()
    p2.join()

    assert not q.empty()
    excres = q.get()
    assert not q.empty()
    excres2 = q.get()

    #the two do_import_test calls should do the same thing
    assert [nm for nm, e, tb in excres] == [nm for nm, e, tb in excres2]

    #filter out known import failures
    excres = [t for t in excres if t[0] not in known_import_fails]

    modnms = []
    for modnm, exc, tb in excres:
        modnms.append(modnm)
        msg = 'failed to import due to {0}: {1} '
        print 'Module', modnm, msg.format(exc.__class__.__name__, exc)
        print tb

    failedexcstr = ', '.join(modnms)
    msg = ('The following {0} modules did not import sucessfully: {1}'
           ' see stdout for tracebacks'.format(len(modnms), failedexcstr))
    assert not failedexcstr, msg
