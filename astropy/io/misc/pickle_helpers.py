# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains simple input/output related functionality that is not
part of a larger framework or standard.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six
from ...extern.six.moves import range


__all__ = ['fnpickle', 'fnunpickle']


def fnunpickle(fileorname, number=0, usecPickle=True):
    """ Unpickle pickled objects from a specified file and return the contents.

    Parameters
    ----------
    fileorname : str or file-like
        The file name or file from which to unpickle objects. If a file object,
        it should have been opened in binary mode.
    number : int
        If 0, a single object will be returned (the first in the file). If >0,
        this specifies the number of objects to be unpickled, and a list will
        be returned with exactly that many objects. If <0, all objects in the
        file will be unpickled and returned as a list.
    usecPickle : bool
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster). This only applies for python 2.x.

    Raises
    ------
    EOFError
        If ``number`` is >0 and there are fewer than ``number`` objects in the
        pickled file.

    Returns
    -------
    contents : obj or list
        If ``number`` is 0, this is a individual object - the first one unpickled
        from the file. Otherwise, it is a list of objects unpickled from the
        file.

    """

    if usecPickle and six.PY2:
        import cPickle as pickle
    else:
        import pickle

    if isinstance(fileorname, six.string_types):
        f = open(fileorname, 'rb')
        close = True
    else:
        f = fileorname
        close = False

    try:
        if number > 0:  # get that number
            res = []
            for i in range(number):
                res.append(pickle.load(f))
        elif number < 0:  # get all objects
            res = []
            eof = False
            while not eof:
                try:
                    res.append(pickle.load(f))
                except EOFError:
                    eof = True
        else:  # number==0
            res = pickle.load(f)
    finally:
        if close:
            f.close()

    return res


def fnpickle(object, fileorname, usecPickle=True, protocol=None, append=False):
    """Pickle an object to a specified file.

    Parameters
    ----------
    object
        The python object to pickle.
    fileorname : str or file-like
        The filename or file into which the `object` should be pickled. If a
        file object, it should have been opened in binary mode.
    usecPickle : bool
        If True (default), the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster). This only applies for python 2.x.
    protocol : int or None
        Pickle protocol to use - see the :mod:`pickle` module for details on
        these options. If None, the most recent protocol will be used.
    append : bool
        If True, the object is appended to the end of the file, otherwise the
        file will be overwritten (if a file object is given instead of a
        file name, this has no effect).

    """

    if usecPickle and six.PY2:
        import cPickle as pickle
    else:
        import pickle

    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL

    if isinstance(fileorname, six.string_types):
        f = open(fileorname, 'ab' if append else 'wb')
        close = True
    else:
        f = fileorname
        close = False

    try:
        pickle.dump(object, f, protocol=protocol)
    finally:
        if close:
            f.close()
