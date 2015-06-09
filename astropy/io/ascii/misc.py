"""A Collection of useful miscellaneous functions.

misc.py:
  Collection of useful miscellaneous functions.

:Copyright: None -- Hack away at your leisure!!
:Author: Hannes Breytenbach (hannes@saao.ac.za)
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import itertools as itt
import collections as coll

from ...extern.six.moves import zip, map, filter


def nth(iterable, n, default=None):
    '''Returns the nth item or a default value'''
    return next(itt.islice(iterable, n, None), default)

def nthzip(n, *its):
    '''Return the nth component of the zipped sequence'''
    return tuple( nth(it, n) for it in its )
    
def accordion(it):
    ''' interleave an iterator with a -1 rolled version of itself'''
    return list( zip(it, np.roll(it,-1)) )[:-1]
    
def first_true_idx(iterable, pred=None, default=None):
    '''find the first index position for the which the callable pred returns True'''
    if pred is None:    func = lambda x : x[1]
    else:               func = lambda x : pred(x[1])
    ii = next( filter(func, enumerate(iterable)), default )     #either index-item pair or default
    return ii[0] if ii else default

def first_false_idx(iterable, pred=None, default=None):
    '''find the first index position for the which the callable pred returns False'''
    if pred is None:    func = lambda x : not x
    else:               func = lambda x : not pred(x)
    return first_true_idx(iterable, func, default)
    
def where_true( iterable, pred=None ):
    '''Return the indeces of an iterable for which the callable pred evaluates as True'''
    func = lambda x: pred(x[1])
    return nthzip(0, *filter(func, enumerate(iterable)))
    
    
def sortmore(*args, **kw):
    '''
    Sorts any number of lists according to: 
    optionally given item sorting key function(s) and/or a global sorting key function.
    
    Parameters
    ----------
    globalkey : None                revert to sorting by key function
    globalkey : callable            Sort by evaluated value for all items in the lists
                                    (call signiture of this function needs to be such that it accepts an
                                    argument tuple of items from each list.
                                    eg.: globalkey = lambda *l : sum(l) will order all the lists by the 
                                    sum of the items from each list
    
    if key : None                   sorting done by value of first input list 
                                    (in this case the objects in the first iterable need the comparison 
                                    methods __lt__ etc...)
    if key : callable               sorting done by value of key(item) for items in first iterable
    if key : tuple                  sorting done by value of (key(item_0), ...,key(item_n)) for items in 
                                    the first n iterables (where n is the length of the key tuple)
                                    i.e. the first callable is the primary sorting criterion, and the 
                                    rest act as tie-breakers.
    '''
    
    farg = list( args[0] )
    if not len(farg):
        return args
    
    globalkey   =       kw.get('globalkey')
    key         =       kw.get('key')
    if key is None:
        if globalkey:   key = lambda x: 0               #if global sort function given and no local (secondary) key given, ==> no tiebreakers
        else:           key = lambda x: x               #if no global sort and no local sort keys given, sort by item values
    if globalkey is None:         
        globalkey = lambda *x: 0
    
    if not isinstance(globalkey, coll.Callable):
        raise ValueError( 'globalkey needs to be callable' )
        
    if isinstance(key, coll.Callable):
        k = lambda x: (globalkey(*x), key(x[0]))
    elif isinstance(key, tuple):
        key = ( k if k else lambda x: 0 for k in key )
        k = lambda x : (globalkey(*x),) + tuple( f(z) for (f,z) in zip(key,x) )
    else:
        raise KeyError( "kw arg 'key' should be None, callable, or a sequence of callables, not {}".format(type(key)) )
    
    res = sorted( list(zip(*args)), key=k )
    if 'order' in kw:
        if kw['order'].startswith(('descend', 'reverse')):
            res = reversed(res)
    
    return tuple( map(list, zip(*res) ) )
    
    
#====================================================================================================
def groupmore(func=None, *its):
    if not func:        func = lambda x: x
    its = sortmore(*its, key=func)
    nfunc = lambda x : func(x[0])
    zipper = itt.groupby( zip(*its), nfunc )
    unzipper = ((key, zip(*groups)) for key, groups in zipper)
    return unzipper