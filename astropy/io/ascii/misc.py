"""A Collection of useful miscellaneous functions.

misc.py:
  Collection of useful miscellaneous functions.

:Author: Hannes Breytenbach (hannes@saao.ac.za)
"""

import collections.abc
import itertools
import operator


def first_true_index(iterable, pred=None, default=None):
    """find the first index position for the which the callable pred returns True."""
    if pred is None:
        func = operator.itemgetter(1)
    else:
        func = lambda x: pred(x[1])
    # either index-item pair or default
    ii = next(filter(func, enumerate(iterable)), default)
    return ii[0] if ii else default


def first_false_index(iterable, pred=None, default=None):
    """find the first index position for the which the callable pred returns False."""
    if pred is None:
        func = operator.not_
    else:
        func = lambda x: not pred(x)
    return first_true_index(iterable, func, default)


def sortmore(*args, **kw):
    """
    Sorts any number of lists according to:
    optionally given item sorting key function(s) and/or a global sorting key function.

    Parameters
    ----------
    *args : list
        One or more lists to sort.
    globalkey : callable or None, optional
        If `None` (default), sort by the ``key`` function. If a callable, sort
        by its evaluated value for all items in the lists (its call signature
        needs to be such that it accepts an argument tuple of items from each
        list, e.g. ``globalkey = lambda *l: sum(l)`` will order all the lists
        by the sum of the items from each list).
    key : callable or tuple of callable or None, optional
        If `None` (default), sorting is done by value of the first input list
        (in this case the objects in the first iterable need the comparison
        methods ``__lt__`` etc.). If a callable, sorting is done by value of
        ``key(item)`` for items in the first iterable. If a tuple, sorting is
        done by value of ``(key(item_0), ..., key(item_n))`` for items in the
        first n iterables (where n is the length of the key tuple), i.e. the
        first callable is the primary sorting criterion, and the rest act as
        tie-breakers.

    Returns
    -------
    sorted_lists : tuple of list
        The sorted lists.

    Examples
    --------
    Capture sorting indices::

        l = list('CharacterS')
        In [1]: sortmore( l, range(len(l)) )
        Out[1]: (['C', 'S', 'a', 'a', 'c', 'e', 'h', 'r', 'r', 't'],
                 [0, 9, 2, 4, 5, 7, 1, 3, 8, 6])
        In [2]: sortmore( l, range(len(l)), key=str.lower )
        Out[2]: (['a', 'a', 'C', 'c', 'e', 'h', 'r', 'r', 'S', 't'],
                 [2, 4, 0, 5, 7, 1, 3, 8, 9, 6])
    """
    first = list(args[0])
    if not first:
        return args

    globalkey = kw.get("globalkey")
    key = kw.get("key")
    if key is None:
        if globalkey:
            # if global sort function given and no local (secondary) key given, ==> no tiebreakers
            key = lambda x: 0
        else:
            # if no global sort and no local sort keys given, sort by item values
            key = lambda x: x
    if globalkey is None:
        globalkey = lambda *x: 0

    if not isinstance(globalkey, collections.abc.Callable):
        raise ValueError("globalkey needs to be callable")

    if isinstance(key, collections.abc.Callable):
        k = lambda x: (globalkey(*x), key(x[0]))
    elif isinstance(key, tuple):
        key = (k or (lambda x: 0) for k in key)
        k = lambda x: (globalkey(*x),) + tuple(f(z) for (f, z) in zip(key, x))
    else:
        raise KeyError(
            "kw arg 'key' should be None, callable, or a sequence of callables, "
            f"not {type(key)}"
        )

    res = sorted(zip(*args), key=k)
    if "order" in kw:
        if kw["order"].startswith(("descend", "reverse")):
            res = reversed(res)

    return tuple(map(list, zip(*res)))


def groupmore(func=None, *its):
    """Extends the itertools.groupby functionality to arbitrary number of iterators."""
    if not func:
        func = lambda x: x
    its = sortmore(*its, key=func)
    nfunc = lambda x: func(x[0])
    zipper = itertools.groupby(zip(*its), nfunc)
    return ((key, zip(*groups)) for key, groups in zipper)
