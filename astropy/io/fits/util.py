import functools
import itertools
import os
import sys
import tempfile
import warnings

import numpy as np


BLOCK_SIZE = 2880 # the FITS block size


def itersubclasses(cls, _seen=None):
    """
    itersubclasses(cls)

    Generator over all subclasses of a given class, in depth first order.

    >>> list(itersubclasses(int)) == [bool]
    True
    >>> class A(object): pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL (new-style) classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)] #doctest: +ELLIPSIS
    ['type', ...'tuple', ...]

    From http://code.activestate.com/recipes/576949/
    """

    if not isinstance(cls, type):
        raise TypeError('itersubclasses must be called with '
                        'new-style classes, not %.100r' % cls)
    if _seen is None: _seen = set()
    try:
        subs = cls.__subclasses__()
    except TypeError: # fails only when cls is type
        subs = cls.__subclasses__(cls)
    for sub in subs:
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            for sub in itersubclasses(sub, _seen):
                yield sub


class lazyproperty(object):
    """
    Works similarly to property(), but computes the value only once.

    Adapted from the recipe at
    http://code.activestate.com/recipes/363602-lazy-property-evaluation
    """

    def __init__(self, fget, fset=None, fdel=None, doc=None):
        self._fget = fget
        self._fset = fset
        self._fdel = fdel
        if doc is None:
            self.__doc__ = fget.__doc__
        else:
            self.__doc__ = doc

    def __get__(self, obj, owner=None):
        if obj is None:
            return self
        key = self._fget.func_name
        if key not in obj.__dict__:
            val = self._fget(obj)
            obj.__dict__[key] = val
            return val
        else:
            return obj.__dict__[key]

    def __set__(self, obj, val):
        if self._fset:
            self._fset(obj, val)
        obj.__dict__[self._fget.func_name] = val

    def __delete__(self, obj):
        if self._fdel:
            self._fdel(obj)
        key = self._fget.func_name
        if key in obj.__dict__:
            del obj.__dict__[key]

    def getter(self, fget):
        return self.__ter(fget, 0)

    def setter(self, fset):
        return self.__ter(fset, 1)

    def deleter(self, fdel):
        return self.__ter(fdel, 2)

    def __ter(self, f, arg):
        args = [self._fget, self._fset, self._fdel, self.__doc__]
        args[arg] = f
        cls_ns = sys._getframe(1).f_locals
        for k, v in cls_ns.iteritems():
            if v is self:
                property_name = k
                break

        cls_ns[property_name] = lazyproperty(*args)

        return cls_ns[property_name]


def deprecated(message='', name='', alternative='', pending=False):
    """
    Used to mark a function as deprecated.

    TODO: Provide a class deprecation marker as well.

    To mark an attribute as deprecated, replace that attribute with a
    depcrecated property.

    Parameters
    ------------
    message : str, optional
        Override the default deprecation message.  The format specifier
        %(func)s may be used for the name of the function, and %(alternative)s
        may be used in the deprecation message to insert the name of an
        alternative to the deprecated function.

    name : str, optional
        The name of the deprecated function; if not provided the name is
        automatically determined from the passed in function, though this is
        useful in the case of renamed functions, where the new function is just
        assigned to the name of the deprecated function.  For example:
            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative function that the user may use in place of the
        deprecated function.  The deprecation warning will tell the user about
        this alternative if provided.

    pending : bool, optional
        If True, uses a PendingDeprecationWarning instead of a
        DeprecationWarning.

    """

    def deprecate(func):
        if isinstance(func, classmethod):
            try:
                func = func.__func__
            except AttributeError:
                # classmethods in Python2.6 and below lack the __func__
                # attribute so we need to hack around to get it
                method = func.__get__(None, object)
                if hasattr(method, '__func__'):
                    func = method.__func__
                elif hasattr(method, 'im_func'):
                    func = method.im_func
                else:
                    # Nothing we can do really...  just return the original
                    # classmethod
                    return func
            is_classmethod = True
        else:
            is_classmethod = False

        @functools.wraps(func)
        def deprecated_func(*args, **kwargs):
            # _message and _name are necessary; otherwise assignments to name
            # and message will cause the interpreter to treat them as local
            # variables, instead of encapuslated variables
            _message = message
            _name = name
            if not _name:
                _name = func.__name__

            if not _message or type(_message) == type(deprecate):
                if pending:
                    _message = 'The %(func)s function will be deprecated in' \
                               'a future version.'
                else:
                    _message = 'The %(func)s function is deprecated and may ' \
                               'be removed in a future version.'
                if alternative:
                    _message += '  Use %(alternative)s instead.'

            if pending:
                category = DeprecationPendingWarning
            else:
                category = DeprecationWarning

            warnings.warn(
                _message % {'func': _name, 'alternative': alternative},
                category, stacklevel=2)

            return func(*args, **kwargs)
        if is_classmethod:
            deprecated_func = classmethod(deprecated_func)
        return deprecated_func

    if type(message) == type(deprecate):
        return deprecate(message)

    return deprecate


def pairwise(iterable):
    """Return the items of an iterable paired with its next item.

    Ex: s -> (s0,s1), (s1,s2), (s2,s3), ....
    """

    a, b = itertools.tee(iterable)
    for _ in b:
        # Just a little trick to advance b without having to catch
        # StopIter if b happens to be empty
        break
    return itertools.izip(a, b)


def encode_ascii(s):
    """
    In Python 2 this is a no-op.  Strings are left alone.  In Python 3 this
    will be replaced with a function that actually encodes unicode strings to
    ASCII bytes.
    """

    return s


def decode_ascii(s):
    """
    In Python 2 this is a no-op.  Strings are left alone.  In Python 3 this
    will be replaced with a function that actually decodes ascii bytes to
    unicode.
    """

    return s


def isreadable(f):
    """
    Returns True if the file-like object can be read from.  This is a common-
    sense approximation of io.IOBase.readable.
    """

    if hasattr(f, 'closed') and f.closed:
        # This mimics the behavior of io.IOBase.readable
        raise ValueError('I/O operation on closed file')

    if not hasattr(f, 'read'):
        return False

    if hasattr(f, 'mode') and not any((c in f.mode for c in 'r+')):
        return False

    # Not closed, has a 'read()' method, and either has no known mode or a
    # readable mode--should be good enough to assume 'readable'
    return True


def iswritable(f):
    """
    Returns True if the file-like object can be written to.  This is a common-
    sense approximation of io.IOBase.writable.
    """

    if hasattr(f, 'closed') and f.closed:
        # This mimics the behavior of io.IOBase.writable
        raise ValueError('I/O operation on closed file')

    if not hasattr(f, 'write'):
        return False

    if hasattr(f, 'mode') and not any((c in f.mode for c in 'wa+')):
        return False

    # Note closed, has a 'write()' method, and either has no known mode or a
    # mode that supports writing--should be good enough to assume 'writable'
    return True


def isfile(f):
    """
    Returns True if the given object represents an OS-level file (that is,
    isinstance(f, file)).

    On Python 3 this also returns True if the given object is higher level
    wrapper on top of a FileIO object, such as a TextIOWrapper.
    """

    return isinstance(f, file)


def fileobj_name(f):
    """
    Returns the 'name' of file-like object f, if it has anything that could be
    called its name.  Otherwise f's class or type is returned.  If f is a
    string f itself is returned.
    """

    if isinstance(f, basestring):
        return f
    elif hasattr(f, 'name'):
        return f.name
    elif hasattr(f, 'filename'):
        return f.filename
    elif hasattr(f, '__class__'):
        return str(f.__class__)
    else:
        return str(type(f))


def fileobj_closed(f):
    """
    Returns True if the given file-like object is closed or if f is not a
    file-like object.
    """

    if hasattr(f, 'closed'):
        return f.closed
    elif hasattr(f, 'fileobj') and hasattr(f.fileobj, 'closed'):
        return f.fileobj.closed
    elif hasattr(f, 'fp') and hasattr(f.fp.closed):
        return f.fp.closed
    else:
        return True


def fileobj_mode(f):
    """
    Returns the 'mode' string of a file-like object if such a thing exists.
    Otherwise returns None.
    """

    # Go from most to least specific--for example gzip objects have a 'mode'
    # attribute, but it's not analogous to the file.mode attribute
    if hasattr(f, 'fileobj') and hasattr(f.fileobj, 'mode'):
        mode = f.fileobj.mode
    elif hasattr(f, 'fp') and hasattr(f.fp, 'mode'):
        mode = f.fp.mode
    elif hasattr(f, 'mode'):
        mode = f.mode
    else:
        mode = None

    if isinstance(mode, basestring):
        # On Python 3 files opened in 'a' mode actually get opened in 'w' or
        # 'r+' instead of 'a', since the behavior of 'a' mode is not normalized
        # across systems.  So we should normalize in the same way...
        if 'a' in mode:
            if '+' in mode:
                mode = mode.replace('a', 'r')
            else:
                mode = mode.replace('a', 'w')
    else:
        # If mode still turned out to be something other than a string, it's
        # not the droid we were looking for, and we should just toss it
        mode = None
    return mode


def translate(s, table, deletechars):
    """
    This is a version of string/unicode.translate() that can handle string or
    unicode strings the same way using a translation table made with
    string.maketrans.
    """

    if isinstance(s, str):
        return s.translate(table, deletechars)
    elif isinstance(s, unicode):
        table = dict((x, ord(table[x])) for x in range(256)
                     if ord(table[x]) != x)
        for c in deletechars:
            table[ord(c)] = None
        return s.translate(table)



def _array_from_file(infile, dtype, count, sep):
    """Create a numpy array from a file or a file-like object."""

    if isfile(infile):
        return np.fromfile(infile, dtype=dtype, count=count, sep=sep)
    else:
        # treat as file-like object with "read" method; this includes gzip file
        # objects, because numpy.fromfile just reads the compressed bytes from
        # their underlying file object, instead of the decompresed bytes
        read_size = np.dtype(dtype).itemsize * count
        s = infile.read(read_size)
        return np.fromstring(s, dtype=dtype, count=count, sep=sep)


def _array_to_file(arr, outfile):
    """Write a numpy array to a file or a file-like object."""

    if isfile(outfile):
        arr.tofile(outfile)
    else:
        # treat as file-like object with "write" method
        _write_string(outfile, arr.tostring())


def _write_string(f, s):
    """
    Write a string to a file, encoding to ASCII if the file is open in binary
    mode, or decoding if the file is open in text mode.
    """

    # Assume if the file object doesn't have a specific mode, that the mode is
    # binary
    mode = fileobj_mode(f)
    if mode:
        binmode = 'b' in mode
    else:
        binmode = True

    if binmode and isinstance(s, unicode):
        s = encode_ascii(s)
    elif not binmode and not isinstance(f, unicode):
        s = decode_ascii(s)
    f.write(s)


def _convert_array(array, dtype):
    """
    Converts an array to a new dtype--if the itemsize of the new dtype is
    the same as the old dtype, a view is returned.  Otherwise a new array must
    be created.
    """

    if array.dtype == dtype:
        return array
    elif array.dtype.itemsize == dtype.itemsize:
        return array.view(dtype)
    else:
        return array.astype(dtype)


def _unsigned_zero(dtype):
    """
    Given a numpy dtype, finds its "zero" point, which is exactly in the
    middle of its range.
    """

    assert dtype.kind == 'u'
    return 1 << (dtype.itemsize * 8 - 1)


def _is_pseudo_unsigned(dtype):
    return dtype.kind == 'u' and dtype.itemsize >= 2


def _is_int(val):
    return isinstance(val, (int, long, np.integer))


def _str_to_num(val):
    """Converts a given string to either an int or a float if necessary."""

    try:
        num = int(val)
    except ValueError:
        # If this fails then an exception should be raised anyways
        num = float(val)
    return num


def _pad_length(stringlen):
    """Bytes needed to pad the input stringlen to the next FITS block."""

    return (BLOCK_SIZE - (stringlen % BLOCK_SIZE)) % BLOCK_SIZE


def _normalize_slice(input, naxis):
    """
    Set the slice's start/stop in the regular range.
    """

    def _normalize(indx, npts):
        if indx < -npts:
            indx = 0
        elif indx < 0:
            indx += npts
        elif indx > npts:
            indx = npts
        return indx

    _start = input.start
    if _start is None:
        _start = 0
    elif _is_int(_start):
        _start = _normalize(_start, naxis)
    else:
        raise IndexError('Illegal slice %s; start must be integer.' % input)

    _stop = input.stop
    if _stop is None:
        _stop = naxis
    elif _is_int(_stop):
        _stop = _normalize(_stop, naxis)
    else:
        raise IndexError('Illegal slice %s; stop must be integer.' % input)

    if _stop < _start:
        raise IndexError('Illegal slice %s; stop < start.' % input)

    _step = input.step
    if _step is None:
        _step = 1
    elif _is_int(_step):
        if _step <= 0:
            raise IndexError('Illegal slice %s; step must be positive.'
                             % input)
    else:
        raise IndexError('Illegal slice %s; step must be integer.' % input)

    return slice(_start, _stop, _step)


def _tmp_name(input):
    """
    Create a temporary file name which should not already exist.  Use the
    directory of the input file as the base name of the mkstemp() output.
    """

    if input is not None:
        input = os.path.dirname(input)
    f, fn = tempfile.mkstemp(dir=input)
    os.close(f)
    return fn
