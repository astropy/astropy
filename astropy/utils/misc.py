# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A "grab bag" of relatively small general-purpose utilities that don't have
a clear module/package to live in.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import abc
import contextlib
import difflib
import inspect
import json
import os
import re
import signal
import sys
import traceback
import unicodedata

from collections import defaultdict

from ..extern import six
from ..extern.six.moves import urllib
from ..utils.compat.odict import OrderedDict


__all__ = ['isiterable', 'silence', 'format_exception', 'NumpyRNGContext',
           'find_api_page', 'is_path_hidden', 'walk_skip_hidden',
           'JsonCustomEncoder', 'indent', 'InheritDocstrings',
           'OrderedDescriptor', 'OrderedDescriptorContainer']


def isiterable(obj):
    """Returns `True` if the given object is iterable."""

    try:
        iter(obj)
        return True
    except TypeError:
        return False


def indent(s, shift=1, width=4):
    """Indent a block of text.  The indentation is applied to each line."""

    indented = '\n'.join(' ' * (width * shift) + l if l else ''
                         for l in s.splitlines())
    if s[-1] == '\n':
        indented += '\n'

    return indented


class _DummyFile(object):
    """A noop writeable object."""

    def write(self, s):
        pass


@contextlib.contextmanager
def silence():
    """A context manager that silences sys.stdout and sys.stderr."""

    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = _DummyFile()
    sys.stderr = _DummyFile()
    yield
    sys.stdout = old_stdout
    sys.stderr = old_stderr


def format_exception(msg, *args, **kwargs):
    """
    Given an exception message string, uses new-style formatting arguments
    ``{filename}``, ``{lineno}``, ``{func}`` and/or ``{text}`` to fill in
    information about the exception that occurred.  For example:

        try:
            1/0
        except:
            raise ZeroDivisionError(
                format_except('A divide by zero occurred in {filename} at '
                              'line {lineno} of function {func}.'))

    Any additional positional or keyword arguments passed to this function are
    also used to format the message.

    .. note::
        This uses `sys.exc_info` to gather up the information needed to fill
        in the formatting arguments. Python 2.x and 3.x have slightly
        different behavior regarding `sys.exc_info` (the latter will not
        carry it outside a handled exception), so it's not wise to use this
        outside of an ``except`` clause - if it is, this will substitute
        '<unkown>' for the 4 formatting arguments.
    """

    tb = traceback.extract_tb(sys.exc_info()[2], limit=1)
    if len(tb) > 0:
        filename, lineno, func, text = tb[0]
    else:
        filename = lineno = func = text = '<unknown>'

    return msg.format(*args, filename=filename, lineno=lineno, func=func,
                      text=text, **kwargs)


class NumpyRNGContext(object):
    """
    A context manager (for use with the ``with`` statement) that will seed the
    numpy random number generator (RNG) to a specific value, and then restore
    the RNG state back to whatever it was before.

    This is primarily intended for use in the astropy testing suit, but it
    may be useful in ensuring reproducibility of Monte Carlo simulations in a
    science context.

    Parameters
    ----------
    seed : int
        The value to use to seed the numpy RNG

    Examples
    --------
    A typical use case might be::

        with NumpyRNGContext(<some seed value you pick>):
            from numpy import random

            randarr = random.randn(100)
            ... run your test using `randarr` ...

        #Any code using numpy.random at this indent level will act just as it
        #would have if it had been before the with statement - e.g. whatever
        #the default seed is.


    """
    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):
        from numpy import random

        self.startstate = random.get_state()
        random.seed(self.seed)

    def __exit__(self, exc_type, exc_value, traceback):
        from numpy import random

        random.set_state(self.startstate)


def find_api_page(obj, version=None, openinbrowser=True, timeout=None):
    """
    Determines the URL of the API page for the specified object, and
    optionally open that page in a web browser.

    .. note::
        You must be connected to the internet for this to function even if
        ``openinbrowser`` is `False`, unless you provide a local version of
        the documentation to ``version`` (e.g., ``file:///path/to/docs``).

    Parameters
    ----------
    obj
        The object to open the docs for or its fully-qualified name
        (as a str).
    version : str
        The doc version - either a version number like '0.1', 'dev' for
        the development/latest docs, or a URL to point to a specific
        location that should be the *base* of the documentation. Defaults to
        latest if you are on aren't on a release, otherwise, the version you
        are on.
    openinbrowser : bool
        If `True`, the `webbrowser` package will be used to open the doc
        page in a new web browser window.
    timeout : number, optional
        The number of seconds to wait before timing-out the query to
        the astropy documentation.  If not given, the default python
        stdlib timeout will be used.

    Returns
    -------
    url : str
        The loaded URL

    Raises
    ------
    ValueError
        If the documentation can't be found

    """
    import webbrowser

    from zlib import decompress

    if (not isinstance(obj, six.string_types) and
            hasattr(obj, '__module__') and
            hasattr(obj, '__name__')):
        obj = obj.__module__ + '.' + obj.__name__
    elif inspect.ismodule(obj):
        obj = obj.__name__

    if version is None:
        from .. import version

        if version.release:
            version = 'v' + version.version
        else:
            version = 'dev'

    if '://' in version:
        if version.endswith('index.html'):
            baseurl = version[:-10]
        elif version.endswith('/'):
            baseurl = version
        else:
            baseurl = version + '/'
    elif version == 'dev' or version == 'latest':
        baseurl = 'http://devdocs.astropy.org/'
    else:
        baseurl = 'http://docs.astropy.org/en/{vers}/'.format(vers=version)

    if timeout is None:
        uf = urllib.request.urlopen(baseurl + 'objects.inv')
    else:
        uf = urllib.request.urlopen(baseurl + 'objects.inv', timeout=timeout)

    try:
        # we read these lines so that `oistr` only gets the compressed
        # contents, not the header information
        isvers = uf.readline().rstrip().decode('utf-8')  # intersphinx version line
        proj = uf.readline().rstrip().decode('utf-8')  # project name
        vers = uf.readline().rstrip().decode('utf-8')  # project version
        uf.readline().rstrip().decode('utf-8')
        oistr = uf.read()
    finally:
        uf.close()

    oistr = decompress(oistr).decode('utf-8')

    resurl = None

    for l in oistr.strip().splitlines():
        ls = l.split()
        name = ls[0]
        loc = ls[3]
        if loc.endswith('$'):
            loc = loc[:-1] + name

        if name == obj:
            resurl = baseurl + loc
            break

    if resurl is None:
        raise ValueError('Could not find the docs for the object {obj}'.format(obj=obj))
    elif openinbrowser:
        webbrowser.open(resurl)

    return resurl


def signal_number_to_name(signum):
    """
    Given an OS signal number, returns a signal name.  If the signal
    number is unknown, returns ``'UNKNOWN'``.
    """
    # Since these numbers and names are platform specific, we use the
    # builtin signal module and build a reverse mapping.

    signal_to_name_map = dict(
        (k, v) for v, k in signal.__dict__.iteritems() if v.startswith('SIG'))

    return signal_to_name_map.get(signum, 'UNKNOWN')


if sys.platform == 'win32':
    import ctypes

    def _has_hidden_attribute(filepath):
        """
        Returns True if the given filepath has the hidden attribute on
        MS-Windows.  Based on a post here:
        http://stackoverflow.com/questions/284115/cross-platform-hidden-file-detection
        """
        if isinstance(filepath, bytes):
            filepath = filepath.decode(sys.getfilesystemencoding())
        try:
            attrs = ctypes.windll.kernel32.GetFileAttributesW(filepath)
            assert attrs != -1
            result = bool(attrs & 2)
        except (AttributeError, AssertionError):
            result = False
        return result
else:
    def _has_hidden_attribute(filepath):
        return False


def is_path_hidden(filepath):
    """
    Determines if a given file or directory is hidden.

    Parameters
    ----------
    filepath : str
        The path to a file or directory

    Returns
    -------
    hidden : bool
        Returns `True` if the file is hidden
    """
    name = os.path.basename(os.path.abspath(filepath))
    if isinstance(name, bytes):
        is_dotted = name.startswith(b'.')
    else:
        is_dotted = name.startswith('.')
    return is_dotted or _has_hidden_attribute(filepath)


def walk_skip_hidden(top, onerror=None, followlinks=False):
    """
    A wrapper for `os.walk` that skips hidden files and directories.

    This function does not have the parameter ``topdown`` from
    `os.walk`: the directories must always be recursed top-down when
    using this function.

    See also
    --------
    os.walk : For a description of the parameters
    """
    for root, dirs, files in os.walk(
            top, topdown=True, onerror=onerror,
            followlinks=followlinks):
        # These lists must be updated in-place so os.walk will skip
        # hidden directories
        dirs[:] = [d for d in dirs if not is_path_hidden(d)]
        files[:] = [f for f in files if not is_path_hidden(f)]
        yield root, dirs, files


class JsonCustomEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder
    does not do.

    This includes:

        * Numpy array or number
        * Complex number
        * Set
        * Bytes (Python 3)

    Examples
    --------
    >>> import json
    >>> import numpy as np
    >>> from astropy.utils.misc import JsonCustomEncoder
    >>> json.dumps(np.arange(3), cls=JsonCustomEncoder)
    '[0, 1, 2]'

    """
    def default(self, obj):
        import numpy as np
        if isinstance(obj, (np.ndarray, np.number)):
            return obj.tolist()
        elif isinstance(obj, (complex, np.complex)):
            return [obj.real, obj.imag]
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, bytes):  # pragma: py3
            return obj.decode()
        return json.JSONEncoder.default(self, obj)


def strip_accents(s):
    """
    Remove accents from a Unicode string.

    This helps with matching "ångström" to "angstrom", for example.
    """
    return ''.join(
        c for c in unicodedata.normalize('NFD', s)
        if unicodedata.category(c) != 'Mn')


def did_you_mean(s, candidates, n=3, cutoff=0.8, fix=None):
    """
    When a string isn't found in a set of candidates, we can be nice
    to provide a list of alternatives in the exception.  This
    convenience function helps to format that part of the exception.

    Parameters
    ----------
    s : str

    candidates : sequence of str or dict of str keys

    n : int
        The maximum number of results to include.  See
        `difflib.get_close_matches`.

    cutoff : float
        In the range [0, 1]. Possibilities that don't score at least
        that similar to word are ignored.  See
        `difflib.get_close_matches`.

    fix : callable
        A callable to modify the results after matching.  It should
        take a single string and return a sequence of strings
        containing the fixed matches.

    Returns
    -------
    message : str
        Returns the string "Did you mean X, Y, or Z?", or the empty
        string if no alternatives were found.
    """
    if isinstance(s, six.text_type):
        s = strip_accents(s)
    s_lower = s.lower()

    # Create a mapping from the lower case name to all capitalization
    # variants of that name.
    candidates_lower = {}
    for candidate in candidates:
        candidate_lower = candidate.lower()
        candidates_lower.setdefault(candidate_lower, [])
        candidates_lower[candidate_lower].append(candidate)

    # The heuristic here is to first try "singularizing" the word.  If
    # that doesn't match anything use difflib to find close matches in
    # original, lower and upper case.
    if s_lower.endswith('s') and s_lower[:-1] in candidates_lower:
        matches = [s_lower[:-1]]
    else:
        matches = difflib.get_close_matches(
            s_lower, candidates_lower, n=n, cutoff=cutoff)

    if len(matches):
        capitalized_matches = set()
        for match in matches:
            capitalized_matches.update(candidates_lower[match])
        matches = capitalized_matches

        if fix is not None:
            mapped_matches = []
            for match in matches:
                mapped_matches.extend(fix(match))
            matches = mapped_matches

        matches = list(set(matches))
        matches = sorted(matches)

        if len(matches) == 1:
            matches = matches[0]
        else:
            matches = (', '.join(matches[:-1]) + ' or ' +
                       matches[-1])
        return 'Did you mean {0}?'.format(matches)

    return ''


class InheritDocstrings(type):
    """
    This metaclass makes methods of a class automatically have their
    docstrings filled in from the methods they override in the base
    class.

    If the class uses multiple inheritance, the docstring will be
    chosen from the first class in the bases list, in the same way as
    methods are normally resolved in Python.  If this results in
    selecting the wrong docstring, the docstring will need to be
    explicitly included on the method.

    For example::

        >>> from astropy.utils.misc import InheritDocstrings
        >>> from astropy.extern import six
        >>> @six.add_metaclass(InheritDocstrings)
        ... class A(object):
        ...     def wiggle(self):
        ...         "Wiggle the thingamajig"
        ...         pass
        >>> class B(A):
        ...     def wiggle(self):
        ...         pass
        >>> B.wiggle.__doc__
        u'Wiggle the thingamajig'
    """

    def __init__(cls, name, bases, dct):
        def is_public_member(key):
            return (
                (key.startswith('__') and key.endswith('__')
                 and len(key) > 4) or
                not key.startswith('_'))

        for key, val in six.iteritems(dct):
            if (inspect.isfunction(val) and
                is_public_member(key) and
                val.__doc__ is None):
                for base in cls.__mro__[1:]:
                    super_method = getattr(base, key, None)
                    if super_method is not None:
                        val.__doc__ = super_method.__doc__
                        break

        super(InheritDocstrings, cls).__init__(name, bases, dct)


@six.add_metaclass(abc.ABCMeta)
class OrderedDescriptor(object):
    """
    Base class for descriptors whose order in the class body should be
    preserved.  Intended for use in concert with the
    `OrderedDescriptorContainer` metaclass.

    Subclasses of `OrderedDescriptor` must define a value for a class attribute
    called ``_class_attribute_``.  This is the name of a class attribute on the
    *container* class for these descriptors, which will be set to an
    `~collections.OrderedDict` at class creation time.  This
    `~collections.OrderedDict` will contain a mapping of all class attributes
    that were assigned instances of the `OrderedDescriptor` subclass, to the
    instances themselves.  See the documentation for
    `OrderedDescriptorContainer` for a concrete example.

    Optionally, subclasses of `OrderedDescriptor` may define a value for a
    class attribute called ``_name_attribute_``.  This should be the name of
    an attribute on instances of the subclass.  When specified, during
    creation of a class containing these descriptors, the name attribute on
    each instance will be set to the name of the class attribute it was
    assigned to on the class.

    .. note::

        Although this class is intended for use with *descriptors* (i.e.
        classes that define any of the ``__get__``, ``__set__``, or
        ``__delete__`` magic methods), this base class is not itself a
        descriptor, and technically this could be used for classes that are
        not descriptors too.  However, use with descriptors is the original
        intended purpose.
    """

    # This id increments for each OrderedDescriptor instance created, so they
    # are always ordered in the order they were created.  Class bodies are
    # guaranteed to be executed from top to bottom.  Not sure if this is
    # thread-safe though.
    _nextid = 1

    _class_attribute_ = abc.abstractproperty()
    """
    Subclasses should define this attribute to the name of an attribute on
    classes containing this subclass.  That attribute will contain the mapping
    of all instances of that `OrderedDescriptor` subclass defined in the class
    body.  If the same descriptor needs to be used with different classes,
    each with different names of this attribute, multiple subclasses will be
    needed.
    """

    _name_attribute_ = None
    """
    Subclasses may optionally define this attribute to specify the name of an
    attribute on instances of the class that should be filled with the
    instance's attribute name at class creation time.
    """

    def __init__(self, *args, **kwargs):
        # The _nextid attribute is shared across all subclasses so that
        # different subclasses of OrderedDescriptors can be sorted correctly
        # between themselves
        self.__order = OrderedDescriptor._nextid
        OrderedDescriptor._nextid += 1
        super(OrderedDescriptor, self).__init__()

    def __lt__(self, other):
        """
        Defined for convenient sorting of `OrderedDescriptor` instances, which
        are defined to sort in their creation order.
        """

        if (isinstance(self, OrderedDescriptor) and
                isinstance(other, OrderedDescriptor)):
            try:
                return self.__order < other.__order
            except AttributeError:
                raise RuntimeError(
                    'Could not determine ordering for {0} and {1}; at least '
                    'one of them is not calling super().__init__ in its '
                    '__init__.'.format(self, other))
        else:
            return NotImplemented


class OrderedDescriptorContainer(type):
    """
    Classes should use this metaclass if they wish to use `OrderedDescriptor`
    attributes, which are class attributes that "remember" the order in which
    they were defined in the class body.

    Every subclass of `OrderedDescriptor` has an attribute called
    ``_class_attribute_``.  For example, if we have

    .. code:: python

        class ExampleDecorator(OrderedDescriptor):
            _class_attribute_ = '_examples_'

    Then when a class with the `OrderedDescriptorContainer` metaclass is
    created, it will automatically be assigned a class attribute ``_examples_``
    referencing an `~collections.OrderedDict` containing all instances of
    ``ExampleDecorator`` defined in the class body, mapped to by the names of
    the attributes they were assigned to.

    When subclassing a class with this metaclass, the descriptor dict (i.e.
    ``_examples_`` in the above example) will *not* contain descriptors
    inherited from the base class.  That is, this only works by default with
    decorators explicitly defined in the class body.  However, the subclass
    *may* define an attribute ``_inherit_decorators_`` which lists
    `OrderedDescriptor` classes that *should* be added from base classes.
    See the examples section below for an example of this.

    Examples
    --------

    >>> from astropy.extern import six
    >>> from astropy.utils import OrderedDescriptor, OrderedDescriptorContainer
    >>> class TypedAttribute(OrderedDescriptor):
    ...     \"\"\"
    ...     Attributes that may only be assigned objects of a specific type,
    ...     or subclasses thereof.  For some reason we care about their order.
    ...     \"\"\"
    ...
    ...     _class_attribute_ = 'typed_attributes'
    ...     _name_attribute_ = 'name'
    ...     # A default name so that instances not attached to a class can
    ...     # still be repr'd; useful for debugging
    ...     name = '<unbound>'
    ...
    ...     def __init__(self, type):
    ...         # Make sure not to forget to call the super __init__
    ...         super(TypedAttribute, self).__init__()
    ...         self.type = type
    ...
    ...     def __get__(self, obj, objtype=None):
    ...         if obj is None:
    ...             return self
    ...         if self.name in obj.__dict__:
    ...             return obj.__dict__[self.name]
    ...         else:
    ...             raise AttributeError(self.name)
    ...
    ...     def __set__(self, obj, value):
    ...         if not isinstance(value, self.type):
    ...             raise ValueError('{0}.{1} must be of type {2!r}'.format(
    ...                 obj.__class__.__name__, self.name, self.type))
    ...         obj.__dict__[self.name] = value
    ...
    ...     def __delete__(self, obj):
    ...         if self.name in obj.__dict__:
    ...             del obj.__dict__[self.name]
    ...         else:
    ...             raise AttributeError(self.name)
    ...
    ...     def __repr__(self):
    ...         if isinstance(self.type, tuple) and len(self.type) > 1:
    ...             typestr = '({0})'.format(
    ...                 ', '.join(t.__name__ for t in self.type))
    ...         else:
    ...             typestr = self.type.__name__
    ...         return '<{0}(name={1}, type={2})>'.format(
    ...                 self.__class__.__name__, self.name, typestr)
    ...

    Now let's create an example class that uses this ``TypedAttribute``::

        >>> @six.add_metaclass(OrderedDescriptorContainer)
        ... class Point2D(object):
        ...     x = TypedAttribute((float, int))
        ...     y = TypedAttribute((float, int))
        ...
        ...     def __init__(self, x, y):
        ...         self.x, self.y = x, y
        ...
        >>> p1 = Point2D(1.0, 2.0)
        >>> p1.x
        1.0
        >>> p1.y
        2.0
        >>> p2 = Point2D('a', 'b')  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: Point2D.x must be of type (float, int>)

    We see that ``TypedAttribute`` works more or less as advertised, but
    there's nothing special about that.  Let's see what
    `OrderedDescriptorContainer` did for us::

        >>> Point2D.typed_attributes
        OrderedDict([('x', <TypedAttribute(name=x, type=(float, int))>),
        ('y', <TypedAttribute(name=y, type=(float, int))>)])

    If we create a subclass, it does *not* by default add inherited descriptors
    to ``typed_attributes``::

        >>> class Point3D(Point2D):
        ...     z = TypedAttribute((float, int))
        ...
        >>> Point3D.typed_attributes
        OrderedDict([('z', <TypedAttribute(name=z, type=(float, int))>)])

    However, if we specify ``_inherit_descriptors_`` from ``Point2D`` then
    it will do so::

        >>> class Point3D(Point2D):
        ...     _inherit_descriptors_ = (TypedAttribute,)
        ...     z = TypedAttribute((float, int))
        ...
        >>> Point3D.typed_attributes
        OrderedDict([('x', <TypedAttribute(name=x, type=(float, int))>),
        ('y', <TypedAttribute(name=y, type=(float, int))>),
        ('z', <TypedAttribute(name=z, type=(float, int))>)])

    .. note::

        Hopefully it is clear from these examples that this construction
        also allows a class of type `OrderedDescriptorContainer` to use
        multiple different `OrderedDescriptor` classes simultaneously.
    """

    _inherit_descriptors_ = ()

    def __init__(cls, cls_name, bases, members):
        descriptors = defaultdict(list)
        seen = set()
        inherit_descriptors = ()
        descr_bases = {}

        for mro_cls in cls.__mro__:
            for name, obj in mro_cls.__dict__.items():
                if name in seen:
                    # Checks if we've already seen an attribute of the given
                    # name (if so it will override anything of the same name in
                    # any base class)
                    continue

                seen.add(name)

                if (not isinstance(obj, OrderedDescriptor) or
                        (inherit_descriptors and
                            not isinstance(obj, inherit_descriptors))):
                    # The second condition applies when checking any
                    # subclasses, to see if we can inherit any descriptors of
                    # the given type from subclasses (by default inheritance is
                    # disabled unless the class has _inherit_descriptors_
                    # defined)
                    continue

                if obj._name_attribute_ is not None:
                    setattr(obj, obj._name_attribute_, name)

                # Don't just use the descriptor's class directly; instead go
                # through its MRO and find the class on which _class_attribute_
                # is defined directly.  This way subclasses of some
                # OrderedDescriptor *may* override _class_attribute_ and have
                # its own _class_attribute_, but by default all subclasses of
                # some OrderedDescriptor are still grouped together
                # TODO: It might be worth clarifying this in the docs
                if obj.__class__ not in descr_bases:
                    for obj_cls_base in obj.__class__.__mro__:
                        if '_class_attribute_' in obj_cls_base.__dict__:
                            descr_bases[obj.__class__] = obj_cls_base
                            descriptors[obj_cls_base].append((obj, name))
                            break
                else:
                    # Make sure to put obj first for sorting purposes
                    obj_cls_base = descr_bases[obj.__class__]
                    descriptors[obj_cls_base].append((obj, name))

            if not (isinstance(mro_cls, type(cls)) and
                        mro_cls._inherit_descriptors_):
                # If _inherit_descriptors_ is undefined then we don't inherit
                # any OrderedDescriptors from any of the base classes, and
                # there's no reason to continue through the MRO
                break
            else:
                inherit_descriptors = mro_cls._inherit_descriptors_

        for descriptor_cls, instances in descriptors.items():
            instances.sort()
            instances = OrderedDict((key, value) for value, key in instances)
            setattr(cls, descriptor_cls._class_attribute_, instances)

        super(OrderedDescriptorContainer, cls).__init__(cls_name, bases,
                                                        members)
