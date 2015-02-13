# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sundry function and class decorators."""

from __future__ import print_function


import functools
import inspect
import sys
import textwrap
import types
import warnings

from .codegen import make_function_with_signature
from .exceptions import (AstropyDeprecationWarning,
                         AstropyPendingDeprecationWarning)
from ..extern import six


__all__ = ['deprecated', 'deprecated_attribute', 'lazyproperty',
           'sharedmethod', 'wraps']


def deprecated(since, message='', name='', alternative='', pending=False,
               obj_type=None):
    """
    Used to mark a function or class as deprecated.

    To mark an attribute as deprecated, use `deprecated_attribute`.

    Parameters
    ------------
    since : str
        The release at which this API became deprecated.  This is
        required.

    message : str, optional
        Override the default deprecation message.  The format
        specifier ``func`` may be used for the name of the function,
        and ``alternative`` may be used in the deprecation message
        to insert the name of an alternative to the deprecated
        function. ``obj_type`` may be used to insert a friendly name
        for the type of object being deprecated.

    name : str, optional
        The name of the deprecated function or class; if not provided
        the name is automatically determined from the passed in
        function or class, though this is useful in the case of
        renamed functions, where the new function is just assigned to
        the name of the deprecated function.  For example::

            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative function or class name that the user may use in
        place of the deprecated object.  The deprecation warning will
        tell the user about this alternative if provided.

    pending : bool, optional
        If True, uses a AstropyPendingDeprecationWarning instead of a
        AstropyDeprecationWarning.

    obj_type : str, optional
        The type of this object, if the automatically determined one
        needs to be overridden.
    """

    method_types = (classmethod, staticmethod, types.MethodType)

    def deprecate_doc(old_doc, message):
        """
        Returns a given docstring with a deprecation message prepended
        to it.
        """
        if not old_doc:
            old_doc = ''
        old_doc = textwrap.dedent(old_doc).strip('\n')
        new_doc = (('\n.. deprecated:: %(since)s'
                    '\n    %(message)s\n\n' %
                    {'since': since, 'message': message.strip()}) + old_doc)
        if not old_doc:
            # This is to prevent a spurious 'unexpected unindent' warning from
            # docutils when the original docstring was blank.
            new_doc += r'\ '
        return new_doc

    def get_function(func):
        """
        Given a function or classmethod (or other function wrapper type), get
        the function object.
        """
        if isinstance(func, method_types):
            try:
                func = func.__func__
            except AttributeError:
                # classmethods in Python2.6 and below lack the __func__
                # attribute so we need to hack around to get it
                method = func.__get__(None, object)
                if isinstance(method, types.FunctionType):
                    # For staticmethods anyways the wrapped object is just a
                    # plain function (not a bound method or anything like that)
                    func = method
                elif hasattr(method, '__func__'):
                    func = method.__func__
                elif hasattr(method, 'im_func'):
                    func = method.im_func
                else:
                    # Nothing we can do really...  just return the original
                    # classmethod, etc.
                    return func
        return func

    def deprecate_function(func, message):
        """
        Returns a wrapped function that displays an
        ``AstropyDeprecationWarning`` when it is called.
        """

        if isinstance(func, method_types):
            func_wrapper = type(func)
        else:
            func_wrapper = lambda f: f

        func = get_function(func)

        def deprecated_func(*args, **kwargs):
            if pending:
                category = AstropyPendingDeprecationWarning
            else:
                category = AstropyDeprecationWarning

            warnings.warn(message, category, stacklevel=2)

            return func(*args, **kwargs)

        # If this is an extension function, we can't call
        # functools.wraps on it, but we normally don't care.
        # This crazy way to get the type of a wrapper descriptor is
        # straight out of the Python 3.3 inspect module docs.
        if type(func) != type(str.__dict__['__add__']):
            deprecated_func = functools.wraps(func)(deprecated_func)

        deprecated_func.__doc__ = deprecate_doc(
            deprecated_func.__doc__, message)

        return func_wrapper(deprecated_func)

    def deprecate_class(cls, message):
        """
        Returns a wrapper class with the docstrings updated and an
        __init__ function that will raise an
        ``AstropyDeprectationWarning`` warning when called.
        """
        # Creates a new class with the same name and bases as the
        # original class, but updates the dictionary with a new
        # docstring and a wrapped __init__ method.  __module__ needs
        # to be manually copied over, since otherwise it will be set
        # to *this* module (astropy.utils.misc).

        # This approach seems to make Sphinx happy (the new class
        # looks enough like the original class), and works with
        # extension classes (which functools.wraps does not, since
        # it tries to modify the original class).

        # We need to add a custom pickler or you'll get
        #     Can't pickle <class ..>: it's not found as ...
        # errors. Picklability is required for any class that is
        # documented by Sphinx.

        members = cls.__dict__.copy()

        members.update({
            '__doc__': deprecate_doc(cls.__doc__, message),
            '__init__': deprecate_function(get_function(cls.__init__),
                                           message),
        })

        return type(cls.__name__, cls.__bases__, members)

    def deprecate(obj, message=message, name=name, alternative=alternative,
                  pending=pending):
        if obj_type is None:
            if isinstance(obj, type):
                obj_type_name = 'class'
            elif inspect.isfunction(obj):
                obj_type_name = 'function'
            elif inspect.ismethod(obj) or isinstance(obj, method_types):
                obj_type_name = 'method'
            else:
                obj_type_name = 'object'
        else:
            obj_type_name = obj_type

        if not name:
            name = get_function(obj).__name__

        altmessage = ''
        if not message or type(message) == type(deprecate):
            if pending:
                message = ('The %(func)s %(obj_type)s will be deprecated in a '
                           'future version.')
            else:
                message = ('The %(func)s %(obj_type)s is deprecated and may '
                           'be removed in a future version.')
            if alternative:
                altmessage = '\n        Use %s instead.' % alternative

        message = ((message % {
            'func': name,
            'name': name,
            'alternative': alternative,
            'obj_type': obj_type_name}) +
            altmessage)

        if isinstance(obj, type):
            return deprecate_class(obj, message)
        else:
            return deprecate_function(obj, message)

    if type(message) == type(deprecate):
        return deprecate(message)

    return deprecate


def deprecated_attribute(name, since, message=None, alternative=None,
                         pending=False):
    """
    Used to mark a public attribute as deprecated.  This creates a
    property that will warn when the given attribute name is accessed.
    To prevent the warning (i.e. for internal code), use the private
    name for the attribute by prepending an underscore
    (i.e. ``self._name``).

    Parameters
    ----------
    name : str
        The name of the deprecated attribute.

    since : str
        The release at which this API became deprecated.  This is
        required.

    message : str, optional
        Override the default deprecation message.  The format
        specifier ``name`` may be used for the name of the attribute,
        and ``alternative`` may be used in the deprecation message
        to insert the name of an alternative to the deprecated
        function.

    alternative : str, optional
        An alternative attribute that the user may use in place of the
        deprecated attribute.  The deprecation warning will tell the
        user about this alternative if provided.

    pending : bool, optional
        If True, uses a AstropyPendingDeprecationWarning instead of a
        AstropyDeprecationWarning.

    Examples
    --------

    ::

        class MyClass:
            # Mark the old_name as deprecated
            old_name = misc.deprecated_attribute('old_name', '0.1')

            def method(self):
                self._old_name = 42
    """
    private_name = '_' + name

    @deprecated(since, name=name, obj_type='attribute')
    def get(self):
        return getattr(self, private_name)

    @deprecated(since, name=name, obj_type='attribute')
    def set(self, val):
        setattr(self, private_name, val)

    @deprecated(since, name=name, obj_type='attribute')
    def delete(self):
        delattr(self, private_name)

    return property(get, set, delete)


class lazyproperty(object):
    """
    Works similarly to property(), but computes the value only once.

    This essentially memorizes the value of the property by storing the result
    of its computation in the ``__dict__`` of the object instance.  This is
    useful for computing the value of some property that should otherwise be
    invariant.  For example::

        >>> class LazyTest(object):
        ...     @lazyproperty
        ...     def complicated_property(self):
        ...         print('Computing the value for complicated_property...')
        ...         return 42
        ...
        >>> lt = LazyTest()
        >>> lt.complicated_property
        Computing the value for complicated_property...
        42
        >>> lt.complicated_property
        42

    As the example shows, the second time ``complicated_property`` is accessed,
    the ``print`` statement is not executed.  Only the return value from the
    first access off ``complicated_property`` is returned.

    If a setter for this property is defined, it will still be possible to
    manually update the value of the property, if that capability is desired.

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
        self._key = self._fget.__name__

    def __get__(self, obj, owner=None):
        if obj is None:
            return self
        try:
            return obj.__dict__[self._key]
        except KeyError:
            val = self._fget(obj)
            obj.__dict__[self._key] = val
            return val

    def __set__(self, obj, val):
        obj_dict = obj.__dict__
        if self._fset:
            ret = self._fset(obj, val)
            if ret is not None and obj_dict.get(self._key) is ret:
                # By returning the value set the setter signals that it took
                # over setting the value in obj.__dict__; this mechanism allows
                # it to override the input value
                return
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self._fdel:
            self._fdel(obj)
        if self._key in obj.__dict__:
            del obj.__dict__[self._key]

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
        for k, v in six.iteritems(cls_ns):
            if v is self:
                property_name = k
                break

        cls_ns[property_name] = lazyproperty(*args)

        return cls_ns[property_name]


class sharedmethod(classmethod):
    """
    This is a method decorator that allows both an instancemethod and a
    `classmethod` to share the same name.

    When using `sharedmethod` on a method defined in a class's body, it
    may be called on an instance, or on a class.  In the former case it
    behaves like a normal instance method (a reference to the instance is
    automatically passed as the first ``self`` argument of the method)::

        >>> class Example(object):
        ...     @sharedmethod
        ...     def identify(self, *args):
        ...         print('self was', self)
        ...         print('additional args were', args)
        ...
        >>> ex = Example()
        >>> ex.identify(1, 2)
        self was <astropy.utils.decorators.Example object at 0x...>
        additional args were (1, 2)

    In the latter case, when the `sharedmethod` is called directly from a
    class, it behaves like a `classmethod`::

        >>> Example.identify(3, 4)
        self was <class 'astropy.utils.decorators.Example'>
        additional args were (3, 4)

    This also supports a more advanced usage, where the `classmethod`
    implementation can be written separately.  If the class's *metaclass*
    has a method of the same name as the `sharedmethod`, the version on
    the metaclass is delegated to::

        >>> from astropy.extern.six import add_metaclass
        >>> class ExampleMeta(type):
        ...     def identify(self):
        ...         print('this implements the {0}.identify '
        ...               'classmethod'.format(self.__name__))
        ...
        >>> @add_metaclass(ExampleMeta)
        ... class Example(object):
        ...     @sharedmethod
        ...     def identify(self):
        ...         print('this implements the instancemethod')
        ...
        >>> Example().identify()
        this implements the instancemethod
        >>> Example.identify()
        this implements the Example.identify classmethod
    """

    if sys.version_info[:2] < (2, 7):
        # Workaround for Python 2.6 which does not have classmethod.__func__
        @property
        def __func__(self):
            try:
                meth = classmethod.__get__(self, self.__obj__,
                                           self.__objtype__)
            except AttributeError:
                # self.__obj__ not set when called from __get__, but then it
                # doesn't matter anyways
                meth = classmethod.__get__(self, None, object)
            return meth.__func__

        def __getobjwrapper(orig_get):
            """
            Used to temporarily set/unset self.__obj__ and self.__objtype__
            for use by __func__.
            """
            def __get__(self, obj, objtype=None):
                self.__obj__ = obj
                self.__objtype__ = objtype

                try:
                    return orig_get(self, obj, objtype)
                finally:
                    del self.__obj__
                    del self.__objtype__

            return __get__
    else:
        def __getobjwrapper(func):
            return func

    @__getobjwrapper
    def __get__(self, obj, objtype=None):
        if obj is None:
            mcls = type(objtype)
            clsmeth = getattr(mcls, self.__func__.__name__, None)
            if callable(clsmeth):
                if isinstance(clsmeth, types.MethodType):
                    # This case will generally only apply on Python 2, which
                    # uses MethodType for unbound methods; Python 3 has no
                    # particular concept of unbound methods and will just
                    # return a function
                    func = clsmeth.__func__
                else:
                    func = clsmeth
            else:
                func = self.__func__

            return self._make_method(func, objtype)
        else:
            return self._make_method(self.__func__, obj)

    del __getobjwrapper

    if six.PY3:
        # The 'instancemethod' type of Python 2 and the method type of
        # Python 3 have slightly different constructors
        @staticmethod
        def _make_method(func, instance):
            return types.MethodType(func, instance)
    else:
        @staticmethod
        def _make_method(func, instance):
            return types.MethodType(func, instance, type(instance))


def wraps(wrapped, assigned=functools.WRAPPER_ASSIGNMENTS,
          updated=functools.WRAPPER_UPDATES):
    """
    An alternative to `functools.wraps` which also preserves the original
    function's call signature by way of
    `~astropy.utils.codegen.make_function_with_signature`.

    The documentation for the original `functools.wraps` follows:

    """

    def wrapper(func):
        func = make_function_with_signature(func, name=wrapped.__name__,
                                            **_get_function_args(wrapped))
        func = functools.update_wrapper(func, wrapped, assigned=assigned,
                                        updated=updated)
        return func

    return wrapper


wraps.__doc__ += functools.wraps.__doc__


if six.PY3:
    def _get_function_args(func):
        """
        Utility function for `wraps`.

        Reads the argspec for the given function and converts it to arguments
        for `make_function_with_signature`.  This requires different
        implementations on Python 2 versus Python 3.
        """

        argspec = inspect.getfullargspec(func)

        if argspec.defaults:
            args = argspec.args[:-len(argspec.defaults)]
            kwargs = zip(argspec.args[len(args):], argspec.defaults)
        else:
            args = argspec.args
            kwargs = []

        if argspec.kwonlyargs:
            kwargs.extend((argname, argspec.kwonlydefaults[argname])
                          for argname in argspec.kwonlyargs)

        return {'args': args, 'kwargs': kwargs, 'varargs': argspec.varargs,
                'varkwargs': argspec.varkw}
else:
    def _get_function_args(func):
        """
        Utility function for `wraps`.

        Reads the argspec for the given function and converts it to arguments
        for `make_function_with_signature`.  This requires different
        implementations on Python 2 versus Python 3.
        """

        argspec = inspect.getargspec(func)

        if argspec.defaults:
            args = argspec.args[:-len(argspec.defaults)]
            kwargs = zip(argspec.args[len(args):], argspec.defaults)
        else:
            args = argspec.args
            kwargs = {}

        return {'args': args, 'kwargs': kwargs, 'varargs': argspec.varargs,
                'varkwargs': argspec.keywords}
