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


__all__ = ['deprecated', 'deprecated_attribute', 'classproperty',
           'lazyproperty', 'sharedmethod', 'wraps']


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

        return type(cls)(cls.__name__, cls.__bases__, members)

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


# TODO: This can still be made to work for setters by implementing an
# accompanying metaclass that supports it; we just don't need that right this
# second
class classproperty(property):
    """
    Similar to `property`, but allows class-level properties.  That is,
    a property whose getter is like a `classmethod`.

    The wrapped method may explicitly use the `classmethod` decorator (which
    must become before this decorator), or the `classmethod` may be omitted
    (it is implicit through use of this decorator).

    .. note::

        classproperty only works for *read-only* properties.  It does not
        currently allow writeable/deleteable properties, due to subtleties of how
        Python descriptors work.  In order to implement such properties on a class
        a metaclass for that class must be implemented.

    Parameters
    ----------
    fget : callable
        The function that computes the value of this property (in particular,
        the function when this is used as a decorator) a la `property`.

    doc : str, optional
        The docstring for the property--by default inherited from the getter
        function.

    lazy : bool, optional
        If True, caches the value returned by the first call to the getter
        function, so that it is only called once (used for lazy evaluation
        of an attribute).  This is analogous to `lazyproperty`.  The ``lazy``
        argument can also be used when `classproperty` is used as a decorator
        (see the third example below).  When used in the decorator syntax this
        *must* be passed in as a keyword argument.

    Examples
    --------

    ::

        >>> class Foo(object):
        ...     _bar_internal = 1
        ...     @classproperty
        ...     def bar(cls):
        ...         return cls._bar_internal + 1
        ...
        >>> Foo.bar
        2
        >>> foo_instance = Foo()
        >>> foo_instance.bar
        2
        >>> foo_instance._bar_internal = 2
        >>> foo_instance.bar  # Ignores instance attributes
        2

    As previously noted, a `classproperty` is limited to implementing
    read-only attributes::

        >>> class Foo(object):
        ...     _bar_internal = 1
        ...     @classproperty
        ...     def bar(cls):
        ...         return cls._bar_internal
        ...     @bar.setter
        ...     def bar(cls, value):
        ...         cls._bar_internal = value
        ...
        Traceback (most recent call last):
        ...
        NotImplementedError: classproperty can only be read-only; use a
        metaclass to implement modifiable class-level properties

    When the ``lazy`` option is used, the getter is only called once::

        >>> class Foo(object):
        ...     @classproperty(lazy=True)
        ...     def bar(cls):
        ...         print("Performing complicated calculation")
        ...         return 1
        ...
        >>> Foo.bar
        Performing complicated calculation
        1
        >>> Foo.bar
        1

    If a subclass inherits a lazy `classproperty` the property is still
    re-evaluated for the subclass::

        >>> class FooSub(Foo):
        ...     pass
        ...
        >>> FooSub.bar
        Performing complicated calculation
        1
        >>> FooSub.bar
        1
    """

    def __new__(cls, fget=None, doc=None, lazy=False):
        if fget is None:
            # Being used as a decorator--return a wrapper that implements
            # decorator syntax
            def wrapper(func):
                return cls(func, lazy=lazy)

            return wrapper

        return super(classproperty, cls).__new__(cls)

    def __init__(self, fget, doc=None, lazy=False):
        self._lazy = lazy
        if lazy:
            self._cache = {}
        fget = self._wrap_fget(fget)

        super(classproperty, self).__init__(fget=fget, doc=doc)

        # There is a buglet in Python where self.__doc__ doesn't
        # get set properly on instances of property subclasses if
        # the doc argument was used rather than taking the docstring
        # from fget
        if doc is not None:
            self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if self._lazy and objtype in self._cache:
            return self._cache[objtype]

        if objtype is not None:
            # The base property.__get__ will just return self here;
            # instead we pass objtype through to the original wrapped
            # function (which takes the class as its sole argument)
            val = self.fget.__wrapped__(objtype)
        else:
            val = super(classproperty, self).__get__(obj, objtype=objtype)

        if self._lazy:
            if objtype is None:
                objtype = obj.__class__

            self._cache[objtype] = val

        return val

    def getter(self, fget):
        return super(classproperty, self).getter(self._wrap_fget(fget))

    def setter(self, fset):
        raise NotImplementedError(
            "classproperty can only be read-only; use a metaclass to "
            "implement modifiable class-level properties")

    def deleter(self, fdel):
        raise NotImplementedError(
            "classproperty can only be read-only; use a metaclass to "
            "implement modifiable class-level properties")

    @staticmethod
    def _wrap_fget(orig_fget):
        if isinstance(orig_fget, classmethod):
            orig_fget = orig_fget.__func__

        # Using stock functools.wraps instead of the fancier version
        # found later in this module, which is overkill for this purpose

        @wraps(orig_fget)
        def fget(obj):
            return orig_fget(obj.__class__)

        # Set the __wrapped__ attribute manually for support on Python 2
        fget.__wrapped__ = orig_fget

        return fget


class lazyproperty(property):
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
        super(lazyproperty, self).__init__(fget, fset, fdel, doc)
        self._key = self.fget.__name__

    def __get__(self, obj, owner=None):
        try:
            return obj.__dict__[self._key]
        except KeyError:
            val = self.fget(obj)
            obj.__dict__[self._key] = val
            return val
        except AttributeError:
            if obj is None:
                return self
            raise

    def __set__(self, obj, val):
        obj_dict = obj.__dict__
        if self.fset:
            ret = self.fset(obj, val)
            if ret is not None and obj_dict.get(self._key) is ret:
                # By returning the value set the setter signals that it took
                # over setting the value in obj.__dict__; this mechanism allows
                # it to override the input value
                return
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        if self._key in obj.__dict__:
            del obj.__dict__[self._key]


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
          updated=functools.WRAPPER_UPDATES, exclude_args=()):
    """
    An alternative to `functools.wraps` which also preserves the original
    function's call signature by way of
    `~astropy.utils.codegen.make_function_with_signature`.

    This also adds an optional ``exclude_args`` argument.  If given it should
    be a sequence of argument names that should not be copied from the wrapped
    function (either positional or keyword arguments).

    The documentation for the original `functools.wraps` follows:

    """

    wrapped_args = _get_function_args(wrapped, exclude_args=exclude_args)

    def wrapper(func):
        if '__name__' in assigned:
            name = wrapped.__name__
        else:
            name = func.__name__

        func = make_function_with_signature(func, name=name, **wrapped_args)
        func = functools.update_wrapper(func, wrapped, assigned=assigned,
                                        updated=updated)
        return func

    return wrapper


if isinstance(wraps.__doc__, six.string_types):
    wraps.__doc__ += functools.wraps.__doc__


if six.PY3:
    def _get_function_args_internal(func):
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
    def _get_function_args_internal(func):
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


def _get_function_args(func, exclude_args=()):
    all_args = _get_function_args_internal(func)

    if exclude_args:
        exclude_args = set(exclude_args)

        for arg_type in ('args', 'kwargs'):
            all_args[arg_type] = [arg for arg in all_args[arg_type]
                                  if arg not in exclude_args]

        for arg_type in ('varargs', 'varkwargs'):
            if all_args[arg_type] in exclude_args:
                all_args[arg_type] = None

    return all_args
