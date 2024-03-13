# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sundry function and class decorators."""

import functools
import inspect
import textwrap
import threading
import types
import warnings
from inspect import signature

from .exceptions import (
    AstropyDeprecationWarning,
    AstropyPendingDeprecationWarning,
    AstropyUserWarning,
)

__all__ = [
    "classproperty",
    "deprecated",
    "deprecated_attribute",
    "deprecated_renamed_argument",
    "format_doc",
    "lazyproperty",
    "sharedmethod",
]

_NotFound = object()


def deprecated(
    since,
    message="",
    name="",
    alternative="",
    pending=False,
    obj_type=None,
    warning_type=AstropyDeprecationWarning,
):
    """
    Used to mark a function or class as deprecated.

    To mark an attribute as deprecated, use `deprecated_attribute`.

    Parameters
    ----------
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
        ``warning_type``.

    obj_type : str, optional
        The type of this object, if the automatically determined one
        needs to be overridden.

    warning_type : Warning
        Warning to be issued.
        Default is `~astropy.utils.exceptions.AstropyDeprecationWarning`.
    """
    method_types = (classmethod, staticmethod, types.MethodType)

    def deprecate_doc(old_doc, message):
        """
        Returns a given docstring with a deprecation message prepended
        to it.
        """
        if not old_doc:
            old_doc = ""
        old_doc = textwrap.dedent(old_doc).strip("\n")
        new_doc = f"\n.. deprecated:: {since}\n    {message.strip()}\n\n" + old_doc
        if not old_doc:
            # This is to prevent a spurious 'unexpected unindent' warning from
            # docutils when the original docstring was blank.
            new_doc += r"\ "
        return new_doc

    def get_function(func):
        """
        Given a function or classmethod (or other function wrapper type), get
        the function object.
        """
        if isinstance(func, method_types):
            func = func.__func__
        return func

    def deprecate_function(func, message, warning_type=warning_type):
        """
        Returns a wrapped function that displays ``warning_type``
        when it is called.
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
                category = warning_type

            warnings.warn(message, category, stacklevel=2)

            return func(*args, **kwargs)

        # If this is an extension function, we can't call
        # functools.wraps on it, but we normally don't care.
        # This crazy way to get the type of a wrapper descriptor is
        # straight out of the Python 3.3 inspect module docs.
        if type(func) is not type(str.__dict__["__add__"]):
            deprecated_func = functools.wraps(func)(deprecated_func)

        deprecated_func.__doc__ = deprecate_doc(deprecated_func.__doc__, message)
        deprecated_func.__deprecated__ = message

        return func_wrapper(deprecated_func)

    def deprecate_class(cls, message, warning_type=warning_type):
        """
        Update the docstring and wrap the ``__init__`` in-place (or ``__new__``
        if the class or any of the bases overrides ``__new__``) so it will give
        a deprecation warning when an instance is created.

        This won't work for extension classes because these can't be modified
        in-place and the alternatives don't work in the general case:

        - Using a new class that looks and behaves like the original doesn't
          work because the __new__ method of extension types usually makes sure
          that it's the same class or a subclass.
        - Subclassing the class and return the subclass can lead to problems
          with pickle and will look weird in the Sphinx docs.
        """
        cls.__doc__ = deprecate_doc(cls.__doc__, message)
        cls.__deprecated__ = message
        if cls.__new__ is object.__new__:
            cls.__init__ = deprecate_function(
                get_function(cls.__init__), message, warning_type
            )
        else:
            cls.__new__ = deprecate_function(
                get_function(cls.__new__), message, warning_type
            )
        return cls

    def deprecate(
        obj,
        message=message,
        name=name,
        alternative=alternative,
        pending=pending,
        warning_type=warning_type,
    ):
        if obj_type is None:
            if isinstance(obj, type):
                obj_type_name = "class"
            elif inspect.isfunction(obj):
                obj_type_name = "function"
            elif inspect.ismethod(obj) or isinstance(obj, method_types):
                obj_type_name = "method"
            else:
                obj_type_name = "object"
        else:
            obj_type_name = obj_type

        if not name:
            name = get_function(obj).__name__

        altmessage = ""
        if not message or type(message) is type(deprecate):
            if pending:
                message = (
                    "The {func} {obj_type} will be deprecated in a future version."
                )
            else:
                message = (
                    "The {func} {obj_type} is deprecated and may "
                    "be removed in a future version."
                )
            if alternative:
                altmessage = f"\n        Use {alternative} instead."

        message = (
            message.format(
                func=name,
                name=name,
                alternative=alternative,
                obj_type=obj_type_name,
            )
        ) + altmessage

        if isinstance(obj, type):
            return deprecate_class(obj, message, warning_type)
        else:
            return deprecate_function(obj, message, warning_type)

    if type(message) is type(deprecate):
        return deprecate(message)

    return deprecate


def deprecated_attribute(
    name,
    since,
    message=None,
    alternative=None,
    pending=False,
    warning_type=AstropyDeprecationWarning,
):
    """
    Used to mark a public attribute as deprecated.  This creates a
    property that will warn when the given attribute name is accessed.
    To prevent the warning (i.e. for internal code), use the private
    name for the attribute by prepending an underscore
    (i.e. ``self._name``), or set an alternative explicitly.

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
        If True, uses a AstropyPendingDeprecationWarning instead of
        ``warning_type``.

    warning_type : Warning
        Warning to be issued.
        Default is `~astropy.utils.exceptions.AstropyDeprecationWarning`.

    Examples
    --------
    ::

        class MyClass:
            # Mark the old_name as deprecated
            old_name = deprecated_attribute("old_name", "0.1")

            def method(self):
                self._old_name = 42

        class MyClass2:
            old_name = deprecated_attribute(
                "old_name", "1.2", alternative="new_name"
            )

            def method(self):
                self.new_name = 24
    """
    private_name = alternative or "_" + name

    specific_deprecated = deprecated(
        since,
        name=name,
        obj_type="attribute",
        message=message,
        alternative=alternative,
        pending=pending,
        warning_type=warning_type,
    )

    @specific_deprecated
    def get(self):
        return getattr(self, private_name)

    @specific_deprecated
    def set(self, val):
        setattr(self, private_name, val)

    @specific_deprecated
    def delete(self):
        delattr(self, private_name)

    return property(get, set, delete)


def deprecated_renamed_argument(
    old_name,
    new_name,
    since,
    arg_in_kwargs=False,
    relax=False,
    pending=False,
    warning_type=AstropyDeprecationWarning,
    alternative="",
    message="",
):
    """Deprecate a _renamed_ or _removed_ function argument.

    The decorator assumes that the argument with the ``old_name`` was removed
    from the function signature and the ``new_name`` replaced it at the
    **same position** in the signature.  If the ``old_name`` argument is
    given when calling the decorated function the decorator will catch it and
    issue a deprecation warning and pass it on as ``new_name`` argument.

    Parameters
    ----------
    old_name : str or sequence of str
        The old name of the argument.

    new_name : str or sequence of str or None
        The new name of the argument. Set this to `None` to remove the
        argument ``old_name`` instead of renaming it.

    since : str or number or sequence of str or number
        The release at which the old argument became deprecated.

    arg_in_kwargs : bool or sequence of bool, optional
        If the argument is not a named argument (for example it
        was meant to be consumed by ``**kwargs``) set this to
        ``True``.  Otherwise the decorator will throw an Exception
        if the ``new_name`` cannot be found in the signature of
        the decorated function.
        Default is ``False``.

    relax : bool or sequence of bool, optional
        If ``False`` a ``TypeError`` is raised if both ``new_name`` and
        ``old_name`` are given.  If ``True`` the value for ``new_name`` is used
        and a Warning is issued.
        Default is ``False``.

    pending : bool or sequence of bool, optional
        If ``True`` this will hide the deprecation warning and ignore the
        corresponding ``relax`` parameter value.
        Default is ``False``.

    warning_type : Warning
        Warning to be issued.
        Default is `~astropy.utils.exceptions.AstropyDeprecationWarning`.

    alternative : str, optional
        An alternative function or class name that the user may use in
        place of the deprecated object if ``new_name`` is None. The deprecation
        warning will tell the user about this alternative if provided.

    message : str, optional
        A custom warning message. If provided then ``since`` and
        ``alternative`` options will have no effect.

    Raises
    ------
    TypeError
        If the new argument name cannot be found in the function
        signature and arg_in_kwargs was False or if it is used to
        deprecate the name of the ``*args``-, ``**kwargs``-like arguments.
        At runtime such an Error is raised if both the new_name
        and old_name were specified when calling the function and
        "relax=False".

    Notes
    -----
    The decorator should be applied to a function where the **name**
    of an argument was changed but it applies the same logic.

    .. warning::
        If ``old_name`` is a list or tuple the ``new_name`` and ``since`` must
        also be a list or tuple with the same number of entries. ``relax`` and
        ``arg_in_kwarg`` can be a single bool (applied to all) or also a
        list/tuple with the same number of entries like ``new_name``, etc.

    Examples
    --------
    The deprecation warnings are not shown in the following examples.

    To deprecate a positional or keyword argument::

        >>> from astropy.utils.decorators import deprecated_renamed_argument
        >>> @deprecated_renamed_argument('sig', 'sigma', '1.0')
        ... def test(sigma):
        ...     return sigma

        >>> test(2)
        2
        >>> test(sigma=2)
        2
        >>> test(sig=2)  # doctest: +SKIP
        2

    To deprecate an argument caught inside the ``**kwargs`` the
    ``arg_in_kwargs`` has to be set::

        >>> @deprecated_renamed_argument('sig', 'sigma', '1.0',
        ...                             arg_in_kwargs=True)
        ... def test(**kwargs):
        ...     return kwargs['sigma']

        >>> test(sigma=2)
        2
        >>> test(sig=2)  # doctest: +SKIP
        2

    By default providing the new and old keyword will lead to an Exception. If
    a Warning is desired set the ``relax`` argument::

        >>> @deprecated_renamed_argument('sig', 'sigma', '1.0', relax=True)
        ... def test(sigma):
        ...     return sigma

        >>> test(sig=2)  # doctest: +SKIP
        2

    It is also possible to replace multiple arguments. The ``old_name``,
    ``new_name`` and ``since`` have to be `tuple` or `list` and contain the
    same number of entries::

        >>> @deprecated_renamed_argument(['a', 'b'], ['alpha', 'beta'],
        ...                              ['1.0', 1.2])
        ... def test(alpha, beta):
        ...     return alpha, beta

        >>> test(a=2, b=3)  # doctest: +SKIP
        (2, 3)

    In this case ``arg_in_kwargs`` and ``relax`` can be a single value (which
    is applied to all renamed arguments) or must also be a `tuple` or `list`
    with values for each of the arguments.

    """
    cls_iter = (list, tuple)
    if isinstance(old_name, cls_iter):
        n = len(old_name)
        # Assume that new_name and since are correct (tuple/list with the
        # appropriate length) in the spirit of the "consenting adults". But the
        # optional parameters may not be set, so if these are not iterables
        # wrap them.
        if not isinstance(arg_in_kwargs, cls_iter):
            arg_in_kwargs = [arg_in_kwargs] * n
        if not isinstance(relax, cls_iter):
            relax = [relax] * n
        if not isinstance(pending, cls_iter):
            pending = [pending] * n
        if not isinstance(message, cls_iter):
            message = [message] * n
    else:
        # To allow a uniform approach later on, wrap all arguments in lists.
        n = 1
        old_name = [old_name]
        new_name = [new_name]
        since = [since]
        arg_in_kwargs = [arg_in_kwargs]
        relax = [relax]
        pending = [pending]
        message = [message]

    def decorator(function):
        # The named arguments of the function.
        arguments = signature(function).parameters
        keys = list(arguments.keys())
        position = [None] * n

        for i in range(n):
            # Determine the position of the argument.
            if arg_in_kwargs[i]:
                pass
            else:
                if new_name[i] is None:
                    param = arguments[old_name[i]]
                elif new_name[i] in arguments:
                    param = arguments[new_name[i]]
                # In case the argument is not found in the list of arguments
                # the only remaining possibility is that it should be caught
                # by some kind of **kwargs argument.
                # This case has to be explicitly specified, otherwise throw
                # an exception!
                else:
                    raise TypeError(
                        f'"{new_name[i]}" was not specified in the function '
                        "signature. If it was meant to be part of "
                        '"**kwargs" then set "arg_in_kwargs" to "True"'
                    )

                # There are several possibilities now:

                # 1.) Positional or keyword argument:
                if param.kind == param.POSITIONAL_OR_KEYWORD:
                    if new_name[i] is None:
                        position[i] = keys.index(old_name[i])
                    else:
                        position[i] = keys.index(new_name[i])

                # 2.) Keyword only argument:
                elif param.kind == param.KEYWORD_ONLY:
                    # These cannot be specified by position.
                    position[i] = None

                # 3.) positional-only argument, varargs, varkwargs or some
                #     unknown type:
                else:
                    raise TypeError(
                        f'cannot replace argument "{new_name[i]}" '
                        f"of kind {repr(param.kind)}."
                    )

        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            for i in range(n):
                msg = message[i] or (
                    f'"{old_name[i]}" was deprecated in '
                    f"version {since[i]} and will be removed "
                    "in a future version. "
                )
                # The only way to have oldkeyword inside the function is
                # that it is passed as kwarg because the oldkeyword
                # parameter was renamed to newkeyword.
                if old_name[i] in kwargs:
                    value = kwargs.pop(old_name[i])
                    # Display the deprecation warning only when it's not
                    # pending.
                    if not pending[i]:
                        if not message[i]:
                            if new_name[i] is not None:
                                msg += f'Use argument "{new_name[i]}" instead.'
                            elif alternative:
                                msg += f"\n        Use {alternative} instead."
                        warnings.warn(msg, warning_type, stacklevel=2)

                    # Check if the newkeyword was given as well.
                    newarg_in_args = position[i] is not None and len(args) > position[i]
                    newarg_in_kwargs = new_name[i] in kwargs

                    if newarg_in_args or newarg_in_kwargs:
                        if not pending[i]:
                            # If both are given print a Warning if relax is
                            # True or raise an Exception is relax is False.
                            if relax[i]:
                                warnings.warn(
                                    f'"{old_name[i]}" and "{new_name[i]}" '
                                    "keywords were set. "
                                    f'Using the value of "{new_name[i]}".',
                                    AstropyUserWarning,
                                )
                            else:
                                raise TypeError(
                                    f'cannot specify both "{old_name[i]}" and '
                                    f'"{new_name[i]}".'
                                )
                    else:
                        # Pass the value of the old argument with the
                        # name of the new argument to the function
                        if new_name[i] is not None:
                            kwargs[new_name[i]] = value
                        # If old argument has no replacement, cast it back.
                        # https://github.com/astropy/astropy/issues/9914
                        else:
                            kwargs[old_name[i]] = value

                # Deprecated keyword without replacement is given as
                # positional argument.
                elif (
                    not pending[i]
                    and not new_name[i]
                    and position[i]
                    and len(args) > position[i]
                ):
                    if alternative and not message[i]:
                        msg += f"\n        Use {alternative} instead."
                    warnings.warn(msg, warning_type, stacklevel=2)

            return function(*args, **kwargs)

        return wrapper

    return decorator


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
        currently allow writeable/deletable properties, due to subtleties of how
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

        >>> class Foo:
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

        >>> class Foo:
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

        >>> class Foo:
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

        return super().__new__(cls)

    def __init__(self, fget, doc=None, lazy=False):
        self._lazy = lazy
        if lazy:
            self._lock = threading.RLock()  # Protects _cache
            self._cache = {}
        fget = self._wrap_fget(fget)

        super().__init__(fget=fget, doc=doc)

        # There is a buglet in Python where self.__doc__ doesn't
        # get set properly on instances of property subclasses if
        # the doc argument was used rather than taking the docstring
        # from fget
        # Related Python issue: https://bugs.python.org/issue24766
        if doc is not None:
            self.__doc__ = doc

    def __get__(self, obj, objtype):
        if self._lazy:
            val = self._cache.get(objtype, _NotFound)
            if val is _NotFound:
                with self._lock:
                    # Check if another thread initialised before we locked.
                    val = self._cache.get(objtype, _NotFound)
                    if val is _NotFound:
                        val = self.fget.__wrapped__(objtype)
                        self._cache[objtype] = val
        else:
            # The base property.__get__ will just return self here;
            # instead we pass objtype through to the original wrapped
            # function (which takes the class as its sole argument)
            val = self.fget.__wrapped__(objtype)
        return val

    def getter(self, fget):
        return super().getter(self._wrap_fget(fget))

    def setter(self, fset):
        raise NotImplementedError(
            "classproperty can only be read-only; use a metaclass to "
            "implement modifiable class-level properties"
        )

    def deleter(self, fdel):
        raise NotImplementedError(
            "classproperty can only be read-only; use a metaclass to "
            "implement modifiable class-level properties"
        )

    @staticmethod
    def _wrap_fget(orig_fget):
        if isinstance(orig_fget, classmethod):
            orig_fget = orig_fget.__func__

        # Using stock functools.wraps instead of the fancier version
        # found later in this module, which is overkill for this purpose

        @functools.wraps(orig_fget)
        def fget(obj):
            return orig_fget(obj.__class__)

        return fget


# Adapted from the recipe at
# http://code.activestate.com/recipes/363602-lazy-property-evaluation
class lazyproperty(property):
    """
    Works similarly to property(), but computes the value only once.

    This essentially memorizes the value of the property by storing the result
    of its computation in the ``__dict__`` of the object instance.  This is
    useful for computing the value of some property that should otherwise be
    invariant.  For example::

        >>> class LazyTest:
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

    By default, a setter and deleter are used which simply overwrite and
    delete, respectively, the value stored in ``__dict__``. Any user-specified
    setter or deleter is executed before executing these default actions.
    The one exception is that the default setter is not run if the user setter
    already sets the new value in ``__dict__`` and returns that value and the
    returned value is not ``None``.

    """

    def __init__(self, fget, fset=None, fdel=None, doc=None):
        super().__init__(fget, fset, fdel, doc)
        self._key = self.fget.__name__
        self._lock = threading.RLock()

    def __get__(self, obj, owner=None):
        try:
            obj_dict = obj.__dict__
            val = obj_dict.get(self._key, _NotFound)
            if val is _NotFound:
                with self._lock:
                    # Check if another thread beat us to it.
                    val = obj_dict.get(self._key, _NotFound)
                    if val is _NotFound:
                        val = self.fget(obj)
                        obj_dict[self._key] = val
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
                # By returning the value set the setter signals that it
                # took over setting the value in obj.__dict__; this
                # mechanism allows it to override the input value
                return
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        obj.__dict__.pop(self._key, None)  # Delete if present


class sharedmethod(classmethod):
    """
    This is a method decorator that allows both an instance method and a
    `classmethod` to share the same name.

    When using `sharedmethod` on a method defined in a class's body, it
    may be called on an instance, or on a class.  In the former case it
    behaves like a normal instance method (a reference to the instance is
    automatically passed as the first ``self`` argument of the method)::

        >>> class Example:
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
    implementation can be written separately.  If the class' *metaclass*
    has a method of the same name as the `sharedmethod`, the version on
    the metaclass is delegated to::

        >>> class ExampleMeta(type):
        ...     def identify(self):
        ...         print('this implements the {0}.identify '
        ...               'classmethod'.format(self.__name__))
        ...
        >>> class Example(metaclass=ExampleMeta):
        ...     @sharedmethod
        ...     def identify(self):
        ...         print('this implements the instancemethod')
        ...
        >>> Example().identify()
        this implements the instancemethod
        >>> Example.identify()
        this implements the Example.identify classmethod
    """

    def __get__(self, obj, objtype=None):
        if obj is None:
            mcls = type(objtype)
            clsmeth = getattr(mcls, self.__func__.__name__, None)
            if callable(clsmeth):
                func = clsmeth
            else:
                func = self.__func__

            return self._make_method(func, objtype)
        else:
            return self._make_method(self.__func__, obj)

    @staticmethod
    def _make_method(func, instance):
        return types.MethodType(func, instance)


def format_doc(docstring, *args, **kwargs):
    """
    Replaces the docstring of the decorated object and then formats it.

    The formatting works like :meth:`str.format` and if the decorated object
    already has a docstring this docstring can be included in the new
    documentation if you use the ``{__doc__}`` placeholder.
    Its primary use is for reusing a *long* docstring in multiple functions
    when it is the same or only slightly different between them.

    Parameters
    ----------
    docstring : str or object or None
        The docstring that will replace the docstring of the decorated
        object. If it is an object like a function or class it will
        take the docstring of this object. If it is a string it will use the
        string itself. One special case is if the string is ``None`` then
        it will use the decorated functions docstring and formats it.

    args :
        passed to :meth:`str.format`.

    kwargs :
        passed to :meth:`str.format`. If the function has a (not empty)
        docstring the original docstring is added to the kwargs with the
        keyword ``'__doc__'``.

    Raises
    ------
    ValueError
        If the ``docstring`` (or interpreted docstring if it was ``None``
        or not a string) is empty.

    IndexError, KeyError
        If a placeholder in the (interpreted) ``docstring`` was not filled. see
        :meth:`str.format` for more information.

    Notes
    -----
    Using this decorator allows, for example Sphinx, to parse the
    correct docstring.

    Examples
    --------
    Replacing the current docstring is very easy::

        >>> from astropy.utils.decorators import format_doc
        >>> @format_doc('''Perform num1 + num2''')
        ... def add(num1, num2):
        ...     return num1+num2
        ...
        >>> help(add) # doctest: +SKIP
        Help on function add in module __main__:
        <BLANKLINE>
        add(num1, num2)
            Perform num1 + num2

    sometimes instead of replacing you only want to add to it::

        >>> doc = '''
        ...       {__doc__}
        ...       Parameters
        ...       ----------
        ...       num1, num2 : Numbers
        ...       Returns
        ...       -------
        ...       result: Number
        ...       '''
        >>> @format_doc(doc)
        ... def add(num1, num2):
        ...     '''Perform addition.'''
        ...     return num1+num2
        ...
        >>> help(add) # doctest: +SKIP
        Help on function add in module __main__:
        <BLANKLINE>
        add(num1, num2)
            Perform addition.
            Parameters
            ----------
            num1, num2 : Numbers
            Returns
            -------
            result : Number

    in case one might want to format it further::

        >>> doc = '''
        ...       Perform {0}.
        ...       Parameters
        ...       ----------
        ...       num1, num2 : Numbers
        ...       Returns
        ...       -------
        ...       result: Number
        ...           result of num1 {op} num2
        ...       {__doc__}
        ...       '''
        >>> @format_doc(doc, 'addition', op='+')
        ... def add(num1, num2):
        ...     return num1+num2
        ...
        >>> @format_doc(doc, 'subtraction', op='-')
        ... def subtract(num1, num2):
        ...     '''Notes: This one has additional notes.'''
        ...     return num1-num2
        ...
        >>> help(add) # doctest: +SKIP
        Help on function add in module __main__:
        <BLANKLINE>
        add(num1, num2)
            Perform addition.
            Parameters
            ----------
            num1, num2 : Numbers
            Returns
            -------
            result : Number
                result of num1 + num2
        >>> help(subtract) # doctest: +SKIP
        Help on function subtract in module __main__:
        <BLANKLINE>
        subtract(num1, num2)
            Perform subtraction.
            Parameters
            ----------
            num1, num2 : Numbers
            Returns
            -------
            result : Number
                result of num1 - num2
            Notes : This one has additional notes.

    These methods can be combined; even taking the docstring from another
    object is possible as docstring attribute. You just have to specify the
    object::

        >>> @format_doc(add)
        ... def another_add(num1, num2):
        ...     return num1 + num2
        ...
        >>> help(another_add) # doctest: +SKIP
        Help on function another_add in module __main__:
        <BLANKLINE>
        another_add(num1, num2)
            Perform addition.
            Parameters
            ----------
            num1, num2 : Numbers
            Returns
            -------
            result : Number
                result of num1 + num2

    But be aware that this decorator *only* formats the given docstring not
    the strings passed as ``args`` or ``kwargs`` (not even the original
    docstring)::

        >>> @format_doc(doc, 'addition', op='+')
        ... def yet_another_add(num1, num2):
        ...    '''This one is good for {0}.'''
        ...    return num1 + num2
        ...
        >>> help(yet_another_add) # doctest: +SKIP
        Help on function yet_another_add in module __main__:
        <BLANKLINE>
        yet_another_add(num1, num2)
            Perform addition.
            Parameters
            ----------
            num1, num2 : Numbers
            Returns
            -------
            result : Number
                result of num1 + num2
            This one is good for {0}.

    To work around it you could specify the docstring to be ``None``::

        >>> @format_doc(None, 'addition')
        ... def last_add_i_swear(num1, num2):
        ...    '''This one is good for {0}.'''
        ...    return num1 + num2
        ...
        >>> help(last_add_i_swear) # doctest: +SKIP
        Help on function last_add_i_swear in module __main__:
        <BLANKLINE>
        last_add_i_swear(num1, num2)
            This one is good for addition.

    Using it with ``None`` as docstring allows to use the decorator twice
    on an object to first parse the new docstring and then to parse the
    original docstring or the ``args`` and ``kwargs``.
    """

    def set_docstring(obj):
        if docstring is None:
            # None means: use the objects __doc__
            doc = obj.__doc__
            # Delete documentation in this case so we don't end up with
            # awkwardly self-inserted docs.
            obj.__doc__ = None
        elif isinstance(docstring, str):
            # String: use the string that was given
            doc = docstring
        else:
            # Something else: Use the __doc__ of this
            doc = docstring.__doc__

        if not doc:
            # In case the docstring is empty it's probably not what was wanted.
            raise ValueError(
                "docstring must be a string or containing a "
                "docstring that is not empty."
            )

        # If the original has a not-empty docstring append it to the format
        # kwargs.
        kwargs["__doc__"] = obj.__doc__ or ""
        obj.__doc__ = doc.format(*args, **kwargs)
        return obj

    return set_docstring
