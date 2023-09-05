# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utilities for generating new Python code at runtime."""


import inspect
import itertools
import keyword
import os
import re
import textwrap

from .introspection import find_current_module

__all__ = ["make_function_with_signature"]


_ARGNAME_RE = re.compile(r"^[A-Za-z][A-Za-z_]*")
"""
Regular expression used my make_func which limits the allowed argument
names for the created function.  Only valid Python variable names in
the ASCII range and not beginning with '_' are allowed, currently.
"""


def make_function_with_signature(
    func, args=(), kwargs={}, varargs=None, varkwargs=None, name=None
):
    """
    Make a new function from an existing function but with the desired
    signature.

    The desired signature must of course be compatible with the arguments
    actually accepted by the input function.

    The ``args`` are strings that should be the names of the positional
    arguments.  ``kwargs`` can map names of keyword arguments to their
    default values.  It may be either a ``dict`` or a list of ``(keyword,
    default)`` tuples.

    If ``varargs`` is a string it is added to the positional arguments as
    ``*<varargs>``.  Likewise ``varkwargs`` can be the name for a variable
    keyword argument placeholder like ``**<varkwargs>``.

    If not specified the name of the new function is taken from the original
    function.  Otherwise, the ``name`` argument can be used to specify a new
    name.

    Note, the names may only be valid Python variable names.
    """
    pos_args = []
    key_args = []

    if isinstance(kwargs, dict):
        iter_kwargs = kwargs.items()
    else:
        iter_kwargs = iter(kwargs)

    # Check that all the argument names are valid
    for item in itertools.chain(args, iter_kwargs):
        if isinstance(item, tuple):
            argname = item[0]
            key_args.append(item)
        else:
            argname = item
            pos_args.append(item)

        if keyword.iskeyword(argname) or not _ARGNAME_RE.match(argname):
            raise SyntaxError(f"invalid argument name: {argname}")

    for item in (varargs, varkwargs):
        if item is not None:
            if keyword.iskeyword(item) or not _ARGNAME_RE.match(item):
                raise SyntaxError(f"invalid argument name: {item}")

    def_signature = [", ".join(pos_args)]

    if varargs:
        def_signature.append(f", *{varargs}")

    call_signature = def_signature[:]

    if name is None:
        name = func.__name__

    global_vars = {f"__{name}__func": func}
    local_vars = {}
    # Make local variables to handle setting the default args
    for idx, item in enumerate(key_args):
        key, value = item
        default_var = f"_kwargs{idx}"
        local_vars[default_var] = value
        def_signature.append(f", {key}={default_var}")
        call_signature.append(f", {key}={key}")

    if varkwargs:
        def_signature.append(f", **{varkwargs}")
        call_signature.append(f", **{varkwargs}")

    def_signature = "".join(def_signature).lstrip(", ")
    call_signature = "".join(call_signature).lstrip(", ")

    mod = find_current_module(2)
    frm = inspect.currentframe().f_back

    if mod:
        filename = mod.__file__
        modname = mod.__name__
        if filename.endswith(".pyc"):
            filename = os.path.splitext(filename)[0] + ".py"
    else:
        filename = "<string>"
        modname = "__main__"

    # Subtract 2 from the line number since the length of the template itself
    # is two lines.  Therefore we have to subtract those off in order for the
    # pointer in tracebacks from __{name}__func to point to the right spot.
    lineno = frm.f_lineno - 2

    # The lstrip is in case there were *no* positional arguments (a rare case)
    # in any context this will actually be used...
    template = textwrap.dedent(
        """{0}\
    def {name}({sig1}):
        return __{name}__func({sig2})
    """.format(
            "\n" * lineno, name=name, sig1=def_signature, sig2=call_signature
        )
    )

    code = compile(template, filename, "single")

    eval(code, global_vars, local_vars)

    new_func = local_vars[name]
    new_func.__module__ = modname
    new_func.__doc__ = func.__doc__

    return new_func
