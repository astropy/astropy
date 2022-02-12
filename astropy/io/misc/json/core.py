# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import inspect
import json

import numpy as np

import astropy.units as u
import astropy.coordinates as coord
from astropy.utils.decorators import format_doc
from astropy.utils.misc import indent

from . import __doc__ as _example_doc

__all__ = ['JSONExtendedEncoder', 'JSONExtendedDecoder']


QUALNAME_SUBSTITUTIONS = {}


class JSONExtendedDecodeError(json.JSONDecodeError):
    """Decoding error for JSON-Extended."""


@format_doc(None, examples=indent(_example_doc[77:])[4:])  # (cut off 1st line)
class JSONExtendedEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder does not do.

    This includes (but is not limited to):

    - :mod:`builtins`: `bytes`, `complex`, `set`, `NotImplemented`, `Ellipsis`
    - :mod:`numpy`: `~numpy.array`, `~numpy.number`, `~numpy.dtype`
    - Astropy:

        - :mod:`astropy.units`: `~astropy.units.Unit`, `~astropy.units.Quantity`
        - :mod:`astropy.coordinates`: `~astropy.coordinates.Longitude`
        - :mod:`~astropy.cosmology`: `~astropy.cosmology.Cosmology`

    Parameters
    ----------
    skipkeys : bool, optional keyword-only
        See `json.JSONEncoder` for details.
    ensure_ascii : bool, optional keyword-only
        See `json.JSONEncoder` for details.
    check_circular : bool, optional keyword-only
        See `json.JSONEncoder` for details.
    allow_nan : bool, optional keyword-only
        See `json.JSONEncoder` for details.
    sort_keys : bool, optional keyword-only
        See `json.JSONEncoder` for details.
    indent : int or str or None, optional keyword-only
        See `json.JSONEncoder` for details.
    separators : tuple or None, optional keyword-only
        See `json.JSONEncoder` for details.
    default : None, optional keyword-only
        Must be `None`. A `ValueError` is raised otherwise.

    Raises
    ------
    ValueError
        If ``default`` is not `None`.

    Examples
    --------
    {examples}
    """

    _registry = []
    # _registry = np.array([("", object())],
    #                      dtype=np.dtype([("key", "U110"), ("func", object)], align=True))

    def __init__(self, *, skipkeys=False, ensure_ascii=True, check_circular=True,
                 allow_nan=True, sort_keys=False, indent=None, separators=None,
                 default=None):
        if default is not None:
            raise ValidateError("`default` must be None for `JSONExtendedEncoder`.")
        super().__init__(skipkeys=skipkeys, ensure_ascii=ensure_ascii,
                         check_circular=check_circular, allow_nan=allow_nan,
                         sort_keys=sort_keys, indent=indent,
                         separators=separators)

    def default(self, obj):
        """
        Returns a serializable object for ``obj``, or calls the base
        implementation (to raise a `TypeError`).
        """
        for cls, encoder_func in self._registry:
            if isinstance(obj, cls):
                code = encoder_func(obj)
                break
        else:  # Calls the base implementation.
            code = super().default(obj)

        return code

    @classmethod
    def register_encoding(cls, type):
        """Return a function for registering an encoding. Can be used as a decorator."""
        # TODO! check if already registered `type`

        def register(func):
            # inserting subclasses before parent classes, so encountered first
            # when finding the right encoder.
            for i, (key, _) in enumerate(cls._registry):
                if issubclass(type, key):
                    cls._registry.insert(i, (type, func))
                    break
            else:  # put at the end
                cls._registry.append((type, func))
            return func

        return register


@format_doc(None, examples=indent(_example_doc[77:])[4:])  # (cut off 1st line)
class JSONExtendedDecoder(json.JSONDecoder):
    """Support for data types that JSON default decoder does not do.

    This includes (but is not limited to):

    - :mod:`builtins`: `bytes`, `complex`, `set`, `NotImplemented`, `Ellipsis`
    - :mod:`numpy`: `~numpy.array`, `~numpy.number`, `~numpy.dtype`
    - Astropy:

        - :mod:`astropy.units`: `~astropy.units.Unit`, `~astropy.units.Quantity`
        - :mod:`astropy.coordinates`: `~astropy.coordinates.Longitude`
        - :mod:`~astropy.cosmology`: `~astropy.cosmology.Cosmology`

    Parameters
    ----------
    parse_float : callable or None, optional keyword-only
        See ``json.JSONDecoder`` for details.
    parse_int : callable or None, optional keyword-only
        See ``json.JSONDecoder`` for details.
    parse_constant : callable or None, optional keyword-only
        See ``json.JSONDecoder`` for details.
    strict : bool, optional keyword-only
        See ``json.JSONDecoder`` for details.

    Examples
    --------
    {examples}
    """

    _registry = []

    def __init__(self, *, parse_float=None, parse_int=None, parse_constant=None,
                 strict=True):
        super().__init__(object_hook=self.object_hook, parse_float=parse_float,
                         parse_int=parse_int, parse_constant=parse_constant,
                         strict=strict)

    @classmethod
    def object_hook(cls, code):
        """Called with the result of every JSON object decoded.

        If a method exists in the registry, the return value will be used in
        place of the given dict, otherwise the dict is returned unchanged.
        """
        try:  # Import the constructor from the type field ("!")
            qualname = code["!"].split(".")
            module = importlib.import_module(".".join(qualname[:-1]))
            constructor = getattr(module, qualname[-1])
        except ModuleNotFoundError as e:
            raise JSONExtendedDecodeError from e
            return code

        # Iterate through the decoder registry, trying to figure out if there
        # is a way to decode the object. If not, a warning is issued.
        for key, func in cls._registry:  # try to decode
            if isinstance(constructor, key) or (inspect.isclass(constructor) and issubclass(constructor, key)):
                code.pop("!")
                obj = func(constructor, code.pop("value"), code)
                break
        else:  # Looks like a valid JSONExtended, but it is NOT.
            obj = code

        return obj

    @classmethod
    def register_decoding(cls, type):
        """Return a function for registering a decoding. Can be used as a decorator."""
        def register(func):
            # Subclasses are inserted before parent classes, so the more
            # specific subclass is encountered first when finding a decoder.
            for i, (key, _) in enumerate(cls._registry):
                if issubclass(type, key):
                    cls._registry.insert(i, (type, func))
                    break
            else:  # put at the end
                cls._registry.append((type, func))
            return func

        return register


def _json_base_encode(obj):
    qualname = obj.__class__.__module__ + "." + obj.__class__.__qualname__
    qualname = QUALNAME_SUBSTITUTIONS.get(qualname, qualname)  # possibly substitute
    code = {"!": qualname}
    return code
