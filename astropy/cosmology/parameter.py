# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com), Roban Kramer
# (robanhk@gmail.com), and Nathaniel Starkman (n.starkman@mail.utoronto.ca).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["Parameter"]

# registry of cosmology classes with {key=name : value=class}
_COSMOLOGY_CLASSES = dict()


class Parameter(property):
    r"""Cosmological parameter (descriptor).

    Should only be used with a :class:`~astropy.cosmology.Cosmology` subclass.

    Parameters
    ----------
    fget : callable[[type], Any] or None, optional
        Function to get the value from instances of the cosmology class.
        If None (default) returns the corresponding private attribute.
        Often not set here, but as a decorator with ``getter``.
    fset : callable[[object, object, Any], Any] or {'default', 'float'}, optional
        Function to validate the Parameter value from instances of the
        cosmology class. If "default", uses default setter to assign units
        (with equivalencies), if Parameter has units. If "float" will first do
        units, then take the float value.
        Often not set here, but as a decorator with ``setter``.
    doc : str or None, optional
        Parameter description. If 'doc' is None and 'fget' is not, then 'doc'
        is taken from ``fget.__doc__``.
    unit : unit-like or None (optional, keyword-only)
        The `~astropy.units.Unit` for the Parameter. If None (default) no
        unit as assumed.
    equivalencies : `~astropy.units.Equivalency` or sequence thereof
        Unit equivalencies for this Parameter.
    fmt : str (optional, keyword-only)
        `format` specification, used when making string representation
        of the containing Cosmology.
        See https://docs.python.org/3/library/string.html#formatspec

    derived : bool (optional, keyword-only)
        Whether the Parameter is 'derived', default `False`.
        Fixed parameters behave similarly to normal parameters, but are not
        sorted by the |Cosmology| signature (probably not there) and are not
        included in all methods. For reference, see ``Ode0`` in
        ``FlatFLRWMixin``, which removes :math:`\Omega_{de,0}`` as an
        independent parameter (:math:`\Omega_{de,0} \equiv 1 - \Omega_{tot}`).

    Examples
    --------
    The most common use case of ``Parameter`` is to access the corresponding
    private attribute.

        >>> from astropy.cosmology import LambdaCDM
        >>> from astropy.cosmology.core import Parameter
        >>> class Example1(LambdaCDM):
        ...     param = Parameter(doc="example parameter", unit=u.m)
        ...     def __init__(self, param=15 * u.m):
        ...         super().__init__(70, 0.3, 0.7)
        ...         self._param = param << self.__class__.param.unit
        >>> Example1.param
        <Parameter 'param' at ...
        >>> Example1.param.unit
        Unit("m")

        >>> ex = Example1(param=12357)
        >>> ex.param
        <Quantity 12357. m>

    ``Parameter`` also supports custom ``setter`` methods.
    :attr:`~astropy.cosmology.FLRW.m_nu` is a good example for the former.

        >>> import astropy.units as u
        >>> class Example2(LambdaCDM):
        ...     param = Parameter(doc="example parameter", unit="m")
        ...     def __init__(self, param=15):
        ...         super().__init__(70, 0.3, 0.7)
        ...         self.param = param
        ...     @param.setter
        ...     def param(self, param, value):
        ...         return (value << param.unit).to(u.km)

        >>> ex2 = Example2(param=12357)
        >>> ex2.param
        <Quantity 12.357 km>

    .. doctest::
       :hide:

       >>> from astropy.cosmology.core import _COSMOLOGY_CLASSES
       >>> _ = _COSMOLOGY_CLASSES.pop(Example1.__qualname__)
       >>> _ = _COSMOLOGY_CLASSES.pop(Example2.__qualname__)
    """

    _registry_setters = {}

    def __init__(self, fget=None, fset="default", doc=None, *,
                 unit=None, equivalencies=[], fmt=".3g", derived=False):
        # parse registered fset
        if callable(fset):
            pass
        elif fset in self._registry_setters:
            fset = self._registry_setters[fset]
        elif isinstance(fset, str):
            raise ValueError(f"`fset` if str, must be in {self._registry_setters.keys()}")
        else:
            raise TypeError(f"`fset` must be a function or {self._registry_setters.keys()}")

        # modeled after https://docs.python.org/3/howto/descriptor.html#properties
        super().__init__(fget=fget if not hasattr(fget, "fget") else fget.__get__,
                         fset=fset)
        # TODO! better detection if `fget` is a descriptor.
        # Note: setting here b/c @propert(doc=) is broken in subclasses
        self.__doc__ = fget.__doc__ if (doc is None and fget is not None) else doc

        # units stuff
        self._unit = u.Unit(unit) if unit is not None else None
        self._equivalencies = equivalencies

        # misc
        self._fmt = str(fmt)
        self._derived = derived

        # nested descriptor decorator compatibility
        self.__wrapped__ = fget  # so always have access to `fget`
        self.__name__ = getattr(fget, "__name__", None)  # compat with other descriptors

    def __set_name__(self, cosmo_cls, name):
        # attribute name
        self._attr_name = name
        self._attr_name_private = "_" + name

        # update __name__, if not already set
        self.__name__ = self.__name__ or name

    @property
    def name(self):
        """Parameter name."""
        return self._attr_name

    @property
    def unit(self):
        """Parameter unit."""
        return self._unit

    @property
    def equivalencies(self):
        """Equivalencies used when initializing Parameter."""
        return self._equivalencies

    @property
    def format_spec(self):
        """String format specification."""
        return self._fmt

    @property
    def derived(self):
        """Whether the Parameter is derived; true parameters are not."""
        return self._derived

    # -------------------------------------------
    # descriptor methods

    def __get__(self, cosmology, cosmo_cls=None):
        # get from class
        if cosmology is None:
            return self
        return getattr(cosmology, self._attr_name_private)

    def __set__(self, cosmology, value):
        """Allows attribute setting once. Raises AttributeError subsequently."""
        # raise error if setting 2nd time.
        if hasattr(cosmology, self._attr_name_private):
            raise AttributeError("can't set attribute")

        # validate value, generally setting units if present
        value = self.set(cosmology, value)
        setattr(cosmology, self._attr_name_private, value)

    # -------------------------------------------
    # 'property' descriptor overrides

    def getter(self, fget):
        raise AttributeError("can't create custom Parameter getter.")

    def setter(self, fset):
        """Make new Parameter with custom ``fset``.

        Note: ``Parameter.setter`` must be the top-most descriptor decorator.

        Parameters
        ----------
        fset : callable[[type, type, Any], Any]

        Returns
        -------
        `~astropy.cosmology.Parameter`
            Copy of this Parameter but with custom ``fset``.
        """
        desc = type(self)(fget=self.fget, fset=fset,
                          doc=self.__doc__, fmt=self.format_spec,
                          unit=self.unit, equivalencies=self.equivalencies,
                          derived=self.derived)
        # TODO? need to override __wrapped__?
        return desc

    def deleter(self, fdel):
        raise AttributeError("can't create custom Parameter deleter.")

    # -------------------------------------------
    # set value

    def set(self, cosmology, value):
        """Run the setter on this Parameter.

        Note this setter doesn't actually set the value, but returns the value
        that *should* be set.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` instance
        value : Any
            The object to validate and set.

        Returns
        -------
        Any
            The output of calling ``fset(cosmology, self, value)``
            (yes, that parameter order).
        """
        return self.fset(cosmology, self, value)

    @classmethod
    def register_setter(cls, key, setter=None):
        """Decorator to register setter function.

        Parameters
        ----------
        key : str
        setter : callable[[object, object, Any], Any] or None, optional
            Value setter / validation function.

        Returns
        -------
        ``setter`` or callable[``setter``]
            if setter is None returns a function that takes and registers a
            setter. This allows ``register_setter`` to be used as a decorator.
        """
        if key in cls._registry_setters:
            raise KeyError(f"setter {key!r} already registered with Parameter.")

        # setter directly passed
        if setter is not None:
            cls._registry_setters[key] = setter
            return setter

        # for use as a decorator
        def register(setter):
            """Register setter function.

            Parameters
            ----------
            setter : callable[[object, object, Any], Any]
                Value setter / validation function.

            Returns
            -------
            ``setter``
            """
            cls._registry_setters[key] = setter
            return setter

        return register

    # -------------------------------------------

    def __repr__(self):
        return f"<Parameter {self._attr_name!r} at {hex(id(self))}>"


# ===================================================================
# Built-in setters


@Parameter.register_setter("default")
def _set_with_unit(cosmology, param, value):
    """
    Default Parameter value setter.
    Adds/converts units if Parameter has a unit.
    """
    if param.unit is not None:
        with u.add_enabled_equivalencies(param.equivalencies):
            value = u.Quantity(value, param.unit)
    return value


@Parameter.register_setter("float")
def _set_to_float(cosmology, param, value):
    """Parameter value setter with units, and converted to float."""
    value = _set_with_unit(cosmology, param, value)
    return float(value)


@Parameter.register_setter("scalar")
def _set_to_scalar(cosmology, param, value):
    """"""
    value = _set_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@Parameter.register_setter("non-negative")
def _set_non_negative(cosmology, param, value):
    """Parameter value setter where value is a positive float."""
    value = _set_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} can not be negative.")
    return value
