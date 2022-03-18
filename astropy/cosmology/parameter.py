# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u

__all__ = ["Parameter"]


class Parameter:
    r"""Cosmological parameter (descriptor).

    Should only be used with a :class:`~astropy.cosmology.Cosmology` subclass.

    Parameters
    ----------
    derived : bool (optional, keyword-only)
        Whether the Parameter is 'derived', default `False`.
        Derived parameters behave similarly to normal parameters, but are not
        sorted by the |Cosmology| signature (probably not there) and are not
        included in all methods. For reference, see ``Ode0`` in
        ``FlatFLRWMixin``, which removes :math:`\Omega_{de,0}`` as an
        independent parameter (:math:`\Omega_{de,0} \equiv 1 - \Omega_{tot}`).
    unit : unit-like or None (optional, keyword-only)
        The `~astropy.units.Unit` for the Parameter. If None (default) no
        unit as assumed.
    equivalencies : `~astropy.units.Equivalency` or sequence thereof
        Unit equivalencies for this Parameter.
    fvalidate : callable[[object, object, Any], Any] or str (optional, keyword-only)
        Function to validate the Parameter value from instances of the
        cosmology class. If "default", uses default validator to assign units
        (with equivalencies), if Parameter has units.
        For other valid string options, see ``Parameter._registry_validators``.
        'fvalidate' can also be set through a decorator with
        :meth:`~astropy.cosmology.Parameter.validator`.
    fmt : str (optional, keyword-only)
        `format` specification, used when making string representation
        of the containing Cosmology.
        See https://docs.python.org/3/library/string.html#formatspec
    doc : str or None (optional, keyword-only)
        Parameter description.

    Examples
    --------
    For worked examples see :class:`~astropy.cosmology.FLRW`.
    """

    _registry_validators = {}

    def __init__(self, *, derived=False, unit=None, equivalencies=[],
                 fvalidate="default", fmt="", doc=None):

        # attribute name on container cosmology class.
        # really set in __set_name__, but if Parameter is not init'ed as a
        # descriptor this ensures that the attributes exist.
        self._attr_name = self._attr_name_private = None

        self._derived = derived
        self._fmt = str(fmt)  # @property is `format_spec`
        self.__doc__ = doc

        # units stuff
        # can be a unit-able, None, or a structured unit of unknown length.
        # the latter is looks like (unit-able, ...)
        # See `_validate_with_unit` for usage.
        if unit is None:
            self._unit = unit
        elif isinstance(unit, tuple) and (len(unit) == 2) and (unit[1] is Ellipsis):
            self._unit = (u.Unit(unit[0]), ...)
        else:  # this also catches errors
            self._unit = u.Unit(unit)
        self._equivalencies = equivalencies

        # Parse registered `fvalidate`
        self._fvalidate_in = fvalidate  # Always store input fvalidate.
        if callable(fvalidate):
            pass
        elif fvalidate in self._registry_validators:
            fvalidate = self._registry_validators[fvalidate]
        elif isinstance(fvalidate, str):
            raise ValueError("`fvalidate`, if str, must be in "
                             f"{self._registry_validators.keys()}")
        else:
            raise TypeError("`fvalidate` must be a function or "
                            f"{self._registry_validators.keys()}")
        self._fvalidate = fvalidate

    def __set_name__(self, cosmo_cls, name):
        # attribute name on container cosmology class
        self._attr_name = name
        self._attr_name_private = "_" + name

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
    # descriptor and property-like methods

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
        value = self.validate(cosmology, value)
        setattr(cosmology, self._attr_name_private, value)

    # -------------------------------------------
    # validate value

    @property
    def fvalidate(self):
        """Function to validate a potential value of this Parameter.."""
        return self._fvalidate

    def validator(self, fvalidate):
        """Make new Parameter with custom ``fvalidate``.

        Note: ``Parameter.fvalidator`` must be the top-most descriptor decorator.

        Parameters
        ----------
        fvalidate : callable[[type, type, Any], Any]

        Returns
        -------
        `~astropy.cosmology.Parameter`
            Copy of this Parameter but with custom ``fvalidate``.
        """
        return self.clone(fvalidate=fvalidate)

    def validate(self, cosmology, value):
        """Run the validator on this Parameter.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` instance
        value : Any
            The object to validate.

        Returns
        -------
        Any
            The output of calling ``fvalidate(cosmology, self, value)``
            (yes, that parameter order).
        """
        return self.fvalidate(cosmology, self, value)

    @classmethod
    def register_validator(cls, key, fvalidate=None):
        """Decorator to register a new kind of validator function.

        Parameters
        ----------
        key : str
        fvalidate : callable[[object, object, Any], Any] or None, optional
            Value validation function.

        Returns
        -------
        ``validator`` or callable[``validator``]
            if validator is None returns a function that takes and registers a
            validator. This allows ``register_validator`` to be used as a
            decorator.
        """
        if key in cls._registry_validators:
            raise KeyError(f"validator {key!r} already registered with Parameter.")

        # fvalidate directly passed
        if fvalidate is not None:
            cls._registry_validators[key] = fvalidate
            return fvalidate

        # for use as a decorator
        def register(fvalidate):
            """Register validator function.

            Parameters
            ----------
            fvalidate : callable[[object, object, Any], Any]
                Validation function.

            Returns
            -------
            ``validator``
            """
            cls._registry_validators[key] = fvalidate
            return fvalidate

        return register

    # -------------------------------------------

    def _get_init_arguments(self, processed=False):
        """Initialization arguments.

        Parameters
        ----------
        processed : bool
            Whether to more closely reproduce the input arguments (`False`,
            default) or the processed arguments (`True`). The former is better
            for string representations and round-tripping with ``eval(repr())``.

        Returns
        -------
        dict[str, Any]
        """
        # The keys are added in this order because `repr` prints them in order.
        kw = {"derived": self.derived,
              "unit": self.unit,
              "equivalencies": self.equivalencies,
              # Validator is always turned into a function, but for ``repr`` it's nice
              # to know if it was originally a string.
              "fvalidate": self.fvalidate if processed else self._fvalidate_in,
              "fmt": self.format_spec,
              "doc": self.__doc__}
        return kw

    def clone(self, **kw):
        """Clone this `Parameter`, changing any constructor argument.

        Parameters
        ----------
        **kw
            Passed to constructor. The current values, eg. ``fvalidate`` are
            used as the default values, so an empty ``**kw`` is an exact copy.

        Examples
        --------
        >>> p = Parameter()
        >>> p
        Parameter(derived=False, unit=None, equivalencies=[],
                  fvalidate='default', fmt='', doc=None)

        >>> p.clone(unit="km")
        Parameter(derived=False, unit=Unit("km"), equivalencies=[],
                  fvalidate='default', fmt='', doc=None)
        """
        # Start with defaults, update from kw.
        kwargs = {**self._get_init_arguments(), **kw}
        # All initialization failures, like incorrect input are handled by init
        cloned = type(self)(**kwargs)
        # Transfer over the __set_name__ stuff. If `clone` is used to make a
        # new descriptor, __set_name__ will be called again, overwriting this.
        cloned._attr_name = self._attr_name
        cloned._attr_name_private = self._attr_name_private

        return cloned

    def __eq__(self, other):
        """Check Parameter equality. Only equal to other Parameter objects.

        Returns
        -------
        NotImplemented or True
            `True` if equal, `NotImplemented` otherwise. This allows `other` to
            be check for equality with ``other.__eq__``.

        Examples
        --------
        >>> p1, p2 = Parameter(unit="km"), Parameter(unit="km")
        >>> p1 == p2
        True

        >>> p3 = Parameter(unit="km / s")
        >>> p3 == p1
        False

        >>> p1 != 2
        True
        """
        if not isinstance(other, Parameter):
            return NotImplemented
        # Check equality on all `_init_arguments` & `name`.
        # Need to compare the processed arguments because the inputs are many-
        # to-one, e.g. `fvalidate` can be a string or the equivalent function.
        return ((self._get_init_arguments(True) == other._get_init_arguments(True))
                and (self.name == other.name))

    def __repr__(self):
        """String representation.

        ``eval(repr())`` should work, depending if contents like ``fvalidate``
        can be similarly round-tripped.
        """
        return "Parameter({})".format(", ".join(f"{k}={v!r}" for k, v in
                                                self._get_init_arguments().items()))


# ===================================================================
# Built-in validators


@Parameter.register_validator("default")
def _validate_with_unit(cosmology, param, value):
    """
    Default Parameter value validator.
    Adds/converts units if Parameter has a unit.
    """
    if (unit := param.unit) is not None:
        # structured of unknown length
        if isinstance(unit, tuple) and (len(unit) == 2) and (unit[1] is Ellipsis):
            unit = unit[0]  # get unit from (unit, ...)

            # check if any sub-unit is wrong / missing.
            if (isinstance(value, u.Quantity)  # TODO! simplify this logic?
                and (value.unit == unit if not isinstance(value.unit, u.StructuredUnit)
                     else all(vu == unit for vu in value.unit.values()))):
                pass  # all sub-units are correct.
            else:  # need to convert
                with u.add_enabled_equivalencies(param.equivalencies):
                    value = u.Quantity(value, unit=unit, copy=False)

        # normal unit / pre-built structured
        else:
            with u.add_enabled_equivalencies(param.equivalencies):
                value = u.Quantity(value, unit=unit, copy=False)

    return value


@Parameter.register_validator("float")
def _validate_to_float(cosmology, param, value):
    """Parameter value validator with units, and converted to float."""
    value = _validate_with_unit(cosmology, param, value)
    return float(value)


@Parameter.register_validator("scalar")
def _validate_to_scalar(cosmology, param, value):
    """Parameter value validator, where value must be a scalar."""
    value = _validate_with_unit(cosmology, param, value)
    if not value.isscalar:
        raise ValueError(f"{param.name} is a non-scalar quantity")
    return value


@Parameter.register_validator("non-negative")
def _validate_non_negative(cosmology, param, value):
    """Parameter value validator where value is a positive float."""
    value = _validate_to_float(cosmology, param, value)
    if value < 0.0:
        raise ValueError(f"{param.name} cannot be negative.")
    return value
