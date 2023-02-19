# Licensed under a 3-clause BSD style license - see LICENSE.rst
from . import utils


class Base:
    """
    The abstract base class of all unit formats.
    """

    registry = {}
    _space = " "

    def __new__(cls, *args, **kwargs):
        # This __new__ is to make it clear that there is no reason to
        # instantiate a Formatter--if you try to you'll just get back the
        # class
        return cls

    def __init_subclass__(cls, **kwargs):
        # Keep a registry of all formats.  Key by the class name unless a name
        # is explicitly set (i.e., one *not* inherited from a superclass).
        if "name" not in cls.__dict__:
            cls.name = cls.__name__.lower()

        Base.registry[cls.name] = cls
        super().__init_subclass__(**kwargs)

    @classmethod
    def _get_unit_name(cls, unit):
        return unit.get_format_name(cls.name)

    @classmethod
    def parse(cls, s):
        """
        Convert a string to a unit object.
        """
        raise NotImplementedError(f"Can not parse with {cls.__name__} format")

    @classmethod
    def _format_superscript(cls, number):
        return f"^{number}"

    @classmethod
    def _format_unit_power(cls, unit, power=1):
        """Format the unit for this format class raised to the given power.

        This is overridden in Latex where the name of the unit can depend on the power
        (e.g., for degrees).
        """
        name = cls._get_unit_name(unit)
        if power != 1:
            name += cls._format_superscript(utils.format_power(power))
        return name

    @classmethod
    def _format_unit_list(cls, units):
        return cls._space.join(
            cls._format_unit_power(base_, power) for base_, power in units
        )

    @classmethod
    def to_string(cls, u):
        """
        Convert a unit object to a string.
        """
        raise NotImplementedError(f"Can not output in {cls.__name__} format")
