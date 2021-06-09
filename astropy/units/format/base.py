# Licensed under a 3-clause BSD style license - see LICENSE.rst

class Base:
    """
    The abstract base class of all unit formats.
    """
    registry = {}

    def __new__(cls, *args, **kwargs):
        # This __new__ is to make it clear that there is no reason to
        # instantiate a Formatter--if you try to you'll just get back the
        # class
        return cls

    def __init_subclass__(cls, **kwargs):
        # Keep a registry of all formats.  Key by the class name unless a name
        # is explicitly set (i.e., one *not* inherited from a superclass).
        if 'name' not in cls.__dict__:
            cls.name = cls.__name__.lower()

        Base.registry[cls.name] = cls
        super().__init_subclass__(**kwargs)

    @classmethod
    def parse(cls, s):
        """
        Convert a string to a unit object.
        """

        raise NotImplementedError(
            f"Can not parse with {cls.__name__} format")

    @classmethod
    def to_string(cls, u):
        """
        Convert a unit object to a string.
        """

        raise NotImplementedError(
            f"Can not output in {cls.__name__} format")
