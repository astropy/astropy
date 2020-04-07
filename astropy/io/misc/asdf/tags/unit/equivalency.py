# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from astropy.units.equivalencies import Equivalency
from astropy.units import equivalencies
from astropy.units.quantity import Quantity

from astropy.io.misc.asdf.types import AstropyType


class EquivalencyType(AstropyType):
    name = "units/equivalency"
    types = [Equivalency]
    version = '1.0.0'

    @classmethod
    def to_tree(cls, equiv, ctx):
        node = {}
        if not isinstance(equiv, Equivalency):
            raise TypeError(f"'{equiv}' is not a valid Equivalency")

        eqs = []
        for e, kwargs in zip(equiv.name, equiv.kwargs):
            kwarg_names = list(kwargs.keys())
            kwarg_values = list(kwargs.values())
            eq = {'name': e, 'kwargs_names': kwarg_names, 'kwargs_values': kwarg_values}
            eqs.append(eq)
        return eqs

    @classmethod
    def from_tree(cls, node, ctx):
        eqs = []
        for eq in node:
            equiv = getattr(equivalencies, eq['name'])
            kwargs = dict(zip(eq['kwargs_names'], eq['kwargs_values']))
            eqs.append(equiv(**kwargs))
        return sum(eqs[1:], eqs[0])

    @classmethod
    def assert_equal(cls, a, b):
        assert a == b
