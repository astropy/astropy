# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from astropy.units.equivalencies import Equivalency
from astropy.units import equivalencies
from astropy.units.quantity import Quantity
from asdf.yamlutil import custom_tree_to_tagged_tree

from astropy.io.misc.asdf.types import AstropyAsdfType


class EquivalencyType(AstropyAsdfType):
    name = "unit/equivalency"
    types = [Equivalency]
    requires = ['astropy']
    version = '1.1.0'

    @classmethod
    def to_tree(cls, equiv, ctx):
        node = {}
        if isinstance(equiv, Equivalency):
            eqs = []
            for i, e in enumerate(equiv.name):
                args = equiv.args[i]
                args = [custom_tree_to_tagged_tree(i, ctx) if isinstance(i, Quantity) \
                        else i for i in args]
                kwargs = equiv.kwargs[i]
                eq = {'name': e, 'args': args,  'kwargs_names': list(kwargs.keys()),
                      'kwargs_values': list(kwargs.values())}
                eqs.append(eq)
        else:
            raise TypeError("'{0}' is not a valid Equivalency".format(equiv))
        node['equivalencies'] = eqs
        return node
     


    @classmethod
    def from_tree(cls, node, ctx):
        import operator
        e = []
        for eq in node['equivalencies']:
            equiv = getattr(equivalencies, eq['name'])
            args = eq['args']
            kwargs = dict(zip(eq['kwargs_names'], eq['kwargs_values']))
            e.append(equiv(*args, **kwargs))
        return operator.add(*e)
