# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import copy
import inspect
from collections import OrderedDict

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import Cosmology, FlatLambdaCDM, Planck18, realizations
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter
from astropy.cosmology.io.yaml import from_yaml, to_yaml, yaml_constructor, yaml_representer
from astropy.cosmology.parameters import available
from astropy.io.misc.yaml import AstropyDumper, AstropyLoader, dump, load
from astropy.table import QTable, vstack

from .base import IOTestMixinBase, IOFormatTestBase

cosmo_instances = [getattr(realizations, name) for name in available]
# cosmo_instances.append("TestToFromYAML.setup.<locals>.CosmologyWithKwargs")


##############################################################################
# Test Serializer


def test_yaml_representer():
    """Test :func:`~astropy.cosmology.io.yaml.yaml_representer`."""
    # test `yaml_representer`
    assert callable(yaml_representer)

    sig = inspect.signature(yaml_representer)
    assert len(sig.parameters) == 1
    assert "tag" in sig.parameters
    assert sig.parameters["tag"].default is inspect._empty

    # test function `representer`
    representer = yaml_representer("!astropy.cosmology.flrw.LambdaCDM")
    assert callable(representer)

    sig = inspect.signature(representer)
    assert len(sig.parameters) == 2
    assert "dumper" in sig.parameters
    assert "obj" in sig.parameters
    assert sig.parameters["dumper"].default is inspect._empty
    assert sig.parameters["obj"].default is inspect._empty

    # test the normal method of dumping to YAML
    yml = dump(Planck18)
    assert isinstance(yml, str)
    assert yml.startswith("!astropy.cosmology.flrw.FlatLambdaCDM")


def test_yaml_constructor():
    """Test :func:`~astropy.cosmology.io.yaml.yaml_constructor`."""
    # test `yaml_representer`
    assert callable(yaml_constructor)

    sig = inspect.signature(yaml_constructor)
    assert len(sig.parameters) == 1
    assert "cls" in sig.parameters
    assert sig.parameters["cls"].default is inspect._empty

    # test function `constructor`
    constructor = yaml_constructor(FlatLambdaCDM)
    assert callable(constructor)

    sig = inspect.signature(constructor)
    assert len(sig.parameters) == 2
    assert "loader" in sig.parameters
    assert "node" in sig.parameters
    assert sig.parameters["loader"].default is inspect._empty
    assert sig.parameters["node"].default is inspect._empty

    # it's too hard to manually construct a node, so we only test dump/load
    # this is also a good round-trip test
    yml = dump(Planck18)
    cosmo = load(yml)
    assert isinstance(cosmo, FlatLambdaCDM)
    assert cosmo == Planck18
    assert cosmo.meta == Planck18.meta


##############################################################################
# Test Unified I/O


class ToFromYAMLTestMixin(IOTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="yaml"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    def test_to_from_yaml_instance(self, cosmo, to_format, from_format):
        """Test cosmology -> YAML -> cosmology."""
        # ------------
        # To YAML

        yml = to_format('yaml')
        assert isinstance(yml, str)  # test type
        assert yml.startswith("!astropy.cosmology.")

        # ------------
        # From YAML

        # tests are different if the last argument is a **kwarg
        if tuple(cosmo._init_signature.parameters.values())[-1].kind == 4:
            got = from_format(yml, format="yaml")

            assert got.name == cosmo.name
            assert got.meta == cosmo.meta

            return  # don't continue testing

        # it won't error if everything matches up
        got = from_format(yml, format="yaml")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # it auto-identifies 'format'
        got = from_format(yml)
        assert got == cosmo
        assert got.meta == cosmo.meta

    # # TODO! this is a challenging test to write. It's also unlikely to happen.
    # def test_fromformat_subclass_partial_info_yaml(self, cosmo):
    #     """
    #     Test writing from an instance and reading from that class.
    #     This works with missing information.
    #     """


class TestToFromYAML(IOFormatTestBase, ToFromYAMLTestMixin):
    """
    Directly test ``to/from_yaml``.
    These are not public API and are discouraged from use, in favor of
    ``Cosmology.to/from_format(..., format="yaml")``, but should be tested
    regardless b/c 3rd party packages might use these in their Cosmology I/O.
    Also, it's cheap to test.
    """

    def setup_class(self):
        self.functions = {"to": to_yaml, "from": from_yaml}

    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """Setup and teardown for tests."""
        yield  # run tests
