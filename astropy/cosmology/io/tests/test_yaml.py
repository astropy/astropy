# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import copy
import inspect
from collections import OrderedDict

# THIRD PARTY
import pytest

# LOCAL
import astropy.cosmology.units as cu
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
    # test function `representer`
    representer = yaml_representer("!astropy.cosmology.flrw.LambdaCDM")
    assert callable(representer)

    # test the normal method of dumping to YAML
    yml = dump(Planck18)
    assert isinstance(yml, str)
    assert yml.startswith("!astropy.cosmology.flrw.FlatLambdaCDM")


def test_yaml_constructor():
    """Test :func:`~astropy.cosmology.io.yaml.yaml_constructor`."""
    # test function `constructor`
    constructor = yaml_constructor(FlatLambdaCDM)
    assert callable(constructor)

    # it's too hard to manually construct a node, so we only test dump/load
    # this is also a good round-trip test
    yml = dump(Planck18)
    with u.add_enabled_units(cu):  # needed for redshift units
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

    @pytest.fixture
    def registered_with_yaml(self, cosmo_cls):
        """
        YAML I/O only works on registered classes. So the thing to check is
        if this class is registered. If not, skip this test.
        Some of the tests define custom cosmologies. They are not registered.
        """
        return True if cosmo_cls in AstropyDumper.yaml_representers else False

    def test_tofrom_yaml_instance(self, cosmo, to_format, from_format,
                                  registered_with_yaml):
        """Test cosmology -> YAML -> cosmology."""
        if not registered_with_yaml:
            return

        # ------------
        # To YAML

        yml = to_format('yaml')
        assert isinstance(yml, str)  # test type
        assert yml.startswith("!astropy.cosmology.")

        # ------------
        # From YAML

        got = from_format(yml, format="yaml")

        assert got.name == cosmo.name
        assert got.meta == cosmo.meta

        # it won't error if everything matches up
        got = from_format(yml, format="yaml")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # auto-identify test moved because it doesn't work.

    def test_tofrom_yaml_autoidentify(self, cosmo, to_format, from_format,
                                      registered_with_yaml):
        """As a non-path string, it does NOT auto-identifies 'format'.

        TODO! this says there should be different types of I/O registries.
              not just hacking object conversion on top of file I/O.
        """
        if not registered_with_yaml:
            return

        yml = to_format('yaml')
        with pytest.raises((FileNotFoundError, OSError)):  # OSError in Windows
            from_format(yml)

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
        """Set up fixtures to use ``to/from_yaml``, not the I/O abstractions."""
        self.functions = {"to": to_yaml, "from": from_yaml}

    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """
        Setup and teardown for tests.
        This overrides from super because `IOFormatTestBase` adds a custom
        Cosmology ``CosmologyWithKwargs`` that is not registered with YAML.
        """
        yield  # run tests

    def test_tofrom_yaml_autoidentify(self, cosmo, to_format, from_format):
        """
        If directly calling the function there's no auto-identification.
        So this overrides the test from `ToFromYAMLTestMixin`
        """
        pass
