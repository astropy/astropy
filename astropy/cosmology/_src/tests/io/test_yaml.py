# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.units as u
from astropy.cosmology import Cosmology, FlatLambdaCDM, Planck18
from astropy.cosmology import units as cu
from astropy.cosmology._src.io.builtin.yaml import (
    from_yaml,
    to_yaml,
    yaml_constructor,
    yaml_representer,
)
from astropy.io.misc.yaml import AstropyDumper, dump, load

from .base import ToFromDirectTestBase, ToFromTestMixinBase

##############################################################################
# Test Serializer


def test_yaml_representer():
    """Test :func:`~astropy.cosmology._src.io.builtin.yaml.yaml_representer`."""
    # test function `representer`
    representer = yaml_representer("!astropy.cosmology.LambdaCDM")
    assert callable(representer)

    # test the normal method of dumping to YAML
    yml = dump(Planck18)
    assert isinstance(yml, str)
    assert yml.startswith("!astropy.cosmology.FlatLambdaCDM")


def test_yaml_constructor():
    """Test :func:`~astropy.cosmology._src.io.builtin.yaml.yaml_constructor`."""
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


class ToFromYAMLTestMixin(ToFromTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="yaml"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    @pytest.fixture
    def xfail_if_not_registered_with_yaml(self, cosmo_cls):
        """
        YAML I/O only works on registered classes. So the thing to check is
        if this class is registered. If not, :func:`pytest.xfail` this test.
        Some of the tests define custom cosmologies. They are not registered.
        """
        if cosmo_cls not in AstropyDumper.yaml_representers:
            pytest.xfail(
                f"Cosmologies of type {cosmo_cls} are not registered with YAML."
            )

    # ===============================================================

    def test_to_yaml(self, cosmo_cls, to_format, xfail_if_not_registered_with_yaml):
        """Test cosmology -> YAML."""
        yml = to_format("yaml")

        assert isinstance(yml, str)  # test type
        assert yml.startswith("!" + ".".join(cosmo_cls.__module__.split(".")[:2]))
        # e.g. "astropy.cosmology" for built-in cosmologies, or "__main__" for the test
        # SubCosmology class defined in ``astropy.cosmology._src.tests.test_core``.

    def test_from_yaml_default(
        self, cosmo, to_format, from_format, xfail_if_not_registered_with_yaml
    ):
        """Test cosmology -> YAML -> cosmology."""
        yml = to_format("yaml")

        got = from_format(yml, format="yaml")  # (cannot autoidentify)

        assert got.name == cosmo.name
        assert got.meta == cosmo.meta

        # it won't error if everything matches up
        got = from_format(yml, format="yaml")
        assert got == cosmo
        assert got.meta == cosmo.meta

        # auto-identify test moved because it doesn't work.
        # see test_from_yaml_autoidentify

    def test_from_yaml_autoidentify(
        self, cosmo, to_format, from_format, xfail_if_not_registered_with_yaml
    ):
        """As a non-path string, it does NOT auto-identifies 'format'.

        TODO! this says there should be different types of I/O registries.
              not just hacking object conversion on top of file I/O.
        """
        assert self.can_autodentify("yaml") is False

        # Showing the specific error. The str is interpreted as a file location
        # but is too long a file name.
        yml = to_format("yaml")
        with pytest.raises((FileNotFoundError, OSError)):  # OSError in Windows
            from_format(yml)

    # # TODO! this is a challenging test to write. It's also unlikely to happen.
    # def test_fromformat_subclass_partial_info_yaml(self, cosmo):
    #     """
    #     Test writing from an instance and reading from that class.
    #     This works with missing information.
    #     """

    # -----------------------------------------------------

    @pytest.mark.parametrize("format", [True, False, None])
    def test_is_equivalent_to_yaml(
        self, cosmo, to_format, format, xfail_if_not_registered_with_yaml
    ):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        This test checks that Cosmology equivalency can be extended to any
        Python object that can be converted to a Cosmology -- in this case
        a YAML string. YAML can't be identified without "format" specified.
        """
        obj = to_format("yaml")
        assert not isinstance(obj, Cosmology)

        is_equiv = cosmo.is_equivalent(obj, format=format)
        assert is_equiv is False

    def test_is_equivalent_to_yaml_specify_format(
        self, cosmo, to_format, xfail_if_not_registered_with_yaml
    ):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        Same as ``test_is_equivalent_to_yaml`` but with ``format="yaml"``.
        """
        assert cosmo.is_equivalent(to_format("yaml"), format="yaml") is True


class TestToFromYAML(ToFromDirectTestBase, ToFromYAMLTestMixin):
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
        This overrides from super because `ToFromDirectTestBase` adds a custom
        Cosmology ``CosmologyWithKwargs`` that is not registered with YAML.
        """
        return  # run tests

    def test_from_yaml_autoidentify(self, cosmo, to_format, from_format):
        """
        If directly calling the function there's no auto-identification.
        So this overrides the test from `ToFromYAMLTestMixin`
        """
