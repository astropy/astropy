# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology, w0wzCDM
from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.io.tests import (
    test_cosmology,
    test_ecsv,
    test_html,
    test_json,
    test_mapping,
    test_model,
    test_row,
    test_table,
    test_yaml,
)
from astropy.table import QTable, Row
from astropy.utils.compat.optional_deps import HAS_BS4

###############################################################################
# SETUP

cosmo_instances = cosmology.realizations.available

# Collect the registered read/write formats.
#   (format, supports_metadata, has_all_required_dependencies)
readwrite_formats = {
    ("ascii.ecsv", True, True),
    ("ascii.html", False, HAS_BS4),
    ("json", True, True),
}


# Collect all the registered to/from formats. Unfortunately this is NOT
# automatic since the output format class is not stored on the registry.
#                 (format, data type)
tofrom_formats = [
    ("mapping", dict),
    ("yaml", str),
    ("astropy.cosmology", Cosmology),
    ("astropy.row", Row),
    ("astropy.table", QTable),
]


###############################################################################


class ReadWriteTestMixin(
    test_ecsv.ReadWriteECSVTestMixin,
    test_html.ReadWriteHTMLTestMixin,
    test_json.ReadWriteJSONTestMixin,
):
    """
    Tests for a CosmologyRead/Write on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestReadWriteCosmology`` or ``TestCosmology`` for examples.
    """

    @pytest.mark.parametrize("format, metaio, has_deps", readwrite_formats)
    def test_readwrite_complete_info(self, cosmo, tmp_path, format, metaio, has_deps):
        """
        Test writing from an instance and reading from the base class.
        This requires full information.
        The round-tripped metadata can be in a different order, so the
        OrderedDict must be converted to a dict before testing equality.
        """
        if not has_deps:
            pytest.skip("missing a dependency")

        fname = str(tmp_path / f"{cosmo.name}.{format}")
        cosmo.write(fname, format=format)

        # Also test kwarg "overwrite"
        assert os.path.exists(fname)  # file exists
        with pytest.raises(IOError):
            cosmo.write(fname, format=format, overwrite=False)

        assert os.path.exists(fname)  # overwrite file existing file
        cosmo.write(fname, format=format, overwrite=True)

        # Read back
        got = Cosmology.read(fname, format=format)

        assert got == cosmo
        assert (not metaio) ^ (dict(got.meta) == dict(cosmo.meta))

    @pytest.mark.parametrize("format, metaio, has_deps", readwrite_formats)
    def test_readwrite_from_subclass_complete_info(
        self, cosmo_cls, cosmo, tmp_path, format, metaio, has_deps
    ):
        """
        Test writing from an instance and reading from that class, when there's
        full information saved.
        """
        if not has_deps:
            pytest.skip("missing a dependency")

        fname = str(tmp_path / f"{cosmo.name}.{format}")
        cosmo.write(fname, format=format)

        # read with the same class that wrote.
        got = cosmo_cls.read(fname, format=format)
        assert got == cosmo
        assert (not metaio) ^ (dict(got.meta) == dict(cosmo.meta))

        # this should be equivalent to
        got = Cosmology.read(fname, format=format, cosmology=cosmo_cls)
        assert got == cosmo
        assert (not metaio) ^ (dict(got.meta) == dict(cosmo.meta))

        # and also
        got = Cosmology.read(fname, format=format, cosmology=cosmo_cls.__qualname__)
        assert got == cosmo
        assert (not metaio) ^ (dict(got.meta) == dict(cosmo.meta))


class TestCosmologyReadWrite(ReadWriteTestMixin):
    """Test the classes CosmologyRead/Write."""

    @pytest.fixture(scope="class", params=cosmo_instances)
    def cosmo(self, request):
        return getattr(cosmology.realizations, request.param)

    @pytest.fixture(scope="class")
    def cosmo_cls(self, cosmo):
        return cosmo.__class__

    # ==============================================================

    @pytest.mark.parametrize("format, _, has_deps", readwrite_formats)
    def test_write_methods_have_explicit_kwarg_overwrite(self, format, _, has_deps):
        if not has_deps:
            pytest.skip("missing a dependency")

        writer = readwrite_registry.get_writer(format, Cosmology)
        # test in signature
        sig = inspect.signature(writer)
        assert "overwrite" in sig.parameters

        # also in docstring
        assert "overwrite : bool" in writer.__doc__

    @pytest.mark.parametrize("format, _, has_deps", readwrite_formats)
    def test_readwrite_reader_class_mismatch(
        self, cosmo, tmp_path, format, _, has_deps
    ):
        """Test when the reader class doesn't match the file."""
        if not has_deps:
            pytest.skip("missing a dependency")

        fname = tmp_path / f"{cosmo.name}.{format}"
        cosmo.write(fname, format=format)

        # class mismatch
        # when reading directly
        with pytest.raises(TypeError, match="missing 1 required"):
            w0wzCDM.read(fname, format=format)

        with pytest.raises(TypeError, match="missing 1 required"):
            Cosmology.read(fname, format=format, cosmology=w0wzCDM)

        # when specifying the class
        with pytest.raises(ValueError, match="`cosmology` must be either"):
            w0wzCDM.read(fname, format=format, cosmology="FlatLambdaCDM")


###############################################################################
# To/From_Format Tests


class ToFromFormatTestMixin(
    test_cosmology.ToFromCosmologyTestMixin,
    test_mapping.ToFromMappingTestMixin,
    test_model.ToFromModelTestMixin,
    test_row.ToFromRowTestMixin,
    test_table.ToFromTableTestMixin,
    test_yaml.ToFromYAMLTestMixin,
):
    """
    Tests for a Cosmology[To/From]Format on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.mark.parametrize("format, totype", tofrom_formats)
    def test_tofromformat_complete_info(
        self, cosmo, format, totype, xfail_if_not_registered_with_yaml
    ):
        """Read tests happen later."""
        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, totype)

        # test from_format
        got = Cosmology.from_format(obj, format=format)

        # Test autodetect, if enabled
        if self.can_autodentify(format):
            got2 = Cosmology.from_format(obj)
            assert got2 == got  # internal consistency

        assert got == cosmo  # external consistency
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format, totype", tofrom_formats)
    def test_fromformat_subclass_complete_info(
        self, cosmo_cls, cosmo, format, totype, xfail_if_not_registered_with_yaml
    ):
        """
        Test transforming an instance and parsing from that class, when there's
        full information available.
        Partial information tests are handled in the Mixin super classes.
        """
        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, totype)

        # read with the same class that wrote.
        got = cosmo_cls.from_format(obj, format=format)

        if self.can_autodentify(format):
            got2 = Cosmology.from_format(obj)  # and autodetect
            assert got2 == got  # internal consistency

        assert got == cosmo  # external consistency
        assert got.meta == cosmo.meta

        # this should be equivalent to
        got = Cosmology.from_format(obj, format=format, cosmology=cosmo_cls)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and also
        got = Cosmology.from_format(
            obj, format=format, cosmology=cosmo_cls.__qualname__
        )
        assert got == cosmo
        assert got.meta == cosmo.meta


class TestCosmologyToFromFormat(ToFromFormatTestMixin):
    """Test Cosmology[To/From]Format classes."""

    @pytest.fixture(scope="class", params=cosmo_instances)
    def cosmo(self, request):
        return getattr(cosmology.realizations, request.param)

    @pytest.fixture(scope="class")
    def cosmo_cls(self, cosmo):
        return cosmo.__class__

    # ==============================================================

    @pytest.mark.parametrize("format_type", tofrom_formats)
    def test_fromformat_class_mismatch(self, cosmo, format_type):
        format, totype = format_type

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, totype)

        # class mismatch
        with pytest.raises(TypeError):
            w0wzCDM.from_format(obj, format=format)

        with pytest.raises(TypeError):
            Cosmology.from_format(obj, format=format, cosmology=w0wzCDM)

        # when specifying the class
        with pytest.raises(ValueError, match="`cosmology` must be either"):
            w0wzCDM.from_format(obj, format=format, cosmology="FlatLambdaCDM")
