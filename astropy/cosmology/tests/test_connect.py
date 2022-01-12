# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import inspect
import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology, w0wzCDM
from astropy.cosmology.connect import CosmologyRead, readwrite_registry, convert_registry
from astropy.cosmology.core import Cosmology
from astropy.cosmology.io.tests import (test_ecsv, test_json, test_mapping,
                                        test_model, test_row, test_table, test_yaml)
from astropy.modeling import Model
from astropy.table import QTable, Row

###############################################################################
# SETUP

cosmo_instances = cosmology.parameters.available

# Collect all the registered read/write formats. JSON is added because it is
# registered separately, in ``setup_readwrite``.
readwrite_formats = {k for k, cls in readwrite_registry._readers.keys()}.union({"json"})

# Collect all the registered to/from formats. Unfortunately this is NOT
# automatic since the output format class is not stored on the registry.
#                 (format, data type)
# TODO! figure out how to aumate collection.
tofrom_formats = [("mapping", dict), ("yaml", str),
                  # ("astropy.model", Model),
                  ("astropy.row", Row), ("astropy.table", QTable)]
# tofrom_formats += [("mypackage", ...)]


def test_collected_all_expected_tofrom_formats():
    """Test that all the expected ``to/from_format`` will be tested.

    Unfortunately collecting all these registrants is NOT automatic since
    the output format class is not stored on the registry.

    .. todo:: store the output format class on the registry?
    """
    got = {k for (k, _) in tofrom_formats}

    expect = {k for (k, _) in convert_registry._readers.keys()}
    expect -= {"astropy.model"}  # TODO! requires kwarg 'method'
    expect -= {"mypackage"}  # TODO! only works on FLRW subclasses, shouldn't skip those

    assert got == expect


###############################################################################


class ReadWriteTestMixin(test_ecsv.ReadWriteECSVTestMixin, test_json.ReadWriteJSONTestMixin):
    """
    Tests for a CosmologyRead/Write on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestReadWriteCosmology`` or ``TestCosmology`` for examples.
    """

    # TODO! "myformat" only works on FLRW subclasses. Detect and only skip needed!
    @pytest.mark.parametrize("format", readwrite_formats - {"myformat"})
    def test_readwrite_complete_info(self, cosmo, tmpdir, format):
        """
        Test writing from an instance and reading from the base class.
        This requires full information.

        The round-tripped metadata can be in a different order, so the
        OrderedDict must be converted to a dict before testing equality.
        """
        fname = tmpdir / f"{cosmo.name}.{format}"

        # Write the cosmology instance to a file.
        # Format-specific tests should confirm that writing to a file works,
        # so if it doesn't here, we just skip the test.
        try:
            cosmo.write(str(fname), format=format)
        except:
            pytest.xfail(f"failing write for class {cosmo.__class__} with format={format}")

        # Test that the method accepts the kwarg "overwrite" and that if
        # the file exists it will error if "overwrite" is not True.
        assert os.path.exists(fname)  # file exists
        with pytest.raises(IOError):
            cosmo.write(fname, format=format, overwrite=False)

        assert os.path.exists(fname)  # overwrite existing file
        cosmo.write(fname, format=format, overwrite=True)

        # Now test that round-tripping works as expected
        got = Cosmology.read(fname, format=format)
        assert got == cosmo
        # Test metadata equality. Roundtrip does NOT always preserve order
        # so need to convert OrderedDict to dict.
        assert dict(got.meta) == dict(cosmo.meta)

    # TODO! "myformat" only works on FLRW subclasses. Detect and only skip needed!
    @pytest.mark.parametrize("format", readwrite_formats - {"myformat"})
    def test_readwrite_from_subclass_complete_info(self, cosmo, tmpdir, format):
        """
        Test writing from an instance and reading from that class, when there's
        full information saved.
        The roundtripped metadata can be in a different order, so need to
        convert from OrderedDict to dict before testing equality.
        """
        fname = tmpdir / f"{cosmo.name}.{format}"

        # Write the cosmology instance to a file.
        # Format-specific tests should confirm that writing to a file works,
        # so if it doesn't here, we just skip the test.
        try:
            cosmo.write(str(fname), format=format)
        except:
            pytest.xfail(f"failing write for class {cosmo.__class__} with format={format}")

        # read with the same class that wrote.
        got = cosmo.__class__.read(fname, format=format)
        assert got == cosmo
        assert dict(got.meta) == dict(cosmo.meta)

        # this should be equivalent to
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__)
        assert got == cosmo
        assert dict(got.meta) == dict(cosmo.meta)

        # and also
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__.__qualname__)
        assert got == cosmo
        assert dict(got.meta) == dict(cosmo.meta)


class TestCosmologyReadWrite(ReadWriteTestMixin):
    """Test the classes CosmologyRead/Write."""

    @pytest.fixture(scope="class", params=cosmo_instances)
    def cosmo(self, request):
        return getattr(cosmology.realizations, request.param)

    @pytest.fixture(scope="class")
    def cosmo_cls(self, cosmo):
        return cosmo.__class__

    # ==============================================================

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_write_methods_have_explicit_kwarg_overwrite(self, format):
        writer = readwrite_registry.get_writer(format, Cosmology)
        # test in signature
        sig = inspect.signature(writer)
        assert "overwrite" in sig.parameters

        # also in docstring
        assert "overwrite : bool" in writer.__doc__

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_readwrite_reader_class_mismatch(self, cosmo, tmpdir, format):
        """Test when the reader class doesn't match the file."""

        fname = tmpdir / f"{cosmo.name}.{format}"
        cosmo.write(str(fname), format=format)

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


class ToFromFormatTestMixin(test_mapping.ToFromMappingTestMixin, test_model.ToFromModelTestMixin,
                            test_row.ToFromRowTestMixin, test_table.ToFromTableTestMixin,
                            test_yaml.ToFromYAMLTestMixin):
    """
    Tests for a Cosmology[To/From]Format on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.mark.parametrize("format, totype", tofrom_formats)
    def test_tofromformat_complete_info(self, cosmo, format, totype):
        """Read tests happen later."""
        # Mark tests expected to fail
        self.mark_xfail_if_not_registered_with_yaml(cosmo.__class__, format)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, totype)

        # test from_format
        got = Cosmology.from_format(obj, format=format)
        # and autodetect, if enabled
        if self.can_autodentify(format):
            got2 = Cosmology.from_format(obj)
            assert got2 == got  # internal consistency

        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format, totype", tofrom_formats)
    def test_fromformat_subclass_complete_info(self, cosmo, format, totype):
        """
        Test transforming an instance and parsing from that class, when there's
        full information available.
        Partial information tests are handled in the Mixin super classes.
        """
        # Mark tests expected to fail
        self.mark_xfail_if_not_registered_with_yaml(cosmo.__class__, format)

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, totype)

        # read with the same class that wrote.
        got = cosmo.__class__.from_format(obj, format=format)
        # and autodetect, if enabled
        if self.can_autodentify(format):
            got2 = Cosmology.from_format(obj)
            assert got2 == got  # internal consistency

        assert got == cosmo  # consistency
        assert got.meta == cosmo.meta

        # this should be equivalent to
        got = Cosmology.from_format(obj, format=format, cosmology=cosmo.__class__)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and also
        got = Cosmology.from_format(obj, format=format, cosmology=cosmo.__class__.__qualname__)
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

    @pytest.mark.parametrize("format, totype", tofrom_formats)
    def test_fromformat_class_mismatch(self, cosmo, format, totype):
        # Mark tests expected to fail
        self.mark_xfail_if_not_registered_with_yaml(cosmo.__class__, format)

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
