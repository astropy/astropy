# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import inspect
import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology, w0wzCDM
from astropy.cosmology.connect import CosmologyRead, readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.cosmology.io.tests import test_ecsv, test_mapping, test_model, test_table
from astropy.table import QTable

from .conftest import json_identify, read_json, write_json

###############################################################################
# SETUP

cosmo_instances = cosmology.parameters.available
readwrite_formats = ["json"]
tofrom_formats = [("mapping", dict), ("astropy.table", QTable)]
#                 (format, data type)

###############################################################################


class ReadWriteTestMixin:
    """
    Tests for a CosmologyRead/Write on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestReadWriteCosmology`` or ``TestCosmology`` for examples.
    """

    @pytest.fixture(scope="class", autouse=True)
    def setup_readwrite(self):
        """Setup & teardown for read/write tests."""
        # register
        readwrite_registry.register_reader("json", Cosmology, read_json, force=True)
        readwrite_registry.register_writer("json", Cosmology, write_json, force=True)
        readwrite_registry.register_identifier("json", Cosmology, json_identify, force=True)

        yield  # run all tests in class

        # unregister
        readwrite_registry.unregister_reader("json", Cosmology)
        readwrite_registry.unregister_writer("json", Cosmology)
        readwrite_registry.unregister_identifier("json", Cosmology)

    # ===============================================================
    # Method & Attribute Tests

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_readwrite_complete_info(self, cosmo, tmpdir, format):
        """
        Test writing from an instance and reading from the base class.
        This requires full information.
        """
        fname = tmpdir / f"{cosmo.name}.{format}"

        cosmo.write(str(fname), format=format)

        # Also test kwarg "overwrite"
        assert os.path.exists(str(fname))  # file exists
        with pytest.raises(IOError):
            cosmo.write(str(fname), format=format, overwrite=False)

        assert os.path.exists(str(fname))  # overwrite file existing file
        cosmo.write(str(fname), format=format, overwrite=True)

        # Read back
        got = Cosmology.read(fname, format=format)

        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_readwrite_from_subclass_complete_info(self, cosmo, tmpdir, format):
        """
        Test writing from an instance and reading from that class, when there's
        full information saved.
        """
        fname = tmpdir / f"{cosmo.name}.{format}"
        cosmo.write(str(fname), format=format)

        # read with the same class that wrote.
        got = cosmo.__class__.read(fname, format=format)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # this should be equivalent to
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__)
        assert got == cosmo
        assert got.meta == cosmo.meta

        # and also
        got = Cosmology.read(fname, format=format, cosmology=cosmo.__class__.__qualname__)
        assert got == cosmo
        assert got.meta == cosmo.meta

    @pytest.mark.skip("TODO: generalize over all save formats for this test.")
    def test_readwrite_from_subclass_partial_info(self, cosmo, tmpdir):
        """
        Test writing from an instance and reading from that class.
        This requires partial information.

        .. todo::

            - generalize over all save formats for this test.
            - remove the method defined in subclass
        """


class TestCosmologyReadWrite(ReadWriteTestMixin):
    """Test the classes CosmologyRead/Write."""

    @pytest.fixture(params=cosmo_instances)
    def cosmo(self, request):
        return getattr(cosmology.realizations, request.param)

    # ==============================================================

    @pytest.mark.parametrize("format", readwrite_formats)
    def test_write_methods_have_explicit_kwarg_overwrite(self, format):
        writer = readwrite_registry.get_writer(format, Cosmology)
        # test in signature
        sig = inspect.signature(writer)
        assert "overwrite" in sig.parameters

        # also in docstring
        assert "overwrite : bool" in writer.__doc__

    # @pytest.mark.parametrize("format", readwrite_formats)
    @pytest.mark.parametrize("instance", cosmo_instances)
    def test_readwrite_from_subclass_partial_info(self, instance, tmpdir):
        """
        Test writing from an instance and reading from that class.
        This requires partial information.

        .. todo::

            remove when fix method in super
        """
        cosmo = getattr(cosmology.realizations, instance)

        format = "json"
        fname = tmpdir / f"{cosmo.name}.{format}"

        cosmo.write(str(fname), format=format)

        # partial information
        with open(fname, "r") as file:
            L = file.readlines()[0]
        L = L[: L.index('"cosmology":')] + L[L.index(", ") + 2 :]  # remove cosmology
        i = L.index('"Tcmb0":')  # delete Tcmb0
        L = L[:i] + L[L.index(", ", L.index(", ", i) + 1) + 2 :]  # second occurence

        tempfname = tmpdir / f"{cosmo.name}_temp.{format}"
        with open(tempfname, "w") as file:
            file.writelines([L])

        # read with the same class that wrote fills in the missing info with
        # the default value
        got = cosmo.__class__.read(tempfname, format=format)
        got2 = Cosmology.read(tempfname, format=format, cosmology=cosmo.__class__)
        got3 = Cosmology.read(tempfname, format=format, cosmology=cosmo.__class__.__qualname__)

        assert (got == got2) and (got2 == got3)  # internal consistency

        # not equal, because Tcmb0 is changed
        assert got != cosmo
        assert got.Tcmb0 == cosmo.__class__._init_signature.parameters["Tcmb0"].default
        assert got.clone(name=cosmo.name, Tcmb0=cosmo.Tcmb0) == cosmo
        # but the metadata is the same
        assert got.meta == cosmo.meta

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


class ToFromFormatTestMixin(
    # convert
    test_mapping.ToFromMappingTestMixin, test_model.ToFromModelTestMixin,
    test_table.ToFromTableTestMixin,
    # read/write
    test_ecsv.ReadWriteECSVTestMixin,
):
    """
    Tests for a Cosmology[To/From]Format on a |Cosmology|.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.mark.parametrize("format_type", tofrom_formats)
    def test_tofromformat_complete_info(self, cosmo, format_type):
        """Read tests happen later."""
        format, objtype = format_type

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # test from_format
        got = Cosmology.from_format(obj, format=format)
        # and autodetect
        got2 = Cosmology.from_format(obj)

        assert got2 == got  # internal consistency
        assert got == cosmo  # external consistency
        assert got.meta == cosmo.meta

    @pytest.mark.parametrize("format_type", tofrom_formats)
    def test_fromformat_subclass_complete_info(self, cosmo, format_type):
        """
        Test transforming an instance and parsing from that class, when there's
        full information available.
        Partial information tests are handled in the Mixin super classes.
        """
        format, objtype = format_type

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # read with the same class that wrote.
        got = cosmo.__class__.from_format(obj, format=format)
        got2 = Cosmology.from_format(obj)  # and autodetect

        assert got2 == got  # internal consistency
        assert got == cosmo  # external consistency
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

    @pytest.fixture(params=cosmo_instances)
    def cosmo(self, request):
        return getattr(cosmology.realizations, request.param)

    @pytest.fixture
    def cosmo_cls(self, cosmo):
        return cosmo.__class__

    # ==============================================================

    @pytest.mark.parametrize("format_type", tofrom_formats)
    def test_fromformat_class_mismatch(self, cosmo, format_type):
        format, objtype = format_type

        # test to_format
        obj = cosmo.to_format(format)
        assert isinstance(obj, objtype)

        # class mismatch
        with pytest.raises(TypeError, match="missing 1 required"):
            w0wzCDM.from_format(obj, format=format)

        with pytest.raises(TypeError, match="missing 1 required"):
            Cosmology.from_format(obj, format=format, cosmology=w0wzCDM)

        # when specifying the class
        with pytest.raises(ValueError, match="`cosmology` must be either"):
            w0wzCDM.from_format(obj, format=format, cosmology="FlatLambdaCDM")
