# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.units as u
from astropy.cosmology import Cosmology, Parameter, realizations
from astropy.cosmology import units as cu
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.realizations import available

cosmo_instances = [getattr(realizations, name) for name in available]


##############################################################################


class IOTestBase:
    """Base class for Cosmology I/O tests.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """


class ToFromTestMixinBase(IOTestBase):
    """Tests for a Cosmology[To/From]Format with some ``format``.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.fixture(scope="class")
    def from_format(self):
        """Convert to Cosmology using ``Cosmology.from_format()``."""
        return Cosmology.from_format

    @pytest.fixture(scope="class")
    def to_format(self, cosmo):
        """Convert Cosmology instance using ``.to_format()``."""
        return cosmo.to_format

    def can_autodentify(self, format):
        """Check whether a format can auto-identify."""
        return format in Cosmology.from_format.registry._identifiers


class ReadWriteTestMixinBase(IOTestBase):
    """Tests for a Cosmology[Read/Write].

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.fixture(scope="class")
    def read(self):
        """Read Cosmology instance using ``Cosmology.read()``."""
        return Cosmology.read

    @pytest.fixture(scope="class")
    def write(self, cosmo):
        """Write Cosmology using ``.write()``."""
        return cosmo.write

    @pytest.fixture
    def add_cu(self):
        """Add :mod:`astropy.cosmology.units` to the enabled units."""
        # TODO! autoenable 'cu' if cosmology is imported?
        with u.add_enabled_units(cu):
            yield


##############################################################################


class IODirectTestBase(IOTestBase):
    """Directly test Cosmology I/O functions.

    These functions are not public API and are discouraged from public use, in
    favor of the I/O methods on |Cosmology|. They are tested b/c they are used
    internally and because some tests for the methods on |Cosmology| don't need
    to be run in the |Cosmology| class's large test matrix.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass.
    """

    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """Setup and teardown for tests."""

        class CosmologyWithKwargs(Cosmology):
            Tcmb0 = Parameter(default=0, unit=u.K)

            def __init__(
                self, Tcmb0=0, name="cosmology with kwargs", meta=None, **kwargs
            ):
                super().__init__(name=name, meta=meta)
                self._Tcmb0 = Tcmb0 << u.K

        yield  # run tests

        # pop CosmologyWithKwargs from registered classes
        # but don't error b/c it can fail in parallel
        _COSMOLOGY_CLASSES.pop(CosmologyWithKwargs.__qualname__, None)

    @pytest.fixture(scope="class", params=cosmo_instances)
    def cosmo(self, request):
        """Cosmology instance."""
        if isinstance(request.param, str):  # CosmologyWithKwargs
            return _COSMOLOGY_CLASSES[request.param](Tcmb0=3)
        return request.param

    @pytest.fixture(scope="class")
    def cosmo_cls(self, cosmo):
        """Cosmology classes."""
        return cosmo.__class__


class ToFromDirectTestBase(IODirectTestBase, ToFromTestMixinBase):
    """Directly test ``to/from_<format>``.

    These functions are not public API and are discouraged from public use, in
    favor of ``Cosmology.to/from_format(..., format="<format>")``. They are
    tested because they are used internally and because some tests for the
    methods on |Cosmology| don't need to be run in the |Cosmology| class's
    large test matrix.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass.

    Subclasses should have an attribute ``functions`` which is a dictionary
    containing two items: ``"to"=<function for to_format>`` and
    ``"from"=<function for from_format>``.
    """

    @pytest.fixture(scope="class")
    def from_format(self):
        """Convert to Cosmology using function ``from``."""

        def use_from_format(*args, **kwargs):
            kwargs.pop("format", None)  # specific to Cosmology.from_format
            return self.functions["from"](*args, **kwargs)

        return use_from_format

    @pytest.fixture(scope="class")
    def to_format(self, cosmo):
        """Convert Cosmology to format using function ``to``."""

        def use_to_format(*args, **kwargs):
            return self.functions["to"](cosmo, *args, **kwargs)

        return use_to_format


class ReadWriteDirectTestBase(IODirectTestBase, ToFromTestMixinBase):
    """Directly test ``read/write_<format>``.

    These functions are not public API and are discouraged from public use, in
    favor of ``Cosmology.read/write(..., format="<format>")``. They are tested
    because they are used internally and because some tests for the
    methods on |Cosmology| don't need to be run in the |Cosmology| class's
    large test matrix.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass.

    Subclasses should have an attribute ``functions`` which is a dictionary
    containing two items: ``"read"=<function for read>`` and
    ``"write"=<function for write>``.
    """

    @pytest.fixture(scope="class")
    def read(self):
        """Read Cosmology from file using function ``read``."""

        def use_read(*args, **kwargs):
            kwargs.pop("format", None)  # specific to Cosmology.from_format
            return self.functions["read"](*args, **kwargs)

        return use_read

    @pytest.fixture(scope="class")
    def write(self, cosmo):
        """Write Cosmology to file using function ``write``."""

        def use_write(*args, **kwargs):
            return self.functions["write"](cosmo, *args, **kwargs)

        return use_write
