# Licensed under a 3-clause BSD style license - see LICENSE.rst

# THIRD PARTY
import pytest

# LOCAL
import astropy.units as u
from astropy.cosmology import Cosmology, Parameter, realizations
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.parameters import available

cosmo_instances = [getattr(realizations, name) for name in available]


##############################################################################


class IOTestMixinBase:
    """
    Tests for a Cosmology[To/From]Format with some ``format``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmology`` for an example.
    """

    @pytest.fixture
    def from_format(self):
        """Convert to Cosmology using ``Cosmology.from_format()``."""
        return Cosmology.from_format

    @pytest.fixture
    def to_format(self, cosmo):
        """Convert Cosmology instance using ``.to_format()``."""
        return cosmo.to_format

    @pytest.fixture
    def read(self):
        """Read Cosmology instance using ``Cosmology.read()``."""
        return Cosmology.read

    @pytest.fixture
    def write(self, cosmo):
        """Write Cosmology using ``.write()``."""
        return cosmo.write


class IOFormatTestBase(IOTestMixinBase):
    """
    Directly test ``to/from_<format>``.
    These are not public API and are discouraged from public use, in favor of
    ``Cosmology.to/from_format(..., format="<format>")``.

    Subclasses should have an attribute ``functions`` which is a dictionary
    containing two items: ``"to"=<function for to_format>`` and
    ``"from"=<function for from_format>``.
    """

    @pytest.fixture(scope="class", autouse=True)
    def setup(self):
        """Setup and teardown for tests."""

        class CosmologyWithKwargs(Cosmology):
            Tcmb0 = Parameter(unit=u.K)

            def __init__(self, Tcmb0=0, name="cosmology with kwargs", meta=None, **kwargs):
                super().__init__(name=name, meta=meta)
                self._Tcmb0 = Tcmb0 << u.K

        yield  # run tests

        # pop CosmologyWithKwargs from registered classes
        # but don't error b/c it can fail in parallel
        _COSMOLOGY_CLASSES.pop(CosmologyWithKwargs.__qualname__, None)

    @pytest.fixture(params=cosmo_instances)
    def cosmo(self, request):
        """Cosmology instance."""
        if isinstance(request.param, str):  # CosmologyWithKwargs
            return _COSMOLOGY_CLASSES[request.param](Tcmb0=3)
        return request.param

    @pytest.fixture
    def cosmo_cls(self, cosmo):
        """Cosmology classes."""
        return cosmo.__class__

    # -------------------------------------------

    @pytest.fixture
    def from_format(self):
        """Convert to Cosmology using function ``from``."""
        def use_from_format(*args, **kwargs):
            kwargs.pop("format", None)  # specific to Cosmology.from_format
            return self.functions["from"](*args, **kwargs)

        return use_from_format

    @pytest.fixture
    def to_format(self, cosmo):
        """Convert Cosmology to format using function ``to``."""
        return lambda *args, **kwargs: self.functions["to"](cosmo, *args, **kwargs)

    # -------------------------------------------

    @pytest.fixture
    def read(self):
        """Read Cosmology from file using function ``read``."""
        def use_read(*args, **kwargs):
            kwargs.pop("format", None)  # specific to Cosmology.from_format
            return self.functions["read"](*args, **kwargs)

        return use_read

    @pytest.fixture
    def write(self, cosmo):
        """Write Cosmology to file using function ``write``."""
        return lambda *args, **kwargs: self.functions["write"](cosmo, *args, **kwargs)
