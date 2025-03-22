# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.FLRW`."""

from typing import final

import pytest

from astropy.cosmology import FLRW
from astropy.cosmology._src.core import _COSMOLOGY_CLASSES, dataclass_decorator
from astropy.cosmology._src.tests.helper import get_redshift_methods
from astropy.cosmology._src.tests.test_core import invalid_zs
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .test_base import FLRWTest


@dataclass_decorator
class SubFLRW(FLRW):
    def w(self, z):
        return super().w(z)


@final
class TestFLRW(FLRWTest):
    """Test :class:`astropy.cosmology.FLRW`."""

    abstract_w = True

    def setup_class(self):
        """
        Setup for testing.
        FLRW is abstract, so tests are done on a subclass.
        """
        super().setup_class(self)

        # make sure SubCosmology is known
        _COSMOLOGY_CLASSES["SubFLRW"] = SubFLRW

        self.cls = SubFLRW

    def teardown_class(self):
        super().teardown_class(self)
        _COSMOLOGY_CLASSES.pop("SubFLRW", None)

    # ===============================================================
    # Method & Attribute Tests

    # ---------------------------------------------------------------
    # Methods

    def test_w(self, cosmo):
        """Test abstract :meth:`astropy.cosmology.FLRW.w`."""
        with pytest.raises(NotImplementedError, match="not implemented"):
            cosmo.w(1)

    def test_Otot(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.Otot`."""
        exception = NotImplementedError if HAS_SCIPY else ModuleNotFoundError
        with pytest.raises(exception):
            assert cosmo.Otot(1)

    def test_efunc_vs_invefunc(self, cosmo):
        """
        Test that efunc and inv_efunc give inverse values.
        Here they just fail b/c no ``w(z)`` or no scipy.
        """
        exception = NotImplementedError if HAS_SCIPY else ModuleNotFoundError

        with pytest.raises(exception):
            cosmo.efunc(0.5)

        with pytest.raises(exception):
            cosmo.inv_efunc(0.5)

    @pytest.mark.skip(reason="w(z) is abstract")
    def test_luminosity_distance_pandas(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.luminosity_distance`."""

    _FLRW_redshift_methods = get_redshift_methods(
        FLRW, include_private=True, include_z2=False
    ) - {"w"}

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy is not installed")
    @pytest.mark.parametrize("z, exc", invalid_zs)
    @pytest.mark.parametrize("method", sorted(_FLRW_redshift_methods))
    def test_redshift_method_bad_input(self, cosmo, method, z, exc):
        """Test all the redshift methods for bad input."""
        with pytest.raises(exc):
            getattr(cosmo, method)(z)

    # ===============================================================
    # Usage Tests

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    @pytest.mark.parametrize("method", ("Om", "Ode", "w", "de_density_scale"))
    def test_distance_broadcast(self, cosmo, method):
        with pytest.raises(NotImplementedError):
            super().test_distance_broadcast(cosmo, method)

    @pytest.mark.skip(reason="w(z) is abstract")
    def test_comoving_distance_1arg_equal_to_2arg(self, cosmo):
        """Test :meth:`astropy.cosmology.FLRW.luminosity_distance`."""

    @pytest.mark.skipif(not HAS_SCIPY, reason="scipy required for this test.")
    @pytest.mark.parametrize(
        ("args", "kwargs", "expected"),
        [((70, 0.27, 0.73), {"Tcmb0": 3.0, "Ob0": 0.03}, None)],
    )
    def test_comoving_distance_example(self, cosmo_cls, args, kwargs, expected):
        with pytest.raises(NotImplementedError):
            super().test_comoving_distance_example(cosmo_cls, args, kwargs, expected)
