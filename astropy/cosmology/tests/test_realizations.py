# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import pickle

# THIRD PARTY
import pytest

# LOCAL
import astropy.cosmology.units as cu
import astropy.units as u
from astropy import cosmology
from astropy.cosmology import parameters, realizations
from astropy.cosmology.realizations import Planck13, default_cosmology
from astropy.tests.helper import pickle_protocol


def test_realizations_in_toplevel_dir():
    """Test the realizations are in ``dir`` of :mod:`astropy.cosmology`."""
    d = dir(cosmology)

    assert set(d) == set(cosmology.__all__)
    for n in parameters.available:
        assert n in d


def test_realizations_in_realizations_dir():
    """Test the realizations are in ``dir`` of :mod:`astropy.cosmology.realizations`."""
    d = dir(realizations)

    assert set(d) == set(realizations.__all__)
    for n in parameters.available:
        assert n in d


class Test_default_cosmology(object):
    """Tests for :class:`~astropy.cosmology.realizations.default_cosmology`."""

    # -----------------------------------------------------
    # Get

    def test_get_fail(self):
        """Test bad inputs to :meth:`astropy.cosmology.default_cosmology.get`."""
        # a not-valid option, but still a str
        with pytest.raises(ValueError, match="Unknown cosmology"):
            cosmo = default_cosmology.get("fail!")

        # a not-valid type
        with pytest.raises(TypeError, match="'key' must be must be"):
            cosmo = default_cosmology.get(object())

    def test_get_current(self):
        """Test :meth:`astropy.cosmology.default_cosmology.get` current value."""
        cosmo = default_cosmology.get(None)
        assert cosmo is default_cosmology.get(default_cosmology._value)

    def test_get_none(self):
        """Test :meth:`astropy.cosmology.default_cosmology.get` to `None`."""
        cosmo = default_cosmology.get("no_default")
        assert cosmo is None

    @pytest.mark.parametrize("name", parameters.available)
    def test_get_valid(self, name):
        """Test :meth:`astropy.cosmology.default_cosmology.get` from str."""
        cosmo = default_cosmology.get(name)
        assert cosmo is getattr(realizations, name)

    def test_get_cosmology_from_string(self, recwarn):
        """Test method ``get_cosmology_from_string``."""
        cosmo = default_cosmology.get_cosmology_from_string("no_default")
        assert cosmo is None

        cosmo = default_cosmology.get_cosmology_from_string("Planck13")
        assert cosmo is Planck13

        with pytest.raises(ValueError):
            cosmo = default_cosmology.get_cosmology_from_string("fail!")

    # -----------------------------------------------------
    # Validate

    def test_validate_fail(self):
        """Test :meth:`astropy.cosmology.default_cosmology.validate`."""
        # bad input type
        with pytest.raises(TypeError, match="must be a string or Cosmology"):
            default_cosmology.validate(TypeError)

    def test_validate_default(self):
        """Test method ``validate`` for specific values."""
        value = default_cosmology.validate(None)
        assert value is realizations.Planck18

    @pytest.mark.parametrize("name", parameters.available)
    def test_validate_str(self, name):
        """Test method ``validate`` for string input."""
        value = default_cosmology.validate(name)
        assert value is getattr(realizations, name)

    @pytest.mark.parametrize("name", parameters.available)
    def test_validate_cosmo(self, name):
        """Test method ``validate`` for cosmology instance input."""
        cosmo = getattr(realizations, name)
        value = default_cosmology.validate(cosmo)
        assert value is cosmo


@pytest.mark.parametrize("name", parameters.available)
def test_pickle_builtin_realizations(name, pickle_protocol):
    """
    Test in-built realizations can pickle and unpickle.
    Also a regression test for #12008.
    """
    # get class instance
    original = getattr(cosmology, name)

    # pickle and unpickle
    f = pickle.dumps(original, protocol=pickle_protocol)
    with u.add_enabled_units(cu):
        unpickled = pickle.loads(f)

    assert unpickled == original
    assert unpickled.meta == original.meta

    # if the units are not enabled, it isn't equal because redshift units
    # are not equal. This is a weird, known issue.
    unpickled = pickle.loads(f)
    assert unpickled == original
    assert unpickled.meta != original.meta
