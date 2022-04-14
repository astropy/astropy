# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for :mod:`astropy.cosmology.comparison`"""

import re

import numpy as np
import pytest

from astropy.cosmology import (Cosmology, FlatCosmologyMixin, Planck18, cosmology_equal,
                               cosmology_equivalent, cosmology_not_equal)
from astropy.cosmology.comparison import _CosmologyWrapper, _parse_format, _parse_formats
from astropy.cosmology.connect import convert_registry
from astropy.cosmology.io.tests.base import ToFromTestMixinBase


class ComparisonFunctionTestBase(ToFromTestMixinBase):
    """Tests for cosmology comparison functions.

    This class inherits from
    `astropy.cosmology.io.tests.base.ToFromTestMixinBase` because the cosmology
    comparison functions all have a kwarg ``format`` that allow the arguments
    to be converted to a |Cosmology| using the ``to_format`` architecture.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass.
    """

    @pytest.fixture(scope="class")
    def cosmo(self):
        return Planck18

    @pytest.fixture(scope="class",
    params={k for k, _ in convert_registry._readers.keys()} - {"astropy.cosmology"})
    def format(self, request):
        return request.param

    @pytest.fixture(scope="class")
    def xfail_cant_autoidentify(self, format):
        """`pytest.fixture` form of method ``can_autoidentify`."""
        if not self.can_autodentify(format):
            pytest.xfail("cannot autoidentify")

    @pytest.fixture(scope="class")
    def converted(self, to_format, format):
        if format == "astropy.model":  # special case Model
            return to_format(format, method="comoving_distance")
        return to_format(format)

    @pytest.fixture(scope="class")
    def pert_cosmo(self, cosmo):
        # change one parameter
        p = cosmo.__parameters__[0]
        cosmo2 = cosmo.clone(**{p: getattr(cosmo, p) * 1.0001})
        return cosmo2

    @pytest.fixture(scope="class")
    def pert_converted(self, pert_cosmo, format):
        if format == "astropy.model":  # special case Model
            return pert_cosmo.to_format(format, method="comoving_distance")
        return pert_cosmo.to_format(format)

    # ========================================================================

    def test_parse_format_error_wrong_format(self, cosmo):
        with pytest.raises(ValueError, match=re.escape("for a Cosmology, 'format'")):
            _parse_format(cosmo, "mapping")


class Test_parse_format(ComparisonFunctionTestBase):
    """Test functions ``_parse_format``."""

    @pytest.fixture(scope="class")
    def converted(self, to_format, format):
        if format == "astropy.model":  # special case Model
            return to_format(format, method="comoving_distance")

        converted = to_format(format)

        # Some raise a segfault! TODO: figure out why
        if isinstance(converted, _CosmologyWrapper._cantbroadcast):
            converted = _CosmologyWrapper(converted)

        return converted

    # ========================================================================

    def test_shortcut(self, cosmo):
        assert _parse_format(cosmo, None) == cosmo

    def test_convert(self, converted, format, cosmo):
        """Test converting a cosmology-like object"""
        out = _parse_format(converted, format)

        assert isinstance(out, Cosmology)
        assert out == cosmo

    def test_parse_format_error_noncosmology_cant_convert(self):
        """
        Test ``_parse_format`` errors when given a non-Cosmology object
        and format is False.
        """
        notacosmo = object()

        with pytest.raises(TypeError, match=re.escape("if 'format' is False")):
            _parse_format(notacosmo, False)

    def test_parse_format_vectorized(self, cosmo, format, converted):
        # vectorized on cosmos
        out = _parse_format([cosmo, cosmo], None)
        assert len(out) == 2
        assert np.all(out == cosmo)

        # vectorized on formats
        out = _parse_format(cosmo, [None, None])
        assert len(out) == 2
        assert np.all(out == cosmo)

        # more complex broadcast
        out = _parse_format([[cosmo, converted], [converted, cosmo]],
                             [[None, format], [format, None]])
        assert out.shape == (2, 2)
        assert np.all(out == cosmo)

    def test_parse_formats_vectorized(self, cosmo):
        # vectorized on cosmos
        out = _parse_formats(cosmo, cosmo, format=None)
        assert len(out) == 2
        assert np.all(out == cosmo)

        # does NOT vectorize on formats
        with pytest.raises(ValueError, match="operands could not be broadcast"):
            _parse_formats(cosmo, format=[None, None])


class Test_cosmology_equal(ComparisonFunctionTestBase):
    """Test :func:`astropy.cosmology.comparison.cosmology_equal`"""

    def test_cosmology_equal_simple(self, cosmo, pert_cosmo):
        # equality
        assert cosmology_equal(cosmo, cosmo) is True

        # not equal to perturbed cosmology
        assert cosmology_equal(cosmo, pert_cosmo) is False

    def test_cosmology_equal_format_error(self, cosmo, converted):
        # Not converting `converted`
        with pytest.raises(TypeError):
            cosmology_equal(cosmo, converted)

        with pytest.raises(TypeError):
            cosmology_equal(cosmo, converted, format=False)

    def test_cosmology_equal_format_auto(self, cosmo, converted, xfail_cant_autoidentify):
        # These tests only run if the format can autoidentify.
        assert cosmology_equal(cosmo, converted, format=None) is True
        assert cosmology_equal(cosmo, converted, format=True) is True

    def test_cosmology_equal_format_specify(self, cosmo, format, converted, pert_converted):
        # equality
        assert cosmology_equal(cosmo, converted, format=[None, format]) is True
        assert cosmology_equal(converted, cosmo, format=[format, None]) is True

        # non-equality
        assert cosmology_equal(cosmo, pert_converted, format=[None, format]) is False


class Test_cosmology_not_equal(ComparisonFunctionTestBase):
    """Test :func:`astropy.cosmology.comparison.cosmology_not_equal`"""

    def test_cosmology_not_equal_simple(self, cosmo, pert_cosmo):
        # equality
        assert cosmology_not_equal(cosmo, cosmo) is False

        # not equal to perturbed cosmology
        assert cosmology_not_equal(cosmo, pert_cosmo) is True

    def test_cosmology_not_equal_format_error(self, cosmo, converted):
        # Not converting `converted`
        with pytest.raises(TypeError):
            cosmology_not_equal(cosmo, converted)

        with pytest.raises(TypeError):
            cosmology_not_equal(cosmo, converted, format=False)

    def test_cosmology_not_equal_format_auto(self, cosmo, pert_converted, xfail_cant_autoidentify):
        assert cosmology_not_equal(cosmo, pert_converted, format=None) is True
        assert cosmology_not_equal(cosmo, pert_converted, format=True) is True

    def test_cosmology_not_equal_format_specify(self, cosmo, format, converted, pert_converted):
        # specifying the format
        assert cosmology_not_equal(cosmo, pert_converted, format=[None, format]) is True
        assert cosmology_not_equal(pert_converted, cosmo, format=[format, None]) is True

        # equality
        assert cosmology_not_equal(cosmo, converted, format=[None, format]) is False


class Test_cosmology_equivalent(ComparisonFunctionTestBase):
    """Test :func:`astropy.cosmology.comparison.cosmology_equivalent`"""

    @pytest.fixture
    def cosmo_eqvxflat(self, cosmo):
        if isinstance(cosmo, FlatCosmologyMixin):
            return cosmo.equivalent_nonflat

        pytest.skip("cosmology is not flat, "
                    "so does not have an equivalent non-flat cosmology.")

    @pytest.fixture
    def pert_cosmo_eqvxflat(self, pert_cosmo):
        if isinstance(pert_cosmo, FlatCosmologyMixin):
            return pert_cosmo.equivalent_nonflat

        pytest.skip("cosmology is not flat, "
                    "so does not have an equivalent non-flat cosmology.")

    # ========================================================================

    def test_cosmology_equivalent_equality(self, cosmo, pert_cosmo):
        # equality
        assert cosmology_equivalent(cosmo, cosmo) is True

        # not equal to perturbed cosmology
        assert cosmology_equivalent(cosmo, pert_cosmo) is False

    def test_cosmology_equivalent_equiv(self, cosmo, cosmo_eqvxflat,
                                        pert_cosmo, pert_cosmo_eqvxflat):
        # now need to check equivalent, but not equal, cosmologies.
        assert cosmology_equivalent(cosmo, cosmo.equivalent_nonflat) is True
        assert cosmology_equivalent(pert_cosmo, pert_cosmo.equivalent_nonflat) is True

    def test_cosmology_equivalent_format_error(self, cosmo, converted):
        # Not converting `converted`
        with pytest.raises(TypeError):
            cosmology_equivalent(cosmo, converted)

        with pytest.raises(TypeError):
            cosmology_equivalent(cosmo, converted, format=False)

    def test_cosmology_equivalent_format_auto(self, cosmo, converted, cosmo_eqvxflat, xfail_cant_autoidentify):
        assert cosmology_equivalent(cosmo, cosmo_eqvxflat, format=None) is True
        assert cosmology_equivalent(cosmo, cosmo_eqvxflat, format=True) is True

    def test_cosmology_equivalent_format_specify(self, cosmo, format, converted, cosmo_eqvxflat):
        # specifying the format
        assert cosmology_equivalent(cosmo_eqvxflat, converted, format=[None, format]) is True
        assert cosmology_equivalent(converted, cosmo_eqvxflat, format=[format, None]) is True

        # equality
        assert cosmology_equivalent(cosmo_eqvxflat, converted, format=[None, format]) is True
