# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import random

import numpy as np
import pytest

from astropy.cosmology import Cosmology, w0wzCDM
from astropy.cosmology._src.io.builtin.model import (
    _CosmologyModel,
    from_model,
    to_model,
)
from astropy.cosmology._src.tests.helper import get_redshift_methods
from astropy.modeling.models import Gaussian1D
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .base import ToFromDirectTestBase, ToFromTestMixinBase

###############################################################################


class ToFromModelTestMixin(ToFromTestMixinBase):
    """Tests for a Cosmology[To/From]Format with ``format="astropy.model"``.

    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    @pytest.fixture(scope="class")
    def method_name(self, cosmo):
        # get methods, ignoring private and dunder
        methods = get_redshift_methods(cosmo, include_private=False, include_z2=True)

        # dynamically detect ABC and optional dependencies
        for n in tuple(methods):
            params = inspect.signature(getattr(cosmo, n)).parameters.keys()

            ERROR_SEIVE = (NotImplementedError, ValueError)
            #              # ABC                can't introspect for good input
            if not HAS_SCIPY:
                ERROR_SEIVE = ERROR_SEIVE + (ModuleNotFoundError,)

            args = np.arange(len(params)) + 1
            try:
                getattr(cosmo, n)(*args)
            except ERROR_SEIVE:
                methods.discard(n)
            except TypeError:
                # w0wzCDM has numerical instabilities when evaluating at z->inf
                # TODO: a more robust fix in w0wzCDM itself would be better
                if isinstance(cosmo, w0wzCDM):
                    methods.discard(n)

        # TODO! pytest doesn't currently allow multiple yields (`cosmo`) so
        # testing with 1 random method
        # yield from methods
        return random.choice(tuple(methods)) if methods else None

    # ===============================================================

    def test_fromformat_model_wrong_cls(self, from_format):
        """Test when Model is not the correct class."""
        model = Gaussian1D(amplitude=10, mean=14)

        with pytest.raises(AttributeError):
            from_format(model)

    def test_toformat_model_not_method(self, to_format):
        """Test when method is not a method."""
        with pytest.raises(AttributeError):
            to_format("astropy.model", method="this is definitely not a method.")

    def test_toformat_model_not_callable(self, to_format):
        """Test when method is actually an attribute."""
        with pytest.raises(ValueError):
            to_format("astropy.model", method="name")

    def test_toformat_model(self, cosmo, to_format, method_name):
        """Test cosmology -> astropy.model."""
        if method_name is None:  # no test if no method
            return

        model = to_format("astropy.model", method=method_name)
        assert isinstance(model, _CosmologyModel)

        # Parameters
        expect = tuple(k for k, v in cosmo.parameters.items() if v is not None)
        assert model.param_names == expect

        # scalar result
        args = np.arange(model.n_inputs) + 1

        got = model.evaluate(*args)
        expected = getattr(cosmo, method_name)(*args)
        assert np.all(got == expected)

        got = model(*args)
        expected = getattr(cosmo, method_name)(*args)
        np.testing.assert_allclose(got, expected)

        # vector result
        if "scalar" not in method_name:
            args = (np.ones((model.n_inputs, 3)).T + np.arange(model.n_inputs)).T

            got = model.evaluate(*args)
            expected = getattr(cosmo, method_name)(*args)
            assert np.all(got == expected)

            got = model(*args)
            expected = getattr(cosmo, method_name)(*args)
            np.testing.assert_allclose(got, expected)

    def test_tofromformat_model_instance(
        self, cosmo_cls, cosmo, method_name, to_format, from_format
    ):
        """Test cosmology -> astropy.model -> cosmology."""
        if method_name is None:  # no test if no method
            return

        # ------------
        # To Model
        # this also serves as a test of all added methods / attributes
        # in _CosmologyModel.

        model = to_format("astropy.model", method=method_name)

        assert isinstance(model, _CosmologyModel)
        assert model.cosmology_class is cosmo_cls
        assert model.cosmology == cosmo
        assert model.method_name == method_name

        # ------------
        # From Model

        # it won't error if everything matches up
        got = from_format(model, format="astropy.model")
        assert got == cosmo
        assert set(cosmo.meta.keys()).issubset(got.meta.keys())
        # Note: model adds parameter attributes to the metadata

        # also it auto-identifies 'format'
        got = from_format(model)
        assert got == cosmo
        assert set(cosmo.meta.keys()).issubset(got.meta.keys())

    def test_fromformat_model_subclass_partial_info(self) -> None:
        """
        Test writing from an instance and reading from that class.
        This works with missing information.

        There's no partial information with a Model
        """

    @pytest.mark.parametrize("format", [True, False, None, "astropy.model"])
    def test_is_equivalent_to_model(self, cosmo, method_name, to_format, format):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`.

        This test checks that Cosmology equivalency can be extended to any
        Python object that can be converted to a Cosmology -- in this case
        a model.
        """
        if method_name is None:  # no test if no method
            return

        obj = to_format("astropy.model", method=method_name)
        assert not isinstance(obj, Cosmology)

        is_equiv = cosmo.is_equivalent(obj, format=format)
        assert is_equiv is (format is not False)


class TestToFromModel(ToFromDirectTestBase, ToFromModelTestMixin):
    """Directly test ``to/from_model``."""

    def setup_class(self):
        self.functions = {"to": to_model, "from": from_model}
