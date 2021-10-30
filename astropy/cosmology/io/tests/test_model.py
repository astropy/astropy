# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import inspect
import random

# THIRD PARTY
import pytest

import numpy as np

# LOCAL
import astropy.units as u
from astropy import cosmology
from astropy.cosmology import Cosmology, Planck18, realizations
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Parameter
from astropy.cosmology.io.model import _CosmologyModel, from_model, to_model
from astropy.cosmology.parameters import available
from astropy.modeling import FittableModel
from astropy.modeling.models import Gaussian1D
from astropy.utils.compat.optional_deps import HAS_SCIPY

from .base import IOTestMixinBase, IOFormatTestBase

cosmo_instances = [getattr(realizations, name) for name in available]
cosmo_instances.append("TestToFromTable.setup.<locals>.CosmologyWithKwargs")


###############################################################################


class ToFromModelTestMixin(IOTestMixinBase):
    """
    Tests for a Cosmology[To/From]Format with ``format="astropy.model"``.
    This class will not be directly called by :mod:`pytest` since its name does
    not begin with ``Test``. To activate the contained tests this class must
    be inherited in a subclass. Subclasses must define a :func:`pytest.fixture`
    ``cosmo`` that returns/yields an instance of a |Cosmology|.
    See ``TestCosmologyToFromFormat`` or ``TestCosmology`` for examples.
    """

    @pytest.fixture
    def method_name(self, cosmo):
        # get methods, ignoring private and dunder
        methods = {n for n in dir(cosmo)
                   if (callable(getattr(cosmo, n)) and not n.startswith("_"))}
        # sieve out incompatible methods
        for n in tuple(methods):
            # remove non-introspectable methods
            try:
                sig = inspect.signature(getattr(cosmo, n))
            except ValueError:
                methods.discard(n)
                continue

            params = list(sig.parameters.keys())
            # remove non redshift methods
            if len(params) == 0 or not params[0].startswith("z"):
                methods.discard(n)
                continue

            # dynamically detect ABC and optional dependencies
            ERROR_SEIVE = (NotImplementedError, ValueError)
            #              # ABC                can't introspect for good input
            if not HAS_SCIPY:
                ERROR_SEIVE = ERROR_SEIVE + (ModuleNotFoundError, )

            args = np.arange(len(params)) + 1
            try:
                getattr(cosmo, n)(*args)
            except ERROR_SEIVE:
                methods.discard(n)

        # TODO! pytest doesn't currently allow multiple yields (`cosmo`) so
        # testing with 1 random method
        # yield from methods
        return random.choice(tuple(methods)) if methods else None

    # ===============================================================

    def test_fromformat_model_wrong_cls(self, from_format):
        """Test when Model is not the correct class."""
        model = Gaussian1D(amplitude=10, mean=14)

        with pytest.raises(TypeError, match="`model` must be"):
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
        expect = tuple([n for n in cosmo.__parameters__ if getattr(cosmo, n) is not None])
        assert model.param_names == expect

        # scalar result
        args = np.arange(model.n_inputs) + 1

        got = model.evaluate(*args)
        expected = getattr(cosmo, method_name)(*args)
        assert np.all(got == expected)

        got = model(*args)
        expected = getattr(cosmo, method_name)(*args)
        assert np.all(got == expected)

        # vector result
        if "scalar" not in method_name:

            args = (np.ones((model.n_inputs, 3)).T + np.arange(model.n_inputs)).T

            got = model.evaluate(*args)
            expected = getattr(cosmo, method_name)(*args)
            assert np.all(got == expected)

            got = model(*args)
            expected = getattr(cosmo, method_name)(*args)
            assert np.all(got == expected)

    def test_tofromformat_model_instance(self, cosmo_cls, cosmo, method_name,
                                         to_format, from_format):
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

    def test_fromformat_model_subclass_partial_info(self):
        """
        Test writing from an instance and reading from that class.
        This works with missing information.
        """
        pass  # there's no partial information with a Model


class TestToFromModel(IOFormatTestBase, ToFromModelTestMixin):
    """Directly test ``to/from_model``."""

    def setup_class(self):
        self.functions = {"to": to_model, "from": from_model}
