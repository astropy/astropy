# -*- coding: utf-8 -*-

"""Test :mod:`~astropy.cosmology.builtin`."""

import warnings

from unittest import mock

import pytest

try:
    from pkg_resources import EntryPoint

    HAS_PKG = True
except ImportError:
    HAS_PKG = False

from .. import builtin, core
from astropy.utils.state import _StateProxy
from astropy.utils.exceptions import AstropyUserWarning


def test_parameter_registry():
    """
    Test `~astropy.cosmology.builtin._parameter_registry`.

    """
    registry = builtin._parameter_registry

    assert isinstance(registry, dict)

    minkeys = {"parameters", "references", "cosmo"}

    for n, state in registry.items():
        assert isinstance(state, (_StateProxy, dict))
        assert minkeys.issubset(state.keys())

        assert n in builtin.available


def test_default_cosmology():
    """
    Test `~astropy.cosmology.builtin.default_cosmology`.

    """
    # Attributes test

    assert builtin.default_cosmology._default_value == "Planck15"
    assert builtin.default_cosmology._registry == builtin._parameter_registry
    assert builtin.default_cosmology.available == builtin.available

    # ----------------------------------------
    # Get and Validate

    cosmo = builtin.default_cosmology.validate(None)
    assert cosmo == builtin.Planck15

    cosmo = builtin.default_cosmology.get_cosmology_from_string("no_default")
    assert cosmo is None

    builtin._parameter_registry.pop("test", None)  # just to make sure

    # not in registry
    with pytest.raises(ValueError):
        cosmo = builtin.default_cosmology.get_cosmology_from_string("test")

    # ----------------------------------------
    # Registration tests

    test_parameters = dict(
        Oc0=0.231,
        Ob0=0.0459,
        Om0=0.277,
        Ode0=0.4461,
        H0=70.2,
        n=0.962,
        sigma8=0.817,
        tau=0.088,
        z_reion=11.3,
        t0=13.72,
        Tcmb0=2.725,
        Neff=3.04,
        m_nu=0.0,
        flat=False,  # Needed to hit LambdaCDM, not only FlatLambdaCDM
    )

    try:  # register_parameters
        builtin.default_cosmology.register_parameters(
            name="test",
            parameters=test_parameters,
            references=None,
            cosmo=None,
            viewonly=True,
        )
        assert "test" in builtin._parameter_registry

        # test get_from_registry
        state = builtin.default_cosmology.get_from_registry("test")

        assert state["parameters"] == test_parameters
        assert state["references"] is None
        assert state["cosmo"] is None

        with builtin.default_cosmology.set("test"):
            test_cosmo = builtin.default_cosmology.get()

    finally:
        builtin._parameter_registry.pop("test", None)
        assert "test" not in builtin._parameter_registry

    # ----------------------------------------
    # Try failures
    test_parameters["flat"] = False
    test_parameters.pop("Ode0")

    try:  # register_parameters
        with pytest.raises(ValueError):
            builtin.default_cosmology.register_parameters(
                name="test", parameters=test_parameters
            )

    finally:
        builtin._parameter_registry.pop("test", None)
        assert "test" not in builtin._parameter_registry

    # ----------------------------------------

    try:  # register_cosmology_instance
        builtin.default_cosmology.register_cosmology_instance(cosmo=test_cosmo)

        assert "test" in builtin._parameter_registry
        assert builtin._parameter_registry["test"]["references"] is None
        assert builtin._parameter_registry["test"]["cosmo"] == core.LambdaCDM

        with builtin.default_cosmology.set("test"):
            test2_cosmo = builtin.default_cosmology.get()

        assert test2_cosmo == test_cosmo

    finally:
        builtin._parameter_registry.pop("test", None)
        assert "test" not in builtin._parameter_registry


@pytest.mark.skipif("not HAS_PKG")
class TestEntryPoint:
    """Tests population of fitting with entry point fitters"""

    def setup_class(self):
        self.exception_not_thrown = Exception(
            "The test should not have gotten here. There was no exception thrown"
        )

    def successfulimport(self):
        # This should work
        goodcosmo = core.FlatLambdaCDM(H0=68, Om0=0.3, name="test_cosmo")
        return goodcosmo

    def raiseimporterror(self):
        #  This should fail as it raises an Import Error
        raise ImportError

    def returnbadfunc(self):
        def badfunc():
            # This should import but it should fail type check
            pass

        return badfunc

    def returnbadclass(self):
        # This should import But it should fail subclass type check
        class badclass:
            pass

        return badclass

    def test_working(self):
        """This should work fine."""
        mock_entry_working = mock.create_autospec(EntryPoint)
        mock_entry_working.name = "Working"
        mock_entry_working.load = self.successfulimport
        builtin.populate_entry_points([mock_entry_working])

        assert "test_cosmo" in builtin._parameter_registry

        with builtin.default_cosmology.set("test_cosmo"):
            assert builtin.default_cosmology.get().name == "test_cosmo"

    def test_import_error(self):
        """This raises an import error on load to test that it is handled correctly"""
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                mock_entry_importerror = mock.create_autospec(EntryPoint)
                mock_entry_importerror.name = "IErr"
                mock_entry_importerror.load = self.raiseimporterror
                builtin.populate_entry_points([mock_entry_importerror])
            except AstropyUserWarning as w:
                if (
                    "ImportError" in w.args[0]
                ):  # any error for this case should have this in it.
                    pass
                else:
                    raise w
            else:
                raise self.exception_not_thrown

    def test_bad_class(self):
        """This returns a class which doesn't inherit from fitter """
        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                mock_entry_badclass = mock.create_autospec(EntryPoint)
                mock_entry_badclass.name = "BadClass"
                mock_entry_badclass.load = self.returnbadclass
                builtin.populate_entry_points([mock_entry_badclass])
            except AstropyUserWarning as w:
                if (
                    "Cosmology" in w.args[0]
                ):  # any error for this case should have this in it.
                    pass
                else:
                    raise w
            else:
                raise self.exception_not_thrown
