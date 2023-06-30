# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.core`."""

##############################################################################
# IMPORTS

# STDLIB
import abc
import inspect
import pickle

# THIRD PARTY
import numpy as np
import pytest

# LOCAL
import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import Cosmology, FlatCosmologyMixin
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.cosmology.parameter import Parameter
from astropy.table import Column, QTable, Table
from astropy.utils.compat import PYTHON_LT_3_11
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.metadata import MetaData

from .test_connect import ReadWriteTestMixin, ToFromFormatTestMixin
from .test_parameter import ParameterTestMixin

##############################################################################
# SETUP / TEARDOWN


def make_valid_zs(max_z: float = 1e5):
    """Make a list of valid redshifts for testing."""
    # scalar
    scalar_zs = [
        0,
        1,
        min(1100, max_z),  # interesting times
        # FIXME! np.inf breaks some funcs. 0 * inf is an error
        np.float64(min(3300, max_z)),  # different type
        2 * cu.redshift,
        3 * u.one,  # compatible units
    ]
    # array
    _zarr = np.linspace(0, min(1e5, max_z), num=20)
    array_zs = [
        _zarr,  # numpy
        _zarr.tolist(),  # pure python
        Column(_zarr),  # table-like
        _zarr * cu.redshift,  # Quantity
    ]
    return scalar_zs, _zarr, array_zs, scalar_zs + array_zs


scalar_zs, z_arr, array_zs, valid_zs = make_valid_zs()

invalid_zs = [
    (None, TypeError),  # wrong type
    # Wrong units (the TypeError is for the cython, which can differ)
    (4 * u.MeV, (u.UnitConversionError, TypeError)),  # scalar
    ([0, 1] * u.m, (u.UnitConversionError, TypeError)),  # array
]


class SubCosmology(Cosmology):
    """Defined here to be serializable."""

    H0 = Parameter(unit="km/(s Mpc)")
    Tcmb0 = Parameter(unit=u.K)
    m_nu = Parameter(unit=u.eV)

    def __init__(self, H0, Tcmb0=0 * u.K, m_nu=0 * u.eV, name=None, meta=None):
        super().__init__(name=name, meta=meta)
        self.H0 = H0
        self.Tcmb0 = Tcmb0
        self.m_nu = m_nu

    @property
    def is_flat(self):
        return super().is_flat()


##############################################################################
# TESTS
##############################################################################


class MetaTestMixin:
    """Tests for a :class:`astropy.utils.metadata.MetaData` on a Cosmology."""

    def test_meta_on_class(self, cosmo_cls):
        assert isinstance(cosmo_cls.meta, MetaData)

    def test_meta_on_instance(self, cosmo):
        assert isinstance(cosmo.meta, dict)  # test type
        # value set at initialization
        assert cosmo.meta == self.cls_kwargs.get("meta", {})

    def test_meta_mutable(self, cosmo):
        """The metadata is NOT immutable on a cosmology"""
        key = tuple(cosmo.meta.keys())[0]  # select some key
        cosmo.meta[key] = cosmo.meta.pop(key)  # will error if immutable


class CosmologyTest(
    ParameterTestMixin,
    MetaTestMixin,
    ReadWriteTestMixin,
    ToFromFormatTestMixin,
    metaclass=abc.ABCMeta,
):
    """
    Test subclasses of :class:`astropy.cosmology.Cosmology`.
    """

    @abc.abstractmethod
    def setup_class(self):
        """Setup for testing."""

    def teardown_class(self):
        pass

    @property
    def cls_args(self):
        return tuple(self._cls_args.values())

    @pytest.fixture(scope="class")
    def cosmo_cls(self):
        """The Cosmology class as a :func:`pytest.fixture`."""
        return self.cls

    @pytest.fixture(scope="function")  # ensure not cached.
    def ba(self):
        """Return filled `inspect.BoundArguments` for cosmology."""
        ba = self.cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)
        ba.apply_defaults()
        return ba

    @pytest.fixture(scope="class")
    def cosmo(self, cosmo_cls):
        """The cosmology instance with which to test."""
        ba = self.cls._init_signature.bind(*self.cls_args, **self.cls_kwargs)
        ba.apply_defaults()
        return cosmo_cls(*ba.args, **ba.kwargs)

    # ===============================================================
    # Method & Attribute Tests

    # ---------------------------------------------------------------
    # class-level

    def test_init_subclass(self, cosmo_cls):
        """Test creating subclasses registers classes and manages Parameters."""

        class InitSubclassTest(cosmo_cls):
            pass

        # test parameters
        assert InitSubclassTest.__parameters__ == cosmo_cls.__parameters__

        # test and cleanup registry
        registrant = _COSMOLOGY_CLASSES.pop(InitSubclassTest.__qualname__)
        assert registrant is InitSubclassTest

    def test_init_signature(self, cosmo_cls, cosmo):
        """Test class-property ``_init_signature``."""
        # test presence
        assert hasattr(cosmo_cls, "_init_signature")
        assert hasattr(cosmo, "_init_signature")

        # test internal consistency, so following tests can use either cls or instance.
        assert cosmo_cls._init_signature == cosmo._init_signature

        # test matches __init__, but without 'self'
        sig = inspect.signature(cosmo.__init__)  # (instances don't have self)
        assert set(sig.parameters.keys()) == set(
            cosmo._init_signature.parameters.keys()
        )
        assert all(
            np.all(sig.parameters[k].default == p.default)
            for k, p in cosmo._init_signature.parameters.items()
        )

    # ---------------------------------------------------------------
    # instance-level

    def test_init(self, cosmo_cls):
        """Test initialization."""
        # Cosmology only does name and meta, but this subclass adds H0 & Tcmb0.
        cosmo = cosmo_cls(*self.cls_args, name="test_init", meta={"m": 1})
        assert cosmo.name == "test_init"
        assert cosmo.meta["m"] == 1

        # if meta is None, it is changed to a dict
        cosmo = cosmo_cls(*self.cls_args, name="test_init", meta=None)
        assert cosmo.meta == {}

    def test_name(self, cosmo):
        """Test property ``name``."""
        assert cosmo.name is cosmo._name  # accesses private attribute
        assert cosmo.name is None or isinstance(cosmo.name, str)  # type
        assert cosmo.name == self.cls_kwargs["name"]  # test has expected value

        # immutable
        match = (
            "can't set"
            if PYTHON_LT_3_11
            else f"property 'name' of {cosmo.__class__.__name__!r} object has no setter"
        )
        with pytest.raises(AttributeError, match=match):
            cosmo.name = None

    @abc.abstractmethod
    def test_is_flat(self, cosmo_cls, cosmo):
        """Test property ``is_flat``."""

    # ------------------------------------------------
    # clone

    def test_clone_identical(self, cosmo):
        """Test method ``.clone()`` if no (kw)args."""
        assert cosmo.clone() is cosmo

    def test_clone_name(self, cosmo):
        """Test method ``.clone()`` name argument."""
        # test changing name. clone treats 'name' differently (see next test)
        c = cosmo.clone(name="cloned cosmo")
        assert c.name == "cloned cosmo"  # changed
        # show name is the only thing changed
        c._name = cosmo.name  # first change name back
        assert c == cosmo
        assert c.meta == cosmo.meta

        # now change a different parameter and see how 'name' changes
        c = cosmo.clone(meta={"test_clone_name": True})
        assert c.name == cosmo.name + " (modified)"

    def test_clone_meta(self, cosmo):
        """Test method ``.clone()`` meta argument: updates meta, doesn't clear."""
        # start with no change
        c = cosmo.clone(meta=None)
        assert c.meta == cosmo.meta

        # add something
        c = cosmo.clone(meta=dict(test_clone_meta=True))
        assert c.meta["test_clone_meta"] is True
        c.meta.pop("test_clone_meta")  # remove from meta
        assert c.meta == cosmo.meta  # now they match

    def test_clone_change_param(self, cosmo):
        """
        Test method ``.clone()`` changing a(many) Parameter(s).
        Nothing here b/c no Parameters.
        """

    def test_clone_fail_unexpected_arg(self, cosmo):
        """Test when ``.clone()`` gets an unexpected argument."""
        with pytest.raises(TypeError, match="unexpected keyword argument"):
            cosmo.clone(not_an_arg=4)

    def test_clone_fail_positional_arg(self, cosmo):
        with pytest.raises(TypeError, match="1 positional argument"):
            cosmo.clone(None)

    # ---------------------------------------------------------------
    # comparison methods

    def test_is_equivalent(self, cosmo):
        """Test :meth:`astropy.cosmology.Cosmology.is_equivalent`."""
        # to self
        assert cosmo.is_equivalent(cosmo)

        # same class, different instance
        newclone = cosmo.clone(name="test_is_equivalent")
        assert cosmo.is_equivalent(newclone)
        assert newclone.is_equivalent(cosmo)

        # different class and not convertible to Cosmology.
        assert not cosmo.is_equivalent(2)

    def test_equality(self, cosmo):
        """Test method ``.__eq__()."""
        # wrong class
        assert (cosmo != 2) and (2 != cosmo)
        # correct
        assert cosmo == cosmo
        # different name <= not equal, but equivalent
        newcosmo = cosmo.clone(name="test_equality")
        assert (cosmo != newcosmo) and (newcosmo != cosmo)
        assert cosmo.__equiv__(newcosmo) and newcosmo.__equiv__(cosmo)

    # ---------------------------------------------------------------

    def test_repr(self, cosmo_cls, cosmo):
        """Test method ``.__repr__()``.

        This is a very general test and it is probably good to have a
        hard-coded comparison.
        """
        r = repr(cosmo)

        # class in string rep
        assert cosmo_cls.__qualname__ in r
        assert r.index(cosmo_cls.__qualname__) == 0  # it's the first thing
        r = r[len(cosmo_cls.__qualname__) + 1 :]  # remove

        # name in string rep
        if cosmo.name is not None:
            assert f'name="{cosmo.name}"' in r
            assert r.index("name=") == 0
            r = r[6 + len(cosmo.name) + 3 :]  # remove

        # parameters in string rep
        ps = {k: getattr(cosmo, k) for k in cosmo.__parameters__}
        for k, v in ps.items():
            sv = f"{k}={v}"
            assert sv in r
            assert r.index(k) == 0
            r = r[len(sv) + 2 :]  # remove

    # ------------------------------------------------

    @pytest.mark.parametrize("in_meta", [True, False])
    @pytest.mark.parametrize("table_cls", [Table, QTable])
    def test_astropy_table(self, cosmo, table_cls, in_meta):
        """Test ``astropy.table.Table(cosmology)``."""
        tbl = table_cls(cosmo, cosmology_in_meta=in_meta)

        assert isinstance(tbl, table_cls)
        # the name & all parameters are columns
        for n in ("name", *cosmo.__parameters__):
            assert n in tbl.colnames
            assert np.all(tbl[n] == getattr(cosmo, n))
        # check if Cosmology is in metadata or a column
        if in_meta:
            assert tbl.meta["cosmology"] == cosmo.__class__.__qualname__
            assert "cosmology" not in tbl.colnames
        else:
            assert "cosmology" not in tbl.meta
            assert tbl["cosmology"][0] == cosmo.__class__.__qualname__
        # the metadata is transferred
        for k, v in cosmo.meta.items():
            assert np.all(tbl.meta[k] == v)

    # ===============================================================
    # Usage Tests

    def test_immutability(self, cosmo):
        """
        Test immutability of cosmologies.
        The metadata is mutable: see ``test_meta_mutable``.
        """
        for n in cosmo.__all_parameters__:
            with pytest.raises(AttributeError):
                setattr(cosmo, n, getattr(cosmo, n))

    def test_pickle_class(self, cosmo_cls, pickle_protocol):
        """Test classes can pickle and unpickle."""
        # pickle and unpickle
        f = pickle.dumps(cosmo_cls, protocol=pickle_protocol)
        unpickled = pickle.loads(f)

        # test equality
        assert unpickled == cosmo_cls

    def test_pickle_instance(self, cosmo, pickle_protocol):
        """Test instances can pickle and unpickle."""
        # pickle and unpickle
        f = pickle.dumps(cosmo, protocol=pickle_protocol)
        with u.add_enabled_units(cu):
            unpickled = pickle.loads(f)

        assert unpickled == cosmo
        assert unpickled.meta == cosmo.meta


class TestCosmology(CosmologyTest):
    """Test :class:`astropy.cosmology.Cosmology`.

    Subclasses should define tests for:

    - ``test_clone_change_param()``
    - ``test_repr()``
    """

    def setup_class(self):
        """
        Setup for testing.
        Cosmology should not be instantiated, so tests are done on a subclass.
        """
        # make sure SubCosmology is known
        _COSMOLOGY_CLASSES["SubCosmology"] = SubCosmology

        self.cls = SubCosmology
        self._cls_args = dict(
            H0=70 * (u.km / u.s / u.Mpc), Tcmb0=2.7 * u.K, m_nu=0.6 * u.eV
        )
        self.cls_kwargs = dict(name=self.__class__.__name__, meta={"a": "b"})

    def teardown_class(self):
        """Teardown for testing."""
        super().teardown_class(self)
        _COSMOLOGY_CLASSES.pop("SubCosmology", None)

    # ===============================================================
    # Method & Attribute Tests

    def test_is_flat(self, cosmo_cls, cosmo):
        """Test property ``is_flat``. It's an ABC."""
        with pytest.raises(NotImplementedError, match="is_flat is not implemented"):
            cosmo.is_flat


# -----------------------------------------------------------------------------


class FlatCosmologyMixinTest:
    """Tests for :class:`astropy.cosmology.core.FlatCosmologyMixin` subclasses.

    The test suite structure mirrors the implementation of the tested code.
    Just like :class:`astropy.cosmology.FlatCosmologyMixin` is an abstract
    base class (ABC) that cannot be used by itself, so too is this corresponding
    test class an ABC mixin.

    E.g to use this class::

        class TestFlatSomeCosmology(FlatCosmologyMixinTest, TestSomeCosmology):
            ...
    """

    def test_nonflat_class_(self, cosmo_cls, cosmo):
        """Test :attr:`astropy.cosmology.core.FlatCosmologyMixin.nonflat_cls`."""
        # Test it's a method on the class
        assert issubclass(cosmo_cls, cosmo_cls.__nonflatclass__)

        # It also works from the instance. # TODO! as a "metaclassmethod"
        assert issubclass(cosmo_cls, cosmo.__nonflatclass__)

        # Maybe not the most robust test, but so far all Flat classes have the
        # name of their parent class.
        assert cosmo.__nonflatclass__.__name__ in cosmo_cls.__name__

    def test_is_flat(self, cosmo_cls, cosmo):
        """Test property ``is_flat``."""
        super().test_is_flat(cosmo_cls, cosmo)

        # it's always True
        assert cosmo.is_flat is True

    def test_nonflat(self, cosmo):
        """Test :attr:`astropy.cosmology.core.FlatCosmologyMixin.nonflat`."""
        assert cosmo.nonflat.is_equivalent(cosmo)
        assert cosmo.is_equivalent(cosmo.nonflat)

    # ------------------------------------------------
    # clone

    def test_clone_to_nonflat_equivalent(self, cosmo):
        """Test method ``.clone()``to_nonflat argument."""
        # just converting the class
        nc = cosmo.clone(to_nonflat=True)
        assert isinstance(nc, cosmo.__nonflatclass__)
        assert nc == cosmo.nonflat

    @abc.abstractmethod
    def test_clone_to_nonflat_change_param(self, cosmo):
        """
        Test method ``.clone()`` changing a(many) Parameter(s). No parameters
        are changed here because FlatCosmologyMixin has no Parameters.
        See class docstring for why this test method exists.
        """
        # send to non-flat
        nc = cosmo.clone(to_nonflat=True)
        assert isinstance(nc, cosmo.__nonflatclass__)
        assert nc == cosmo.nonflat

    # ------------------------------------------------

    def test_is_equivalent(self, cosmo):
        """Test :meth:`astropy.cosmology.core.FlatCosmologyMixin.is_equivalent`.

        Normally this would pass up via super(), but ``__equiv__`` is meant
        to be overridden, so we skip super().
        e.g. FlatFLRWMixinTest -> FlatCosmologyMixinTest -> TestCosmology
        vs   FlatFLRWMixinTest -> FlatCosmologyMixinTest -> TestFLRW -> TestCosmology
        """
        CosmologyTest.test_is_equivalent(self, cosmo)

        # See FlatFLRWMixinTest for tests. It's a bit hard here since this class
        # is for an ABC.

    # ===============================================================
    # Usage Tests

    def test_subclassing(self, cosmo_cls):
        """Test when subclassing a flat cosmology."""

        class SubClass1(cosmo_cls):
            pass

        # The classes have the same non-flat parent class
        assert SubClass1.__nonflatclass__ is cosmo_cls.__nonflatclass__

        # A more complex example is when Mixin classes are used.
        class Mixin:
            pass

        class SubClass2(Mixin, cosmo_cls):
            pass

        # The classes have the same non-flat parent class
        assert SubClass2.__nonflatclass__ is cosmo_cls.__nonflatclass__

        # The order of the Mixin should not matter
        class SubClass3(cosmo_cls, Mixin):
            pass

        # The classes have the same non-flat parent class
        assert SubClass3.__nonflatclass__ is cosmo_cls.__nonflatclass__


def test__nonflatclass__multiple_nonflat_inheritance():
    """
    Test :meth:`astropy.cosmology.core.FlatCosmologyMixin.__nonflatclass__`
    when there's more than one non-flat class in the inheritance.
    """

    # Define a non-operable minimal subclass of Cosmology.
    class SubCosmology2(Cosmology):
        def __init__(self, H0, Tcmb0=0 * u.K, m_nu=0 * u.eV, name=None, meta=None):
            super().__init__(name=name, meta=meta)

        @property
        def is_flat(self):
            return False

    # Now make an ambiguous flat cosmology from the two SubCosmologies
    with pytest.raises(TypeError, match="cannot create a consistent non-flat class"):

        class FlatSubCosmology(FlatCosmologyMixin, SubCosmology, SubCosmology2):
            @property
            def nonflat(self):
                pass


# -----------------------------------------------------------------------------


def test_flrw_moved_deprecation():
    """Test the deprecation warning about the move of FLRW classes."""
    from astropy.cosmology import flrw

    # it's deprecated to import `flrw/*` from `core.py`
    with pytest.warns(AstropyDeprecationWarning):
        from astropy.cosmology.core import FLRW

    # but they are the same object
    assert FLRW is flrw.FLRW
