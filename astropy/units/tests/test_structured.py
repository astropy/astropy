# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities.
"""

import copy

import numpy as np
import numpy.lib.recfunctions as rfn
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.tests.helper import check_pickling_recovery, pickle_protocol  # noqa: F401
from astropy.units import Quantity, StructuredUnit, Unit, UnitBase
from astropy.units.quantity import _structured_unit_like_dtype
from astropy.utils.masked import Masked


class StructuredTestBase:
    @classmethod
    def setup_class(cls):
        cls.pv_dtype = np.dtype([("p", "f8"), ("v", "f8")])
        cls.pv_t_dtype = np.dtype([("pv", cls.pv_dtype), ("t", "f8")])
        cls.p_unit = u.km
        cls.v_unit = u.km / u.s
        cls.t_unit = u.s
        cls.pv_dtype = np.dtype([("p", "f8"), ("v", "f8")])
        cls.pv_t_dtype = np.dtype([("pv", cls.pv_dtype), ("t", "f8")])
        cls.pv = np.array([(1.0, 0.25), (2.0, 0.5), (3.0, 0.75)], cls.pv_dtype)
        cls.pv_t = np.array(
            [
                ((4.0, 2.5), 0.0),
                ((5.0, 5.0), 1.0),
                ((6.0, 7.5), 2.0),
            ],
            cls.pv_t_dtype,
        )


class StructuredTestBaseWithUnits(StructuredTestBase):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.pv_unit = StructuredUnit((cls.p_unit, cls.v_unit), ("p", "v"))
        cls.pv_t_unit = StructuredUnit((cls.pv_unit, cls.t_unit), ("pv", "t"))


class TestStructuredUnitBasics(StructuredTestBase):
    def test_initialization_and_keying(self):
        su = StructuredUnit((self.p_unit, self.v_unit), ("p", "v"))
        assert su["p"] is self.p_unit
        assert su["v"] is self.v_unit
        su2 = StructuredUnit((su, self.t_unit), ("pv", "t"))
        assert isinstance(su2["pv"], StructuredUnit)
        assert su2["pv"]["p"] is self.p_unit
        assert su2["pv"]["v"] is self.v_unit
        assert su2["t"] is self.t_unit
        assert su2["pv"] == su
        su3 = StructuredUnit(("AU", "AU/day"), ("p", "v"))
        assert isinstance(su3["p"], UnitBase)
        assert isinstance(su3["v"], UnitBase)
        su4 = StructuredUnit("AU, AU/day", ("p", "v"))
        assert su4["p"] == u.AU
        assert su4["v"] == u.AU / u.day
        su5 = StructuredUnit(("AU", "AU/day"))
        assert su5.field_names == ("f0", "f1")
        assert su5["f0"] == u.AU
        assert su5["f1"] == u.AU / u.day

    def test_recursive_initialization(self):
        su = StructuredUnit(
            ((self.p_unit, self.v_unit), self.t_unit), (("p", "v"), "t")
        )
        assert isinstance(su["pv"], StructuredUnit)
        assert su["pv"]["p"] is self.p_unit
        assert su["pv"]["v"] is self.v_unit
        assert su["t"] is self.t_unit
        su2 = StructuredUnit(
            ((self.p_unit, self.v_unit), self.t_unit), (["p_v", ("p", "v")], "t")
        )
        assert isinstance(su2["p_v"], StructuredUnit)
        assert su2["p_v"]["p"] is self.p_unit
        assert su2["p_v"]["v"] is self.v_unit
        assert su2["t"] is self.t_unit
        su3 = StructuredUnit((("AU", "AU/day"), "yr"), (["p_v", ("p", "v")], "t"))
        assert isinstance(su3["p_v"], StructuredUnit)
        assert su3["p_v"]["p"] == u.AU
        assert su3["p_v"]["v"] == u.AU / u.day
        assert su3["t"] == u.yr
        su4 = StructuredUnit("(AU, AU/day), yr", (("p", "v"), "t"))
        assert isinstance(su4["pv"], StructuredUnit)
        assert su4["pv"]["p"] == u.AU
        assert su4["pv"]["v"] == u.AU / u.day
        assert su4["t"] == u.yr

    def test_extreme_recursive_initialization(self):
        su = StructuredUnit(
            "(yr,(AU,AU/day,(km,(day,day))),m)",
            ("t", ("p", "v", ("h", ("d1", "d2"))), "l"),
        )
        assert su.field_names == (
            't', ['pvhd1d2',
                  ('p', 'v',
                   ['hd1d2',
                    ('h',
                     ['d1d2',
                      ('d1', 'd2')])])],
            'l',
        )  # fmt: skip

        dt = np.dtype(
            [("t", "f8"),
             ("pvhd1d2",
              ([("p", "f8"), ("v", "f8"), ("hd1d2",
                                           [("h", "f8"), ("d1d2",
                                                          [("d1", "f8"), ("d2", "f8")]),
                                            ]),
                ], (5, 5))),  # Note: structured subarray to improve test!
             ("l", "f8")
             ])  # fmt: skip

        su2 = StructuredUnit("(yr,(AU,AU/day,(km,(day,day))),m)", dt)
        assert su2.field_names == su.field_names
        assert su2 == su

    @pytest.mark.parametrize(
        "names, invalid",
        [
            [("t", ["p", "v"]), "['p', 'v']"],
            [("t", ["pv", "p", "v"]), "['pv', 'p', 'v']"],
            [("t", ["pv", ["p", "v"]]), "['pv', ['p', 'v']"],
            [("t", ()), "()"],
            [("t", ("p", None)), "None"],
            [("t", ["pv", ("p", "")]), "''"],
        ],
    )
    def test_initialization_names_invalid_list_errors(self, names, invalid):
        with pytest.raises(ValueError) as exc:
            StructuredUnit("yr,(AU,AU/day)", names)
        assert f"invalid entry {invalid}" in str(exc)

    def test_looks_like_unit(self):
        su = StructuredUnit((self.p_unit, self.v_unit), ("p", "v"))
        assert Unit(su) is su

    def test_initialize_with_float_dtype(self):
        su = StructuredUnit(("AU", "AU/d"), self.pv_dtype)
        assert isinstance(su["p"], UnitBase)
        assert isinstance(su["v"], UnitBase)
        assert su["p"] == u.AU
        assert su["v"] == u.AU / u.day
        su = StructuredUnit((("km", "km/s"), "yr"), self.pv_t_dtype)
        assert isinstance(su["pv"], StructuredUnit)
        assert isinstance(su["pv"]["p"], UnitBase)
        assert isinstance(su["t"], UnitBase)
        assert su["pv"]["v"] == u.km / u.s
        su = StructuredUnit("(km, km/s), yr", self.pv_t_dtype)
        assert isinstance(su["pv"], StructuredUnit)
        assert isinstance(su["pv"]["p"], UnitBase)
        assert isinstance(su["t"], UnitBase)
        assert su["pv"]["v"] == u.km / u.s

    def test_initialize_with_structured_unit_for_names(self):
        su = StructuredUnit(("AU", "AU/d"), names=("p", "v"))
        su2 = StructuredUnit(("km", "km/s"), names=su)
        assert su2.field_names == ("p", "v")
        assert su2["p"] == u.km
        assert su2["v"] == u.km / u.s

    def test_initialize_single_field(self):
        su = StructuredUnit("AU", "p")
        assert isinstance(su, StructuredUnit)
        assert isinstance(su["p"], UnitBase)
        assert su["p"] == u.AU
        su = StructuredUnit("AU")
        assert isinstance(su, StructuredUnit)
        assert isinstance(su["f0"], UnitBase)
        assert su["f0"] == u.AU

    def test_equality(self):
        su = StructuredUnit(("AU", "AU/d"), self.pv_dtype)
        assert su == StructuredUnit(("AU", "AU/d"), self.pv_dtype)
        assert su != StructuredUnit(("m", "AU/d"), self.pv_dtype)
        # Names should be ignored.
        assert su == StructuredUnit(("AU", "AU/d"))
        assert su == StructuredUnit(("AU", "AU/d"), names=("q", "w"))
        assert su != StructuredUnit(("m", "m/s"))

    def test_parsing(self):
        su = Unit("AU, AU/d")
        assert isinstance(su, StructuredUnit)
        assert isinstance(su["f0"], UnitBase)
        assert isinstance(su["f1"], UnitBase)
        assert su["f0"] == u.AU
        assert su["f1"] == u.AU / u.day
        su2 = Unit("AU, AU/d, yr")
        assert isinstance(su2, StructuredUnit)
        assert su2 == StructuredUnit(("AU", "AU/d", "yr"))
        su2a = Unit("(AU, AU/d, yr)")
        assert isinstance(su2a, StructuredUnit)
        assert su2a == su2
        su3 = Unit("(km, km/s), yr")
        assert isinstance(su3, StructuredUnit)
        assert su3 == StructuredUnit((("km", "km/s"), "yr"))
        su4 = Unit("km,")
        assert isinstance(su4, StructuredUnit)
        assert su4 == StructuredUnit((u.km,))
        su5 = Unit("(m,s),")
        assert isinstance(su5, StructuredUnit)
        assert su5 == StructuredUnit(((u.m, u.s),))
        ldbody_unit = Unit("Msun, 0.5rad^2, (au, au/day)")
        assert ldbody_unit == StructuredUnit(
            (u.Msun, Unit(u.rad**2 / 2), (u.AU, u.AU / u.day))
        )

    def test_to_string(self):
        su = StructuredUnit((u.km, u.km / u.s))
        latex_str = r"$(\mathrm{km}, \mathrm{\frac{km}{s}})$"
        assert su.to_string(format="latex") == latex_str
        latex_str = r"$(\mathrm{km}, \mathrm{km\,s^{-1}})$"
        assert su.to_string(format="latex_inline") == latex_str

    def test_str(self):
        su = StructuredUnit(((u.km, u.km / u.s), u.yr))
        assert str(su) == "((km, km / s), yr)"
        assert Unit(str(su)) == su

    def test_repr(self):
        su = StructuredUnit(((u.km, u.km / u.s), u.yr))
        assert repr(su) == 'Unit("((km, km / s), yr)")'
        assert eval(repr(su)) == su


class TestStructuredUnitsCopyPickle(StructuredTestBaseWithUnits):
    def test_copy(self):
        su_copy = copy.copy(self.pv_t_unit)
        assert su_copy is not self.pv_t_unit
        assert su_copy == self.pv_t_unit
        assert su_copy._units is self.pv_t_unit._units

    def test_deepcopy(self):
        su_copy = copy.deepcopy(self.pv_t_unit)
        assert su_copy is not self.pv_t_unit
        assert su_copy == self.pv_t_unit
        assert su_copy._units is not self.pv_t_unit._units

    def test_pickle(self, pickle_protocol):  # noqa: F811
        check_pickling_recovery(self.pv_t_unit, pickle_protocol)


class TestStructuredUnitAsMapping(StructuredTestBaseWithUnits):
    def test_len(self):
        assert len(self.pv_unit) == 2
        assert len(self.pv_t_unit) == 2

    def test_keys(self):
        slv = list(self.pv_t_unit.keys())
        assert slv == ["pv", "t"]

    def test_values(self):
        values = self.pv_t_unit.values()
        assert values == (self.pv_unit, self.t_unit)

    def test_field_names(self):
        field_names = self.pv_t_unit.field_names
        assert isinstance(field_names, tuple)
        assert field_names == (["pv", ("p", "v")], "t")

    @pytest.mark.parametrize("iterable", [list, set])
    def test_as_iterable(self, iterable):
        sl = iterable(self.pv_unit)
        assert isinstance(sl, iterable)
        assert sl == iterable(["p", "v"])

    def test_as_dict(self):
        sd = dict(self.pv_t_unit)
        assert sd == {"pv": self.pv_unit, "t": self.t_unit}

    def test_contains(self):
        assert "p" in self.pv_unit
        assert "v" in self.pv_unit
        assert "t" not in self.pv_unit

    def test_setitem_fails(self):
        with pytest.raises(TypeError, match="item assignment"):
            self.pv_t_unit["t"] = u.Gyr


class TestStructuredUnitMethods(StructuredTestBaseWithUnits):
    def test_physical_type_id(self):
        pv_ptid = self.pv_unit._physical_type_id
        assert len(pv_ptid) == 2
        assert pv_ptid.dtype.names == ("p", "v")
        p_ptid = self.pv_unit["p"]._physical_type_id
        v_ptid = self.pv_unit["v"]._physical_type_id
        # Expected should be (subclass of) void, with structured object dtype.
        expected = np.array((p_ptid, v_ptid), [("p", "O"), ("v", "O")])[()]
        assert pv_ptid == expected
        # Names should be ignored in comparison.
        assert pv_ptid == np.array((p_ptid, v_ptid), "O,O")[()]
        # Should be possible to address by field and by number.
        assert pv_ptid["p"] == p_ptid
        assert pv_ptid["v"] == v_ptid
        assert pv_ptid[0] == p_ptid
        assert pv_ptid[1] == v_ptid
        # More complicated version.
        pv_t_ptid = self.pv_t_unit._physical_type_id
        t_ptid = self.t_unit._physical_type_id
        assert pv_t_ptid == np.array((pv_ptid, t_ptid), "O,O")[()]
        assert pv_t_ptid["pv"] == pv_ptid
        assert pv_t_ptid["t"] == t_ptid
        assert pv_t_ptid["pv"][1] == v_ptid

    def test_physical_type(self):
        pv_pt = self.pv_unit.physical_type
        assert pv_pt == np.array(("length", "speed"), "O,O")[()]

        pv_t_pt = self.pv_t_unit.physical_type
        assert pv_t_pt == np.array((pv_pt, "time"), "O,O")[()]

    def test_si(self):
        pv_t_si = self.pv_t_unit.si
        assert pv_t_si == self.pv_t_unit
        assert pv_t_si["pv"]["v"].scale == 1000

    def test_cgs(self):
        pv_t_cgs = self.pv_t_unit.cgs
        assert pv_t_cgs == self.pv_t_unit
        assert pv_t_cgs["pv"]["v"].scale == 100000

    def test_decompose(self):
        pv_t_decompose = self.pv_t_unit.decompose()
        assert pv_t_decompose["pv"]["v"].scale == 1000

    def test_is_equivalent(self):
        assert self.pv_unit.is_equivalent(("AU", "AU/day"))
        assert not self.pv_unit.is_equivalent("m")
        assert not self.pv_unit.is_equivalent(("AU", "AU"))
        # Names should be ignored.
        pv_alt = StructuredUnit("m,m/s", names=("q", "w"))
        assert pv_alt.field_names != self.pv_unit.field_names
        assert self.pv_unit.is_equivalent(pv_alt)
        # Regular units should work too.
        assert not u.m.is_equivalent(self.pv_unit)

    def test_conversion(self):
        pv1 = self.pv_unit.to(("AU", "AU/day"), self.pv)
        assert isinstance(pv1, np.ndarray)
        assert pv1.dtype == self.pv.dtype
        assert np.all(pv1["p"] * u.AU == self.pv["p"] * self.p_unit)
        assert np.all(pv1["v"] * u.AU / u.day == self.pv["v"] * self.v_unit)
        # Names should be from value.
        su2 = StructuredUnit((self.p_unit, self.v_unit), ("position", "velocity"))
        pv2 = su2.to(("Mm", "mm/s"), self.pv)
        assert pv2.dtype.names == ("p", "v")
        assert pv2.dtype == self.pv.dtype
        # Check recursion.
        pv_t1 = self.pv_t_unit.to((("AU", "AU/day"), "Myr"), self.pv_t)
        assert isinstance(pv_t1, np.ndarray)
        assert pv_t1.dtype == self.pv_t.dtype
        assert np.all(pv_t1["pv"]["p"] * u.AU == self.pv_t["pv"]["p"] * self.p_unit)
        assert np.all(
            pv_t1["pv"]["v"] * u.AU / u.day == self.pv_t["pv"]["v"] * self.v_unit
        )
        assert np.all(pv_t1["t"] * u.Myr == self.pv_t["t"] * self.t_unit)
        # Passing in tuples should work.
        pv_t2 = self.pv_t_unit.to((("AU", "AU/day"), "Myr"), ((1.0, 0.1), 10.0))
        assert pv_t2["pv"]["p"] == self.p_unit.to("AU", 1.0)
        assert pv_t2["pv"]["v"] == self.v_unit.to("AU/day", 0.1)
        assert pv_t2["t"] == self.t_unit.to("Myr", 10.0)
        pv_t3 = self.pv_t_unit.to(
            (("AU", "AU/day"), "Myr"), [((1.0, 0.1), 10.0), ((2.0, 0.2), 20.0)]
        )
        assert np.all(pv_t3["pv"]["p"] == self.p_unit.to("AU", [1.0, 2.0]))
        assert np.all(pv_t3["pv"]["v"] == self.v_unit.to("AU/day", [0.1, 0.2]))
        assert np.all(pv_t3["t"] == self.t_unit.to("Myr", [10.0, 20.0]))


class TestStructuredUnitArithmatic(StructuredTestBaseWithUnits):
    def test_multiplication(self):
        pv_times_au = self.pv_unit * u.au
        assert isinstance(pv_times_au, StructuredUnit)
        assert pv_times_au.field_names == ("p", "v")
        assert pv_times_au["p"] == self.p_unit * u.AU
        assert pv_times_au["v"] == self.v_unit * u.AU
        au_times_pv = u.au * self.pv_unit
        assert au_times_pv == pv_times_au
        pv_times_au2 = self.pv_unit * "au"
        assert pv_times_au2 == pv_times_au
        au_times_pv2 = "AU" * self.pv_unit
        assert au_times_pv2 == pv_times_au
        with pytest.raises(TypeError):
            self.pv_unit * self.pv_unit
        with pytest.raises(TypeError):
            "s,s" * self.pv_unit

    def test_division(self):
        pv_by_s = self.pv_unit / u.s
        assert isinstance(pv_by_s, StructuredUnit)
        assert pv_by_s.field_names == ("p", "v")
        assert pv_by_s["p"] == self.p_unit / u.s
        assert pv_by_s["v"] == self.v_unit / u.s
        pv_by_s2 = self.pv_unit / "s"
        assert pv_by_s2 == pv_by_s
        with pytest.raises(TypeError):
            1.0 / self.pv_unit
        with pytest.raises(TypeError):
            u.s / self.pv_unit


class TestStructuredQuantity(StructuredTestBaseWithUnits):
    def test_initialization_and_keying(self):
        q_pv = Quantity(self.pv, self.pv_unit)
        q_p = q_pv["p"]
        assert isinstance(q_p, Quantity)
        assert isinstance(q_p.unit, UnitBase)
        assert np.all(q_p == self.pv["p"] * self.pv_unit["p"])
        q_v = q_pv["v"]
        assert isinstance(q_v, Quantity)
        assert isinstance(q_v.unit, UnitBase)
        assert np.all(q_v == self.pv["v"] * self.pv_unit["v"])
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q_t = q_pv_t["t"]
        assert np.all(q_t == self.pv_t["t"] * self.pv_t_unit["t"])
        q_pv2 = q_pv_t["pv"]
        assert isinstance(q_pv2, Quantity)
        assert q_pv2.unit == self.pv_unit
        with pytest.raises(ValueError):
            Quantity(self.pv, self.pv_t_unit)
        with pytest.raises(ValueError):
            Quantity(self.pv_t, self.pv_unit)

    def test_initialization_with_unit_tuples(self):
        q_pv_t = Quantity(self.pv_t, (("km", "km/s"), "s"))
        assert isinstance(q_pv_t.unit, StructuredUnit)
        assert q_pv_t.unit == self.pv_t_unit

    def test_initialization_with_string(self):
        q_pv_t = Quantity(self.pv_t, "(km, km/s), s")
        assert isinstance(q_pv_t.unit, StructuredUnit)
        assert q_pv_t.unit == self.pv_t_unit

    def test_initialization_by_multiplication_with_unit(self):
        q_pv_t = self.pv_t * self.pv_t_unit
        assert q_pv_t.unit is self.pv_t_unit
        assert np.all(q_pv_t.value == self.pv_t)
        assert not np.may_share_memory(q_pv_t, self.pv_t)
        q_pv_t2 = self.pv_t_unit * self.pv_t
        assert q_pv_t.unit is self.pv_t_unit
        # Not testing equality of structured Quantity here.
        assert np.all(q_pv_t2.value == q_pv_t.value)

    def test_initialization_by_shifting_to_unit(self):
        q_pv_t = self.pv_t << self.pv_t_unit
        assert q_pv_t.unit is self.pv_t_unit
        assert np.all(q_pv_t.value == self.pv_t)
        assert np.may_share_memory(q_pv_t, self.pv_t)

    def test_initialization_without_unit(self):
        q_pv_t = u.Quantity(self.pv_t, unit=None)

        assert np.all(q_pv_t.value == self.pv_t)

        # Test that unit is a structured unit like the dtype
        expected_unit = _structured_unit_like_dtype(
            u.Quantity._default_unit, self.pv_t.dtype
        )
        assert q_pv_t.unit == expected_unit
        # A more explicit test
        assert q_pv_t.unit == u.StructuredUnit(((u.one, u.one), u.one))

    def test_getitem(self):
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q_pv_t01 = q_pv_t[:2]
        assert isinstance(q_pv_t01, Quantity)
        assert q_pv_t01.unit == q_pv_t.unit
        assert np.all(q_pv_t01["t"] == q_pv_t["t"][:2])
        q_pv_t1 = q_pv_t[1]
        assert isinstance(q_pv_t1, Quantity)
        assert q_pv_t1.unit == q_pv_t.unit
        assert q_pv_t1.shape == ()
        assert q_pv_t1["t"] == q_pv_t["t"][1]

    def test_value(self):
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        value = q_pv_t.value
        assert type(value) is np.ndarray
        assert np.all(value == self.pv_t)
        value1 = q_pv_t[1].value
        assert type(value1) is np.void
        assert np.all(value1 == self.pv_t[1])

    def test_conversion(self):
        q_pv = Quantity(self.pv, self.pv_unit)
        q1 = q_pv.to(("AU", "AU/day"))
        assert isinstance(q1, Quantity)
        assert q1["p"].unit == u.AU
        assert q1["v"].unit == u.AU / u.day
        assert np.all(q1["p"] == q_pv["p"].to(u.AU))
        assert np.all(q1["v"] == q_pv["v"].to(u.AU / u.day))
        q2 = q_pv.to(self.pv_unit)
        assert q2["p"].unit == self.p_unit
        assert q2["v"].unit == self.v_unit
        assert np.all(q2["p"].value == self.pv["p"])
        assert np.all(q2["v"].value == self.pv["v"])
        assert not np.may_share_memory(q2, q_pv)
        pv1 = q_pv.to_value(("AU", "AU/day"))
        assert type(pv1) is np.ndarray
        assert np.all(pv1["p"] == q_pv["p"].to_value(u.AU))
        assert np.all(pv1["v"] == q_pv["v"].to_value(u.AU / u.day))
        pv11 = q_pv[1].to_value(("AU", "AU/day"))
        assert type(pv11) is np.void
        assert pv11 == pv1[1]
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q2 = q_pv_t.to((("kpc", "kpc/Myr"), "Myr"))
        assert q2["pv"]["p"].unit == u.kpc
        assert q2["pv"]["v"].unit == u.kpc / u.Myr
        assert q2["t"].unit == u.Myr
        assert np.all(q2["pv"]["p"] == q_pv_t["pv"]["p"].to(u.kpc))
        assert np.all(q2["pv"]["v"] == q_pv_t["pv"]["v"].to(u.kpc / u.Myr))
        assert np.all(q2["t"] == q_pv_t["t"].to(u.Myr))

    def test_conversion_via_lshift(self):
        q_pv = Quantity(self.pv, self.pv_unit)
        q1 = q_pv << StructuredUnit(("AU", "AU/day"))
        assert isinstance(q1, Quantity)
        assert q1["p"].unit == u.AU
        assert q1["v"].unit == u.AU / u.day
        assert np.all(q1["p"] == q_pv["p"].to(u.AU))
        assert np.all(q1["v"] == q_pv["v"].to(u.AU / u.day))
        q2 = q_pv << self.pv_unit
        assert q2["p"].unit == self.p_unit
        assert q2["v"].unit == self.v_unit
        assert np.all(q2["p"].value == self.pv["p"])
        assert np.all(q2["v"].value == self.pv["v"])
        assert np.may_share_memory(q2, q_pv)
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q2 = q_pv_t << "(kpc,kpc/Myr),Myr"
        assert q2["pv"]["p"].unit == u.kpc
        assert q2["pv"]["v"].unit == u.kpc / u.Myr
        assert q2["t"].unit == u.Myr
        assert np.all(q2["pv"]["p"] == q_pv_t["pv"]["p"].to(u.kpc))
        assert np.all(q2["pv"]["v"] == q_pv_t["pv"]["v"].to(u.kpc / u.Myr))
        assert np.all(q2["t"] == q_pv_t["t"].to(u.Myr))

    def test_inplace_conversion(self):
        # In principle, in-place might be possible, in which case this should be
        # changed -- ie ``q1 is q_link``.
        q_pv = Quantity(self.pv, self.pv_unit)
        q1 = q_pv.copy()
        q_link = q1
        q1 <<= StructuredUnit(("AU", "AU/day"))
        assert q1 is not q_link
        assert q1["p"].unit == u.AU
        assert q1["v"].unit == u.AU / u.day
        assert np.all(q1["p"] == q_pv["p"].to(u.AU))
        assert np.all(q1["v"] == q_pv["v"].to(u.AU / u.day))
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q2 = q_pv_t.copy()
        q_link = q2
        q2 <<= "(kpc,kpc/Myr),Myr"
        assert q2 is not q_link
        assert q2["pv"]["p"].unit == u.kpc
        assert q2["pv"]["v"].unit == u.kpc / u.Myr
        assert q2["t"].unit == u.Myr
        assert np.all(q2["pv"]["p"] == q_pv_t["pv"]["p"].to(u.kpc))
        assert np.all(q2["pv"]["v"] == q_pv_t["pv"]["v"].to(u.kpc / u.Myr))
        assert np.all(q2["t"] == q_pv_t["t"].to(u.Myr))

    def test_si(self):
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q_pv_t_si = q_pv_t.si
        assert_array_equal(q_pv_t_si, q_pv_t.to("(m,m/s),s"))

    def test_cgs(self):
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q_pv_t_cgs = q_pv_t.cgs
        assert_array_equal(q_pv_t_cgs, q_pv_t.to("(cm,cm/s),s"))

    def test_equality(self):
        q_pv = Quantity(self.pv, self.pv_unit)
        equal = q_pv == q_pv
        not_equal = q_pv != q_pv
        assert np.all(equal)
        assert not np.any(not_equal)
        equal2 = q_pv == q_pv[1]
        not_equal2 = q_pv != q_pv[1]
        assert np.all(equal2 == [False, True, False])
        assert np.all(not_equal2 != equal2)
        q1 = q_pv.to(("AU", "AU/day"))
        # Ensure same conversion is done, by placing q1 first.
        assert np.all(q1 == q_pv)
        assert not np.any(q1 != q_pv)
        # Check different names in dtype.
        assert np.all(q1.value * u.Unit("AU, AU/day") == q_pv)
        assert not np.any(q1.value * u.Unit("AU, AU/day") != q_pv)
        assert (q_pv == "b") is False
        assert ("b" != q_pv) is True
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        assert np.all((q_pv_t[2] == q_pv_t) == [False, False, True])
        assert np.all((q_pv_t[2] != q_pv_t) != [False, False, True])
        assert (q_pv == q_pv_t) is False
        assert (q_pv_t != q_pv) is True

    def test_setitem(self):
        q_pv = Quantity(self.pv, self.pv_unit)
        q_pv[1] = (2.0, 2.0) * self.pv_unit
        assert q_pv[1].value == np.array((2.0, 2.0), self.pv_dtype)
        q_pv[1:2] = (1.0, 0.5) * u.Unit("AU, AU/day")
        assert q_pv["p"][1] == 1.0 * u.AU
        assert q_pv["v"][1] == 0.5 * u.AU / u.day
        q_pv["v"] = 1.0 * u.km / u.s
        assert np.all(q_pv["v"] == 1.0 * u.km / u.s)
        with pytest.raises(u.UnitsError):
            q_pv[1] = (1.0, 1.0) * u.Unit("AU, AU")
        with pytest.raises(u.UnitsError):
            q_pv["v"] = 1.0 * u.km
        q_pv_t = Quantity(self.pv_t, self.pv_t_unit)
        q_pv_t[1] = ((2.0, 2.0), 3.0) * self.pv_t_unit
        assert q_pv_t[1].value == np.array(((2.0, 2.0), 3.0), self.pv_t_dtype)
        q_pv_t[1:2] = ((1.0, 0.5), 5.0) * u.Unit("(AU, AU/day), yr")
        assert q_pv_t["pv"][1] == (1.0, 0.5) * u.Unit("AU, AU/day")
        assert q_pv_t["t"][1] == 5.0 * u.yr
        q_pv_t["pv"] = (1.0, 0.5) * self.pv_unit
        assert np.all(q_pv_t["pv"] == (1.0, 0.5) * self.pv_unit)


class TestStructuredQuantityFunctions(StructuredTestBaseWithUnits):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        cls.q_pv = cls.pv << cls.pv_unit
        cls.q_pv_t = cls.pv_t << cls.pv_t_unit

    def test_empty_like(self):
        z = np.empty_like(self.q_pv)
        assert z.dtype == self.pv_dtype
        assert z.unit == self.pv_unit
        assert z.shape == self.pv.shape

    @pytest.mark.parametrize("func", [np.zeros_like, np.ones_like])
    def test_zeros_ones_like(self, func):
        z = func(self.q_pv)
        assert z.dtype == self.pv_dtype
        assert z.unit == self.pv_unit
        assert z.shape == self.pv.shape
        assert_array_equal(z, func(self.pv) << self.pv_unit)

    def test_structured_to_unstructured(self):
        # can't unstructure something with incompatible units
        with pytest.raises(u.UnitConversionError, match="'km / s'"):
            rfn.structured_to_unstructured(self.q_pv)

        # For the other tests of ``structured_to_unstructured``, see
        # ``test_quantity_non_ufuncs.TestRecFunctions.test_structured_to_unstructured``

    def test_unstructured_to_structured(self):
        # can't structure something that's already structured
        dtype = np.dtype([("f1", float), ("f2", float)])
        with pytest.raises(ValueError, match="The length of the last dimension"):
            rfn.unstructured_to_structured(self.q_pv, dtype=self.q_pv.dtype)

        # For the other tests of ``structured_to_unstructured``, see
        # ``test_quantity_non_ufuncs.TestRecFunctions.test_unstructured_to_structured``


class TestStructuredSpecificTypeQuantity(StructuredTestBaseWithUnits):
    def setup_class(self):
        super().setup_class()

        class PositionVelocity(u.SpecificTypeQuantity):
            _equivalent_unit = self.pv_unit

        self.PositionVelocity = PositionVelocity

    def test_init(self):
        pv = self.PositionVelocity(self.pv, self.pv_unit)
        assert isinstance(pv, self.PositionVelocity)
        assert type(pv["p"]) is u.Quantity
        assert_array_equal(pv["p"], self.pv["p"] << self.pv_unit["p"])

        pv2 = self.PositionVelocity(self.pv, "AU,AU/day")
        assert_array_equal(pv2["p"], self.pv["p"] << u.AU)

    def test_error_on_non_equivalent_unit(self):
        with pytest.raises(u.UnitsError):
            self.PositionVelocity(self.pv, "AU")
        with pytest.raises(u.UnitsError):
            self.PositionVelocity(self.pv, "AU,yr")


class TestStructuredLogUnit:
    def setup_class(self):
        self.mag_time_dtype = np.dtype([("mag", "f8"), ("t", "f8")])
        self.mag_time = np.array([(20.0, 10.0), (25.0, 100.0)], self.mag_time_dtype)

    def test_unit_initialization(self):
        mag_time_unit = StructuredUnit((u.STmag, u.s), self.mag_time_dtype)
        assert mag_time_unit["mag"] == u.STmag
        assert mag_time_unit["t"] == u.s

        mag_time_unit2 = u.Unit("mag(ST),s")
        assert mag_time_unit2 == mag_time_unit

    def test_quantity_initialization(self):
        su = u.Unit("mag(ST),s")
        mag_time = self.mag_time << su
        assert isinstance(mag_time["mag"], u.Magnitude)
        assert isinstance(mag_time["t"], u.Quantity)
        assert mag_time.unit == su
        assert_array_equal(mag_time["mag"], self.mag_time["mag"] << u.STmag)
        assert_array_equal(mag_time["t"], self.mag_time["t"] << u.s)

    def test_quantity_si(self):
        mag_time = self.mag_time << u.Unit("mag(ST),yr")
        mag_time_si = mag_time.si
        assert_array_equal(mag_time_si["mag"], mag_time["mag"].si)
        assert_array_equal(mag_time_si["t"], mag_time["t"].si)


class TestStructuredMaskedQuantity(StructuredTestBaseWithUnits):
    """Somewhat minimal tests.  Conversion is most stringent."""

    def setup_class(self):
        super().setup_class()
        self.qpv = self.pv << self.pv_unit
        self.pv_mask = np.array(
            [
                (True, False),
                (False, False),
                (False, True),
            ],
            [("p", bool), ("v", bool)],
        )
        self.mpv = Masked(self.qpv, mask=self.pv_mask)

    def test_init(self):
        assert isinstance(self.mpv, Masked)
        assert isinstance(self.mpv, Quantity)
        assert_array_equal(self.mpv.unmasked, self.qpv)
        assert_array_equal(self.mpv.mask, self.pv_mask)

    def test_slicing(self):
        mp = self.mpv["p"]
        assert isinstance(mp, Masked)
        assert isinstance(mp, Quantity)
        assert_array_equal(mp.unmasked, self.qpv["p"])
        assert_array_equal(mp.mask, self.pv_mask["p"])

    def test_conversion(self):
        mpv = self.mpv.to("AU,AU/day")
        assert isinstance(mpv, Masked)
        assert isinstance(mpv, Quantity)
        assert_array_equal(mpv.unmasked, self.qpv.to("AU,AU/day"))
        assert_array_equal(mpv.mask, self.pv_mask)
        assert np.all(mpv == self.mpv)

    def test_si(self):
        mpv = self.mpv.si
        assert isinstance(mpv, Masked)
        assert isinstance(mpv, Quantity)
        assert_array_equal(mpv.unmasked, self.qpv.si)
        assert_array_equal(mpv.mask, self.pv_mask)
        assert np.all(mpv == self.mpv)
