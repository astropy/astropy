# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test Structured units and quantities specifically with the ERFA ufuncs.
"""

import numpy as np
import pytest
from erfa import ufunc as erfa_ufunc
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.units.quantity_helper.erfa import astrom_unit

ASTROM_UNIT = astrom_unit()


def vvd(val, valok, dval, func, test, status):
    """Mimic routine of erfa/src/t_erfa_c.c (to help copy & paste)"""
    assert u.allclose(val, valok * val.unit, atol=dval * val.unit)


class TestPVUfuncs:
    @classmethod
    def setup_class(cls):
        cls.pv_unit = u.Unit("AU,AU/day")
        cls.pv_value = np.array(
            [
                ([1.0, 0.0, 0.0], [0.0, 0.0125, 0.0]),
                ([0.0, 1.0, 0.0], [-0.0125, 0.0, 0.0]),
            ],
            dtype=erfa_ufunc.dt_pv,
        )
        cls.pv = cls.pv_value << cls.pv_unit

    def test_cpv(self):
        pv_copy = erfa_ufunc.cpv(self.pv)
        assert_array_equal(pv_copy, self.pv)
        assert not np.may_share_memory(pv_copy, self.pv)

    def test_p2pv(self):
        p2pv = erfa_ufunc.p2pv(self.pv["p"])
        assert_array_equal(p2pv["p"], self.pv["p"])
        assert_array_equal(
            p2pv["v"], np.zeros(self.pv.shape + (3,), float) << u.m / u.s
        )

    def test_p2pv_inplace(self):
        # TODO: fix np.zeros_like.
        out = np.zeros_like(self.pv_value) << self.pv_unit
        p2pv = erfa_ufunc.p2pv(self.pv["p"], out=out)
        assert out is p2pv
        assert_array_equal(p2pv["p"], self.pv["p"])
        assert_array_equal(
            p2pv["v"], np.zeros(self.pv.shape + (3,), float) << u.m / u.s
        )

    def test_pv2p(self):
        p = erfa_ufunc.pv2p(self.pv)
        assert_array_equal(p, self.pv["p"])
        out = np.zeros_like(p)
        p2 = erfa_ufunc.pv2p(self.pv, out=out)
        assert out is p2
        assert_array_equal(p2, self.pv["p"])

    def test_pv2s(self):
        theta, phi, r, td, pd, rd = erfa_ufunc.pv2s(self.pv)
        assert theta.unit == u.radian
        assert_quantity_allclose(theta, [0, 90] * u.deg)  # longitude
        assert phi.unit == u.radian
        assert_array_equal(phi.value, np.zeros(self.pv.shape))  # latitude
        assert r.unit == u.AU
        assert_array_equal(r.value, np.ones(self.pv.shape))
        assert td.unit == u.radian / u.day
        assert_array_equal(td.value, np.array([0.0125] * 2))
        assert pd.unit == u.radian / u.day
        assert_array_equal(pd.value, np.zeros(self.pv.shape))
        assert rd.unit == u.AU / u.day
        assert_array_equal(rd.value, np.zeros(self.pv.shape))

    def test_pv2s_non_standard_units(self):
        pv = self.pv_value << u.Unit("Pa,Pa/m")
        theta, phi, r, td, pd, rd = erfa_ufunc.pv2s(pv)
        assert theta.unit == u.radian
        assert_quantity_allclose(theta, [0, 90] * u.deg)  # longitude
        assert phi.unit == u.radian
        assert_array_equal(phi.value, np.zeros(pv.shape))  # latitude
        assert r.unit == u.Pa
        assert_array_equal(r.value, np.ones(pv.shape))
        assert td.unit == u.radian / u.m
        assert_array_equal(td.value, np.array([0.0125] * 2))
        assert pd.unit == u.radian / u.m
        assert_array_equal(pd.value, np.zeros(pv.shape))
        assert rd.unit == u.Pa / u.m
        assert_array_equal(rd.value, np.zeros(pv.shape))

    @pytest.mark.xfail(
        reason=(
            "erfa ufuncs cannot take different names; it is not yet clear whether "
            "this is changeable; see https://github.com/liberfa/pyerfa/issues/77"
        )
    )
    def test_pv2s_non_standard_names_and_units(self):
        pv_value = np.array(self.pv_value, dtype=[("pos", "f8"), ("vel", "f8")])
        pv = pv_value << u.Unit("Pa,Pa/m")
        theta, phi, r, td, pd, rd = erfa_ufunc.pv2s(pv)
        assert theta.unit == u.radian
        assert_quantity_allclose(theta, [0, 90] * u.deg)  # longitude
        assert phi.unit == u.radian
        assert_array_equal(phi.value, np.zeros(pv.shape))  # latitude
        assert r.unit == u.Pa
        assert_array_equal(r.value, np.ones(pv.shape))
        assert td.unit == u.radian / u.m
        assert_array_equal(td.value, np.array([0.0125] * 2))
        assert pd.unit == u.radian / u.m
        assert_array_equal(pd.value, np.zeros(pv.shape))
        assert rd.unit == u.Pa / u.m
        assert_array_equal(rd.value, np.zeros(pv.shape))

    def test_s2pv(self):
        theta, phi, r, td, pd, rd = erfa_ufunc.pv2s(self.pv)
        # On purpose change some of the units away from expected by s2pv.
        pv = erfa_ufunc.s2pv(
            theta.to(u.deg), phi, r.to(u.m), td.to(u.deg / u.day), pd, rd.to(u.m / u.s)
        )
        assert pv.unit == u.StructuredUnit("m, m/s", names=("p", "v"))
        assert_quantity_allclose(pv["p"], self.pv["p"], atol=1 * u.m, rtol=0)
        assert_quantity_allclose(pv["v"], self.pv["v"], atol=1 * u.mm / u.s, rtol=0)

    def test_s2p_not_all_quantity(self):
        # Test for a useful error message - see gh-16873.
        # Non-quantity input should be treated as dimensionless and thus cannot
        # be converted to radians.
        with pytest.raises(
            AttributeError,
            match=(
                "'NoneType' object has no attribute 'get_converter'"
                ".*\n.*treated as dimensionless"
            ),
        ):
            erfa_ufunc.s2p(0.5, 0.5, 4 * u.km)

        # Except if we have the right equivalency in place.
        with u.add_enabled_equivalencies(u.dimensionless_angles()):
            result = erfa_ufunc.s2p(0.5, 0.5, 4 * u.km)

        expected = erfa_ufunc.s2p(0.5 * u.radian, 0.5 * u.radian, 4 * u.km)
        assert_array_equal(result, expected)

    def test_pvstar(self):
        ra, dec, pmr, pmd, px, rv, stat = erfa_ufunc.pvstar(self.pv)
        assert_array_equal(stat, np.zeros(self.pv.shape, dtype="i4"))
        assert ra.unit == u.radian
        assert_quantity_allclose(ra, [0, 90] * u.deg)
        assert dec.unit == u.radian
        assert_array_equal(dec.value, np.zeros(self.pv.shape))  # latitude
        assert pmr.unit == u.radian / u.year
        assert_quantity_allclose(pmr, [0.0125, 0.0125] * u.radian / u.day)
        assert pmd.unit == u.radian / u.year
        assert_array_equal(pmd.value, np.zeros(self.pv.shape))
        assert px.unit == u.arcsec
        assert_quantity_allclose(px, 1 * u.radian)
        assert rv.unit == u.km / u.s
        # RV is non-zero because proper motion induces a small redshift
        # due to second order Doppler shift.
        assert_quantity_allclose(
            rv, np.zeros(self.pv.shape) << (u.km / u.s), atol=1 * u.m / u.s
        )

    def test_starpv(self):
        ra, dec, pmr, pmd, px, rv, stat = erfa_ufunc.pvstar(self.pv)
        pv, stat = erfa_ufunc.starpv(
            ra.to(u.deg), dec.to(u.deg), pmr, pmd, px, rv.to(u.m / u.s)
        )
        assert_array_equal(stat, np.zeros(self.pv.shape, dtype="i4"))
        assert pv.unit == self.pv.unit
        # Roundtrip is not as good as hoped on 32bit, not clear why.
        # But proper motions are ridiculously high...
        assert_quantity_allclose(pv["p"], self.pv["p"], atol=1 * u.m, rtol=0)
        assert_quantity_allclose(pv["v"], self.pv["v"], atol=1 * u.m / u.s, rtol=0)

    def test_pvtob(self):
        pv = erfa_ufunc.pvtob(
            [90, 0] * u.deg,
            0.0 * u.deg,
            100 * u.km,
            0 * u.deg,
            0 * u.deg,
            0 * u.deg,
            90 * u.deg,
        )
        assert pv.unit == u.StructuredUnit("m, m/s", names=("p", "v"))
        assert pv.unit["v"] == u.m / u.s
        assert_quantity_allclose(
            pv["p"], [[-6478, 0, 0], [0, 6478, 0]] * u.km, atol=2 * u.km
        )
        assert_quantity_allclose(
            pv["v"], [[0, -0.5, 0], [-0.5, 0, 0]] * u.km / u.s, atol=0.1 * u.km / u.s
        )

    def test_pvdpv(self):
        pvdpv = erfa_ufunc.pvdpv(self.pv, self.pv)
        assert pvdpv["pdp"].unit == self.pv.unit["p"] ** 2
        assert pvdpv["pdv"].unit == self.pv.unit["p"] * self.pv.unit["v"]
        assert_array_equal(
            pvdpv["pdp"], np.einsum("...i,...i->...", self.pv["p"], self.pv["p"])
        )
        assert_array_equal(
            pvdpv["pdv"], 2 * np.einsum("...i,...i->...", self.pv["p"], self.pv["v"])
        )
        z_axis = u.Quantity(np.array(([0, 0, 1], [0, 0, 0]), erfa_ufunc.dt_pv), "1,1/s")
        pvdpv2 = erfa_ufunc.pvdpv(self.pv, z_axis)
        assert pvdpv2["pdp"].unit == self.pv.unit["p"]
        assert pvdpv2["pdv"].unit == self.pv.unit["v"]
        assert_array_equal(pvdpv2["pdp"].value, np.zeros(self.pv.shape))
        assert_array_equal(pvdpv2["pdv"].value, np.zeros(self.pv.shape))

    def test_pvxpv(self):
        pvxpv = erfa_ufunc.pvxpv(self.pv, self.pv)
        assert pvxpv["p"].unit == self.pv.unit["p"] ** 2
        assert pvxpv["v"].unit == self.pv.unit["p"] * self.pv.unit["v"]
        assert_array_equal(pvxpv["p"].value, np.zeros(self.pv["p"].shape))
        assert_array_equal(pvxpv["v"].value, np.zeros(self.pv["v"].shape))
        z_axis = u.Quantity(np.array(([0, 0, 1], [0, 0, 0]), erfa_ufunc.dt_pv), "1,1/s")
        pvxpv2 = erfa_ufunc.pvxpv(self.pv, z_axis)
        assert pvxpv2["p"].unit == self.pv.unit["p"]
        assert pvxpv2["v"].unit == self.pv.unit["v"]
        assert_array_equal(pvxpv2["p"], [[0.0, -1, 0.0], [1.0, 0.0, 0.0]] * u.AU)
        assert_array_equal(
            pvxpv2["v"], [[0.0125, 0.0, 0.0], [0.0, 0.0125, 0.0]] * u.AU / u.day
        )

    def test_pvm(self):
        pm, vm = erfa_ufunc.pvm(self.pv)
        assert pm.unit == self.pv.unit["p"]
        assert vm.unit == self.pv.unit["v"]
        assert_array_equal(pm, np.linalg.norm(self.pv["p"], axis=-1))
        assert_array_equal(vm, np.linalg.norm(self.pv["v"], axis=-1))

    def test_pvmpv(self):
        pvmpv = erfa_ufunc.pvmpv(self.pv, self.pv)
        assert pvmpv.unit == self.pv.unit
        assert_array_equal(pvmpv["p"], 0 * self.pv["p"])
        assert_array_equal(pvmpv["v"], 0 * self.pv["v"])

    def test_pvppv(self):
        pvppv = erfa_ufunc.pvppv(self.pv, self.pv)
        assert pvppv.unit == self.pv.unit
        assert_array_equal(pvppv["p"], 2 * self.pv["p"])
        assert_array_equal(pvppv["v"], 2 * self.pv["v"])

    def test_pvu(self):
        pvu = erfa_ufunc.pvu(86400 * u.s, self.pv)
        assert pvu.unit == self.pv.unit
        assert_array_equal(pvu["p"], self.pv["p"] + 1 * u.day * self.pv["v"])
        assert_array_equal(pvu["v"], self.pv["v"])

    def test_pvup(self):
        pvup = erfa_ufunc.pvup(86400 * u.s, self.pv)
        assert pvup.unit == self.pv.unit["p"]
        assert_array_equal(pvup, self.pv["p"] + 1 * u.day * self.pv["v"])

    def test_sxpv(self):
        # Not a realistic example!!
        sxpv = erfa_ufunc.sxpv(10.0, self.pv)
        assert sxpv.unit == self.pv.unit
        assert_array_equal(sxpv["p"], self.pv["p"] * 10)
        assert_array_equal(sxpv["v"], self.pv["v"] * 10)
        sxpv2 = erfa_ufunc.sxpv(30.0 * u.s, self.pv)
        assert sxpv2.unit == u.StructuredUnit("AU s,AU s/d", names=("p", "v"))
        assert_array_equal(sxpv2["p"], self.pv["p"] * 30 * u.s)
        assert_array_equal(sxpv2["v"], self.pv["v"] * 30 * u.s)

    def test_s2xpv(self):
        # Not a realistic example!!
        s2xpv = erfa_ufunc.s2xpv(10.0, 1 * u.s, self.pv)
        assert s2xpv.unit == u.StructuredUnit("AU,AU s/d", names=("p", "v"))
        assert_array_equal(s2xpv["p"], self.pv["p"] * 10)
        assert_array_equal(s2xpv["v"], self.pv["v"] * u.s)

    @pytest.mark.parametrize(
        "r",
        [
            np.eye(3),
            np.array(
                [
                    [0.0, -1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.eye(3) / u.s,
        ],
    )
    def test_rxpv(self, r):
        result = erfa_ufunc.rxpv(r, self.pv)
        assert_array_equal(result["p"], np.einsum("...ij,...j->...i", r, self.pv["p"]))
        assert_array_equal(result["v"], np.einsum("...ij,...j->...i", r, self.pv["v"]))

    @pytest.mark.parametrize(
        "r",
        [
            np.eye(3),
            np.array(
                [
                    [0.0, -1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.eye(3) / u.s,
        ],
    )
    def test_trxpv(self, r):
        result = erfa_ufunc.trxpv(r, self.pv)
        assert_array_equal(
            result["p"], np.einsum("...ij,...j->...i", r.T, self.pv["p"])
        )
        assert_array_equal(
            result["v"], np.einsum("...ij,...j->...i", r.T, self.pv["v"])
        )


class TestEraStructUfuncs:
    @classmethod
    def setup_class(cls):
        ldbody = np.array(
            [
                (0.00028574, 3e-10, ([-7.81014427, -5.60956681, -1.98079819],
                                     [0.0030723249, -0.00406995477, -0.00181335842])),
                (0.00095435, 3e-9, ([0.738098796, 4.63658692, 1.9693136],
                                    [-0.00755816922, 0.00126913722, 0.000727999001])),
                (1.0, 6e-6, ([-0.000712174377, -0.00230478303, -0.00105865966],
                             [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]))
            ],
            dtype=erfa_ufunc.dt_eraLDBODY
        )  # fmt: skip
        ldbody_unit = u.StructuredUnit("Msun,radian,(AU,AU/day)", ldbody.dtype)
        cls.ldbody = ldbody << ldbody_unit
        cls.ob = [-0.974170437, -0.2115201, -0.0917583114] << u.AU
        cls.sc = np.array([-0.763276255, -0.608633767, -0.216735543])

        # From t_atciq in t_erfa_c.c
        astrom, eo = erfa_ufunc.apci13(2456165.5, 0.401182685)
        cls.astrom = astrom << ASTROM_UNIT
        cls.rc = 2.71 * u.rad
        cls.dc = 0.174 * u.rad
        cls.pr = 1e-5 * u.rad / u.year
        cls.pd = 5e-6 * u.rad / u.year
        cls.px = 0.1 * u.arcsec
        cls.rv = 55.0 * u.km / u.s

    def test_ldn_basic(self):
        sn = erfa_ufunc.ldn(self.ldbody, self.ob, self.sc)
        assert_quantity_allclose(
            sn,
            [-0.7632762579693333866, -0.6086337636093002660, -0.2167355420646328159]
            * u.one,
            atol=1e-12,
            rtol=0,
        )

    def test_ldn_in_other_unit(self):
        ldbody = self.ldbody.to("kg,rad,(m,m/s)")
        ob = self.ob.to("m")
        sn = erfa_ufunc.ldn(ldbody, ob, self.sc)
        assert_quantity_allclose(
            sn,
            [-0.7632762579693333866, -0.6086337636093002660, -0.2167355420646328159]
            * u.one,
            atol=1e-12,
            rtol=0,
        )

    def test_ldn_in_SI(self):
        sn = erfa_ufunc.ldn(self.ldbody.si, self.ob.si, self.sc)
        assert_quantity_allclose(
            sn,
            [-0.7632762579693333866, -0.6086337636093002660, -0.2167355420646328159]
            * u.one,
            atol=1e-12,
            rtol=0,
        )

    def test_aper(self):
        along = self.astrom["along"]
        astrom2 = erfa_ufunc.aper(10 * u.deg, self.astrom)
        assert astrom2["eral"].unit == u.radian
        assert_quantity_allclose(astrom2["eral"], along + 10 * u.deg)
        astrom3 = self.astrom.to("s,km,1,km,1,1,1,deg,deg,deg,deg,1,1,1,rad,rad,rad")
        astrom4 = erfa_ufunc.aper(10 * u.deg, astrom3)
        assert astrom3["eral"].unit == u.rad
        assert astrom4["eral"].unit == u.deg
        assert astrom4.unit == "s,km,1,km,1,1,1,deg,deg,deg,deg,1,1,1,deg,rad,rad"
        assert_quantity_allclose(astrom4["eral"], along + 10 * u.deg)

    def test_atciq_basic(self):
        ri, di = erfa_ufunc.atciq(
            self.rc, self.dc, self.pr, self.pd, self.px, self.rv, self.astrom
        )
        assert_quantity_allclose(ri, 2.710121572968696744 * u.rad)
        assert_quantity_allclose(di, 0.1729371367219539137 * u.rad)

    def test_atciq_in_other_unit(self):
        astrom = self.astrom.to("s,km,1,km,1,1,1,deg,deg,deg,deg,1,1,1,deg,deg,deg")
        ri, di = erfa_ufunc.atciq(
            self.rc.to(u.deg),
            self.dc.to(u.deg),
            self.pr.to(u.mas / u.yr),
            self.pd.to(u.mas / u.yr),
            self.px,
            self.rv.to(u.m / u.s),
            astrom,
        )
        assert_quantity_allclose(ri, 2.710121572968696744 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(di, 0.1729371367219539137 * u.rad, atol=1e-12 * u.rad)

    def test_atciqn(self):
        ri, di = erfa_ufunc.atciqn(
            self.rc.to(u.deg),
            self.dc.to(u.deg),
            self.pr.to(u.mas / u.yr),
            self.pd.to(u.mas / u.yr),
            self.px,
            self.rv.to(u.m / u.s),
            self.astrom.si,
            self.ldbody.si,
        )
        assert_quantity_allclose(ri, 2.710122008104983335 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(di, 0.1729371916492767821 * u.rad, atol=1e-12 * u.rad)

    def test_atciqz(self):
        ri, di = erfa_ufunc.atciqz(self.rc.to(u.deg), self.dc.to(u.deg), self.astrom.si)
        assert_quantity_allclose(ri, 2.709994899247256984 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(di, 0.1728740720984931891 * u.rad, atol=1e-12 * u.rad)

    def test_aticq(self):
        ri = 2.710121572969038991 * u.rad
        di = 0.1729371367218230438 * u.rad
        rc, dc = erfa_ufunc.aticq(ri.to(u.deg), di.to(u.deg), self.astrom.si)
        assert_quantity_allclose(rc, 2.710126504531716819 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(dc, 0.1740632537627034482 * u.rad, atol=1e-12 * u.rad)

    def test_aticqn(self):
        ri = 2.709994899247599271 * u.rad
        di = 0.1728740720983623469 * u.rad
        rc, dc = erfa_ufunc.aticqn(
            ri.to(u.deg), di.to(u.deg), self.astrom.si, self.ldbody.si
        )
        assert_quantity_allclose(rc, 2.709999575033027333 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(dc, 0.1739999656316469990 * u.rad, atol=1e-12 * u.rad)

    def test_atioq_atoiq(self):
        astrom, _ = erfa_ufunc.apio13(
            2456384.5,
            0.969254051,
            0.1550675,
            -0.527800806,
            -1.2345856,
            2738.0,
            2.47230737e-7,
            1.82640464e-6,
            731.0,
            12.8,
            0.59,
            0.55,
        )
        astrom = astrom << ASTROM_UNIT

        ri = 2.710121572969038991 * u.rad
        di = 0.1729371367218230438 * u.rad
        aob, zob, hob, dob, rob = erfa_ufunc.atioq(
            ri.to(u.deg), di.to(u.deg), astrom.si
        )
        assert_quantity_allclose(
            aob, 0.9233952224895122499e-1 * u.rad, atol=1e-12 * u.rad
        )
        assert_quantity_allclose(zob, 1.407758704513549991 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(
            hob, -0.9247619879881698140e-1 * u.rad, atol=1e-12 * u.rad
        )
        assert_quantity_allclose(dob, 0.1717653435756234676 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(rob, 2.710085107988480746 * u.rad, atol=1e-12 * u.rad)

        # Sadly does not just use the values from above.
        ob1 = 2.710085107986886201 * u.rad
        ob2 = 0.1717653435758265198 * u.rad
        ri2, di2 = erfa_ufunc.atoiq("R", ob1.to(u.deg), ob2.to(u.deg), astrom.si)
        assert_quantity_allclose(ri2, 2.710121574447540810 * u.rad, atol=1e-12 * u.rad)
        assert_quantity_allclose(
            di2, 0.17293718391166087785 * u.rad, atol=1e-12 * u.rad
        )


class TestAp:
    @classmethod
    def setup_class(cls):
        # Values from apco13 and apio test cases in t_erfa_c.c
        # with units changed so we can check conversion.
        cls.utc1 = 2456384.5
        cls.utc2 = 0.969254051
        cls.dut1 = (0.1550675 * u.s).to(u.ms)
        cls.sp = (-3.01974337e-11 * u.rad).to(u.deg)
        cls.theta = (3.14540971 * u.rad).to(u.mdeg)
        cls.elong = (-0.527800806 * u.rad).to(u.Marcsec)
        cls.phi = (-1.2345856 * u.rad).to(u.arcmin)
        cls.hm = (2738.0 * u.m).to(u.km)
        cls.phpa = (731.0 * u.hPa).to(u.Pa)
        cls.tc = (12.8 * u.deg_C).to(u.K, equivalencies=u.temperature())
        cls.rh = (0.59 * u.one).to(u.percent)
        cls.wl = (0.55 * u.micron).to(u.AA)
        # For apio.
        cls.xp = (2.47230737e-7 * u.rad).to(u.arcsec)
        cls.yp = (1.82640464e-6 * u.rad).to(u.arcmin)
        cls.refa = (0.000201418779 * u.rad).to(u.mas)
        cls.refb = (-2.36140831e-7 * u.rad).to(u.uas)
        cls.apco13_args = [
            getattr(cls, name)
            for name in [
                "utc1",
                "utc2",
                "dut1",
                "elong",
                "phi",
                "hm",
                "xp",
                "yp",
                "phpa",
                "tc",
                "rh",
                "wl",
            ]
        ]
        cls.apio_args = [
            getattr(cls, name)
            for name in [
                "sp",
                "theta",
                "elong",
                "phi",
                "hm",
                "xp",
                "yp",
                "refa",
                "refb",
            ]
        ]

    def test_apco13(self):
        astrom, eo, status = erfa_ufunc.apco13(*self.apco13_args)
        assert status == 0
        vvd(eo, -0.003020548354802412839, 1e-14, "eraApco13", "eo", status)
        assert astrom.unit == ASTROM_UNIT
        for name, expected in {
            "pmt": 13.25248468622475727,
            "eb": [-0.9741827107320875162,
                   -0.2115130190489716682,
                   -0.09179840189496755339],
            "eh": [-0.9736425572586935247,
                   -0.2092452121603336166,
                   -0.09075578153885665295],
            "em": 0.9998233240913898141,
            "v": [0.2078704994520489246e-4,
                  -0.8955360133238868938e-4,
                  -0.3863338993055887398e-4],
            "bm1": 0.9999999950277561004,
            "bpn": np.array(
                [[ 0.9999991390295147999,
                   0.4978650075315529277e-7,
                   0.001312227200850293372],
                 [-0.1136336652812486604e-7,
                   0.9999999995713154865,
                  -0.2928086230975367296e-4],
                 [-0.001312227201745553566,
                   0.2928082218847679162e-4,
                   0.9999991386008312212],
                ]).T,
            "along": -0.5278008060295995733,
            "xpl": 0.1133427418130752958e-5,
            "ypl": 0.1453347595780646207e-5,
            "sphi": -0.9440115679003211329,
            "cphi": 0.3299123514971474711,
            "diurab": 0,
            "eral": 2.617608909189664000,
            "refa": 0.2014187785940396921e-3,
            "refb": -0.2361408314943696227e-6,
        }.items():  # fmt: skip
            assert_quantity_allclose(astrom[name], expected * ASTROM_UNIT[name])

    def test_apco13_no_time_units(self):
        msg = "cannot pass in units for 2-part time in apco13"
        with pytest.raises(TypeError, match=msg):
            erfa_ufunc.apco13(self.utc1 * u.day, *self.apco13_args[1:])

        with pytest.raises(TypeError, match=msg):
            erfa_ufunc.apco13(self.utc1, self.utc2 * u.day, *self.apco13_args[2:])

        with pytest.raises(TypeError, match=msg):
            erfa_ufunc.apco13(
                self.utc1 * u.day, self.utc2 * u.day, *self.apco13_args[2:]
            )

    def test_apio(self):
        astrom = erfa_ufunc.apio(*self.apio_args)
        assert astrom.unit == ASTROM_UNIT
        for name, expected in {
            "along": -0.5278008060295995734,
            "xpl": 0.1133427418130752958e-5,
            "ypl": 0.1453347595780646207e-5,
            "sphi": -0.9440115679003211329,
            "cphi": 0.3299123514971474711,
            "diurab": 0.5135843661699913529e-6,
            "eral": 2.617608903970400427,
            "refa": 0.2014187790000000000e-3,
            "refb": -0.2361408310000000000e-6,
        }.items():
            assert_quantity_allclose(astrom[name], expected * ASTROM_UNIT[name])


class TestGeodetic:
    @classmethod
    def setup_class(cls):
        cls.ellipsoid = 1
        cls.length_unit = u.Unit("m")
        cls.equatorial_radius_value = 6378136.0
        cls.equatorial_radius = (cls.equatorial_radius_value << cls.length_unit).to(
            u.km
        )
        cls.flattening = 0.0033528 * u.dimensionless_unscaled
        cls.lon_value = 0.9827937232473290680
        cls.lon_unit = u.Unit("rad")
        cls.lon = (cls.lon_value << cls.lon_unit).to(u.deg)
        cls.lat_value = 0.9716018377570411532
        cls.lat_unit = u.Unit("rad")
        cls.lat = (cls.lat_value << cls.lat_unit).to(u.deg)
        cls.height_value = 332.36862495764397
        cls.height = cls.height_value << cls.length_unit
        cls.x_value = 2e6
        cls.x = (cls.x_value << cls.length_unit).to(u.cm)
        cls.y_value = 3e6
        cls.y = (cls.y_value << cls.length_unit).to(u.km)
        cls.z_value = 5.244e6
        cls.z = (cls.z_value << cls.length_unit).to(u.Mm)
        cls.xyz = np.stack([cls.x, cls.y, cls.z])

    def test_unit_errors(self):
        """Test unit errors when dimensionless parameters are used"""

        msg = (
            "'NoneType' object has no attribute 'get_converter'"
            ".*\n.*treated as dimensionless"
        )
        with pytest.raises(AttributeError, match=msg):
            erfa_ufunc.gc2gde(self.equatorial_radius_value, self.flattening, self.xyz)
        with pytest.raises(AttributeError, match=msg):
            erfa_ufunc.gd2gce(
                self.equatorial_radius_value,
                self.flattening,
                self.lon,
                self.lat,
                self.height,
            )
        with pytest.raises(AttributeError, match=msg):
            erfa_ufunc.gc2gde(self.equatorial_radius, self.flattening, self.xyz.value)
        with pytest.raises(AttributeError, match=msg):
            erfa_ufunc.gd2gce(
                self.equatorial_radius,
                self.flattening,
                self.lon_value,
                self.lat,
                self.height,
            )
        with pytest.raises(AttributeError, match=msg):
            erfa_ufunc.gd2gce(
                self.equatorial_radius,
                self.flattening,
                self.lon,
                self.lat,
                self.height_value,
            )

    def test_gc2gde(self):
        """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""

        e, p, h, status = erfa_ufunc.gc2gde(
            self.equatorial_radius, self.flattening, self.xyz
        )
        assert status == 0
        vvd(e, self.lon_value, 1e-14, "eraGc2gde", "e", status)
        vvd(p, self.lat_value, 1e-14, "eraGc2gde", "p", status)
        vvd(h, self.height_value, 1e-8, "eraGc2gde", "h", status)

    def test_gd2gce(self):
        """Test that we reproduce erfa/src/t_erfa_c.c t_gc2gd"""

        xyz, status = erfa_ufunc.gd2gce(
            self.equatorial_radius, self.flattening, self.lon, self.lat, self.height
        )
        assert status == 0
        vvd(xyz[0], self.x_value, 1e-7, "eraGd2gce", "e", status)
        vvd(xyz[1], self.y_value, 1e-7, "eraGd2gce", "p", status)
        vvd(xyz[2], self.z_value, 1e-7, "eraGd2gce", "h", status)
