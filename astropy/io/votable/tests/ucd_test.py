# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.io.votable import ucd


def test_none():
    assert ucd.check_ucd(None)


examples = {
    'phys.temperature':
        [('ivoa', 'phys.temperature')],
    'pos.eq.ra;meta.main':
        [('ivoa', 'pos.eq.ra'), ('ivoa', 'meta.main')],
    'meta.id;src':
        [('ivoa', 'meta.id'), ('ivoa', 'src')],
    'phot.flux;em.radio;arith.ratio':
        [('ivoa', 'phot.flux'), ('ivoa', 'em.radio'), ('ivoa', 'arith.ratio')],
    'PHot.Flux;EM.Radio;ivoa:arith.Ratio':
        [('ivoa', 'phot.flux'), ('ivoa', 'em.radio'), ('ivoa', 'arith.ratio')],
    'pos.galactic.lat':
        [('ivoa', 'pos.galactic.lat')],
    'meta.code;phot.mag':
        [('ivoa', 'meta.code'), ('ivoa', 'phot.mag')],
    'stat.error;phot.mag':
        [('ivoa', 'stat.error'), ('ivoa', 'phot.mag')],
    'phys.temperature;instr;stat.max':
        [('ivoa', 'phys.temperature'), ('ivoa', 'instr'),
         ('ivoa', 'stat.max')],
    'stat.error;phot.mag;em.opt.V':
        [('ivoa', 'stat.error'), ('ivoa', 'phot.mag'), ('ivoa', 'em.opt.V')],
    'phot.color;em.opt.B;em.opt.V':
        [('ivoa', 'phot.color'), ('ivoa', 'em.opt.B'), ('ivoa', 'em.opt.V')],
    'stat.error;phot.color;em.opt.B;em.opt.V':
        [('ivoa', 'stat.error'), ('ivoa', 'phot.color'), ('ivoa', 'em.opt.B'),
         ('ivoa', 'em.opt.V')],
}


def test_check():
    for s, p in examples.items():
        assert ucd.parse_ucd(s, True, True) == p
        assert ucd.check_ucd(s, True, True)


def test_too_many_colons():
    with pytest.raises(ValueError):
        ucd.parse_ucd("ivoa:stsci:phot", True, True)


def test_invalid_namespace():
    with pytest.raises(ValueError):
        ucd.parse_ucd("_ivoa:phot.mag", True, True)


def test_invalid_word():
    with pytest.raises(ValueError):
        ucd.parse_ucd("-pho")
