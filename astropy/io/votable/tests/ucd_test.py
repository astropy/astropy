# Licensed under a 3-clause BSD style license - see LICENSE.rst




from ....tests.helper import raises

# LOCAL
from .. import ucd


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
}


def test_check():
    for s, p in examples.items():
        assert ucd.parse_ucd(s, True, True) == p
        assert ucd.check_ucd(s, True, True)


@raises(ValueError)
def test_too_many_colons():
    ucd.parse_ucd("ivoa:stsci:phot", True, True)


@raises(ValueError)
def test_invalid_namespace():
    ucd.parse_ucd("_ivoa:phot.mag", True, True)


@raises(ValueError)
def test_invalid_word():
    ucd.parse_ucd("-pho")
