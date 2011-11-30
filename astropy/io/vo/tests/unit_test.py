from .. import unit


def test_units():
    strings = [
        # These first few are from the spec document itself
        ('mW/m2'       , True),
        ('0.1nm'       , True),
        ('solMass3/2'  , False),
        ('[solMass]'   , True),
        ('0.1 nm'      , False),
        ('km/s'        , True),
        ('km / s'      , False),
        ('km s-1.'     , False),
        ('10pix/nm'    , True),
        ('pix/0.1nm'   , False),
        ('pix/(0.1nm)' , False),
        ('1.5x10+11m'  , True),
        ('kW.h'        , True),
        ('kWh'         , False),
        ('m2'          , True),
        ('uarcsec'     , True),

        # These are additional tests
        ('w'                  , False),
        ('+1.0m/s'            , True),
        ('1.0x10+11m/s+10+11' , True),
        ('mm'                 , True),
        ('Om'                 , False)
        ]

    def run(args):
        s, correct = args
        print s, correct
        assert unit.is_unit(s) == correct

    for s, correct in strings:
        yield run, (s, correct)


def test_basic_units():
    def run(args):
        prefix, u = args
        print prefix, u
        assert unit.is_unit(prefix + u)
        assert unit.is_unit('[' + prefix + u + ']')

    for u in unit.unit_names:
        yield run, ('', u)

    for p in unit.unit_prefixes:
        for u in unit.unit_names:
            yield run, (p, u)

