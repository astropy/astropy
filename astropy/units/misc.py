###########################################################################
#UNITS MOVED FROM ASTROPHY.PY TO THIS FILE

###########################################################################
# AREAS

def_unit(['barn', 'barn'], 10 ** -28 * si.m ** 2, namespace=_ns, prefixes=True,
    doc="barn: unit of area used in HEP")    

##########################################################################
# PRESSURE
def_unit(['bar'], 1e5 * si.Pa, namespace=_ns,
         prefixes=[(['m'], ['milli'], 1.e-3)],
         doc="bar: pressure")

##########################################################################
def_unit((['bit', 'b'], ['bit']), namespace=_ns,
         prefixes=si_prefixes + binary_prefixes)
def_unit((['byte', 'B'], ['byte']), 8 * bit, namespace=_ns,
         format={'vounit': 'byte'},
         prefixes=si_prefixes + binary_prefixes,
         exclude_prefixes=['d'])
def_unit(['bin'], namespace=_ns, prefixes=True)

###########################################################################
# ANGULAR MEASUREMENTS

def_unit(['cycle', 'cy'], 2.0 * _numpy.pi * si.rad,
         namespace=_ns, prefixes=False,
         doc="cycle: angular measurement, a full turn or rotation")
def_unit(['spat', 'sp'], 4.0 * _numpy.pi * si.sr,
         namespace=_ns, prefixes=False,
         doc="spat: the solid angle of the sphere, 4pi sr")


###########################################################################
# The torr is almost the same as mmHg but not quite.
# See https://en.wikipedia.org/wiki/Torr
# It may be moved if more similar units are created later.
def_unit(['Torr', 'torr'], _si.atm.value/760. * si.Pa, namespace=_ns,
         prefixes=[(['m'], ['milli'], 1.e-3)],
         doc="Unit of pressure based on an absolute scale, now defined as "
             "exactly 1/760 of a standard atmosphere")

# Unified atomic mass unit
def_unit(['u', 'Da', 'Dalton'], _si.u, namespace=_ns,
         prefixes=True, exclude_prefixes=['a', 'da'],
         doc="Unified atomic mass unit")

######################################################
#EVENTS

def_unit((['pix', 'pixel'], ['pixel']),
         format={'ogip': 'pixel', 'vounit': 'pixel'},
         namespace=_ns, prefixes=True)    