# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines constants used in `astropy.vo.samp`."""

__all__ = ['SAMP_STATUS_OK', 'SAMP_STATUS_WARNING', 'SAMP_STATUS_ERROR',
           'SAMP_HUB_SINGLE_INSTANCE', 'SAMP_HUB_MULTIPLE_INSTANCE',
           'SAMP_RESTRICT_GROUP', 'SAMP_RESTRICT_OWNER',
           'SAFE_MTYPES', 'SAMPY_ICON']

__release__ = "1.3.0"
__profile_version__ = "1.3"

#: General constant for samp.ok status string
SAMP_STATUS_OK = "samp.ok"
#: General constant for samp.warning status string
SAMP_STATUS_WARNING = "samp.warning"
#: General constant for samp.error status string
SAMP_STATUS_ERROR = "samp.error"

#: General constant to specify single instance Hub running mode
SAMP_HUB_SINGLE_INSTANCE = "single"
#: General constant to specify multiple instance Hub running mode
SAMP_HUB_MULTIPLE_INSTANCE = "multiple"

#: General constant to specify the access restriction (through Basic Authentication) to the GROUP
SAMP_RESTRICT_GROUP = "GROUP"
#: General constant to specify the access restriction (through Basic Authentication) to the OWNER
SAMP_RESTRICT_OWNER = "OWNER"

SAFE_MTYPES = ["samp.app.*", "samp.msg.progress", "table.*", "image.*",
               "coord.*", "spectrum.*", "bibcode.*", "voresource.*"]

SAMPY_ICON = \
    b"""iVBORw0KGgoAAAANSUhEUgAAABgAAAAWCAYAAADafVyIAAAABmJLR0QA/wD/AP+gvaeTAAAACXBI
WXMAAA7DAAAOwwHHb6hkAAAAB3RJTUUH2AgFDTkeyUSsmwAAA5FJREFUSMellH9M1GUcx1/PXTfc
MfGGguR1IQeGdZN98ThFSMGmJLIxpFwF0UgXMZXJYm3Qyrm2XLW5wtVGKx1ZOWTgmIgZ5mamaePX
ibsRQQeJ0KEhxK9BePf01133xWOgvv/67nk+3/freT4/HiGlZLZGJkfl0csnaOqxMzo1TvjipZjD
VrIuSiE1Nkks0gWxUInZgN/+6pY5X+7hzthQwB+Cg/Q8t/pZdm98BWtknJi93zXolKuWm0VAgNvj
IaP8VekY6ATAaIhgzRNPo9Pq6B3qo/fvPsamxn3xislCyfOFpMYm+Qzfq/tYFmzKxRRqFACP+dNb
/2z3mW98aj3f7P5MaDUa1QmbeuzyUtc1zjsuYe9zkPdVEYrJIou3FpBoXpvwx51e+kdcmEKN3Afw
mgOkPZPCbHMAW5QibFEKJWmFdA06ZW1LA9XN9eQf2w/QnGi2qlKnAui0uv9r4eqet4CrlptF6fYi
SrcX4RjolACWFbGquqiOmBht9X1XN9Vz3vGTXGi3WFbEitnm9wGiwyJFcowNgBn3DLsq36K09gM5
NDG8YNC8ber657bc9kkOQxPDvrVFuiCy4tPJit9GcoxNPBIAoO9uv3zj67dVRfcqbPFSMpU0MtZs
wRaliIcCAEzNTHO4sUJWXjnJ1Mx0wBijIYL85Jd4eV0WBn2IeCCAf8qOXDhKXds51ZD5y5vCwtTX
iA6LFA8E8GpietJaf72x+fjVGm7c6ggYo9VoyNuwk9L0fQnBQfoWH6Cu7Zysbj7NyOQoANbIOF6w
ZqCYLAGvfeNWhzx+tYa6tu8Dpm/lMhOVu8qJDosUoqb5jCyuOhDwRFstKZSl78P/8fLX0MSwPPxD
BSd+PYXb41HtRSwJ5+z+bxHpn+bKua7svfbOhEwKNuXOCeoadMp3Tn3INWeLaj1vw4uIlI92yJt3
B8hU0vj33gz11xvnhKXGJpGppJEUY8NoiFDB3B4Pe78rkw3tP6pSJcrPfi6P/VzFxXfPYNCHiAsd
l2XJyYOqQZurRVc/HkPEknAA+oddXOz8RRWjmCyIe243rUfel7big8K/NYurDnClu4lH0aHsMjRa
jQZbTAK4XNKvQKLqzQpxKLsMgz7kocwz4raQsz5bzDsHY1PjX5y2NxbUtjbQ1GOf1zg4SM/ezfns
2fy60Go0Cx80L8x+01Hw+6CTrttO2v2678lQI4nmtWTFp6uejf8A94dzfMTZww8AAAAASUVORK5C
YII="""

# TODO: document this global variable. Is this the right place for it?
_THREAD_STARTED_COUNT = 0
