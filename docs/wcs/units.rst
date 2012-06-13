.. include:: references.txt

.. _fits-unit:

FITS unit specification
=======================

Supported units
---------------

The following units are supported by the FITS standard:

**SI base & supplementary units**

====================== ============ =================
Quantity               Unit String  Meaning
====================== ============ =================
length                 m            metre
mass                   kg           kilogram
time                   s            second of time
plane angle            rad          radian
solid angle            sr           steradian
temperature            K            kelvin
electric current       A            ampere
amount of substance    mol          mole
luminous intensity     cd           candela
====================== ============ =================

**IAU-recognized derived units**

====================== ============ ================= ==================
Quantity               Unit String  Meaning           Equivalence
====================== ============ ================= ==================
frequency              Hz           hertz             s\ :sup:`-1`
energy                 J            joule             N m
power                  W            watt              J s\ :sup:`-1`
electric potential     V            volts             J C\ :sup:`-1`
force                  N            newton            kg m s\ :sup:`-2`
pressure, stress       Pa           pascal            N m\ :sup:`-2`
electric charge        C            coulomb           A s
electric resistance    ohm          ohm (Ω)           V A\ :sup:`-1`
electric conductance   S            siemens           A V\ :sup:`-1`
electric capacitance   F            farad             C V\ :sup:`-1`
magnetic flux          Wb           weber             V s
magnetic flux density  T            tesla             Wb m\ :sup:`-2`
inductance             H            henry             Wb A\ :sup:`-1`
luminous flux          lm           lumen             cd sr
illuminance            lx           lux               lm m\ :sup:`-2`
====================== ============ ================= ==================

**Additional units**

====================== ============ ======================== ============================================
Quantity               Unit String  Meaning                  Equivalence
====================== ============ ======================== ============================================
mass                   u            unified atomic mass unit 1.6605387 x 10\ :sup:`-27` kg
mass                   solMass      solar mass               1.9891 x 10\ :sup:`30` kg
plane angle            deg          degree of arc            1.745533 x 10\ :sup:`-2` rad
plane angle            arcsec       second of arc            4.848137 x 10\ :sup:`-6` rad
plane angle            arcmin       minute of arc            2.90888 x 10\ :sup:`-4` rad
time                   min          minute
time                   h            hour
time                   d            day                      8.64 x 10\ :sup:`4` s
time                   yr           year (Julian)            3.15576 x10\ :sup:`-7` s (365.25 d)
energy                 eV           electron volt            1.602177 x 10\ :sup:`-19` J
energy                 erg          erg                      10\ :sup:`-7` J
energy                 Ry           Rydberg                  13.605692 eV
length                 angstrom     angstrom                 10\ :sup:`-10` m
length                 AU           astronomical unit        1.49598 x 10\ :sup:`11` m
length                 lyr          light year               9.460530 x 10\ :sup:`-15` m
length                 pc           parsec                   3.0857 x 10\ :sup:`-16` m
length                 solRad       solar radius             6.9599 x 10\ :sup:`8` m
events                 count        counts
events                 photon       photons
flux density           Jy           jansky                   10\ :sup:`-16` W m\ :sup:`-2` Hz\ :sup:`-1`
flux density           mag          (stellar) magnitude
flux density           Crab         'crab'
flux density           beam         beam                     Jy/beam
flux density           solLum       solar luminosity
magnetic field         G            gauss                    10\ :sup:`-4` T
area                   pixel        (image/detector) pixel
area                   voxel        3-d analog of pixel
area                   barn         barn                     10\ :sup:`-28` m\ :sup:`2`
device                 chan         (detector) channel
device                 byte         (computer) byte
device                 bit          (computer) bits
device                 adu          A/D converter units
misc                   bin          numerous applications
misc                   Sun          wrt. sun
====================== ============ ======================== ============================================

Potentially unsafe translations of ``"D"``, ``"H"``, and ``"S"``, are
optional, using the *translate_units* parameter.

.. _unit-aliases:

Unit aliases
------------

When converting non-standard units to standard ones, a case-sensitive
match is required for the aliases listed below, in particular the only
recognized aliases with metric prefixes are ``"KM"``, ``"KHZ"``,
``"MHZ"``, and ``"GHZ"``.

========== =============================================================
Unit       Recognized aliases
========== =============================================================
Angstrom   angstrom
arcmin     arcmins, ARCMIN, ARCMINS
arcsec     arcsecs, ARCSEC, ARCSECS
beam       BEAM
byte       Byte
count      ct
d          day, days, (D), DAY, DAYS
deg        degree, degrees, DEG, DEGREE, DEGREES
GHz        GHZ
h          hr, (H), HR
Hz         hz, HZ
kHz        KHZ
Jy         JY
K          kelvin, kelvins, Kelvin, Kelvins, KELVIN, KELVINS
km         KM
m          metre, meter, metres, meters, M, METRE, METER, METRES, METERS
min        MIN
MHz        MHZ
Ohm        ohm
Pa         pascal, pascals, Pascal, Pascals, PASCAL, PASCALS
photon     ph
pixel      pixels, PIXEL, PIXELS, pix
rad        radian, radians, RAD, RADIAN, RADIANS
s          sec, second, seconds, (S), SEC, SECOND, SECONDS
V          volt, volts, Volt, Volts, VOLT, VOLTS
yr         year, years, YR, YEAR, YEARS
========== =============================================================

The aliases ``"angstrom"``, ``"ohm"``, and ``"Byte"`` for (Angstrom,
Ohm, and byte) are recognized by astropy.wcs/wcslib itself as an
unofficial extension of the standard, but they are converted to the
standard form here.

Prefixes
--------

The following metric prefixes are supported:

======= ======= ===================
Prefix  String  Magnitude
======= ======= ===================
yocto   y       10\ :sup:`-24`
zepto   z       10\ :sup:`-21`
atto    a       10\ :sup:`-18`
femto   f       10\ :sup:`-15`
pico    p       10\ :sup:`-12`
nano    n       10\ :sup:`-9`
micro   u       10\ :sup:`-6`
milli   m       10\ :sup:`-3`
centi   c       10\ :sup:`-2`
deci    d       10\ :sup:`-1`
deka    da      10\ :sup:`1`
hecto   h       10\ :sup:`2`
kilo    k       10\ :sup:`3`
Mega    M       10\ :sup:`6`
Giga    G       10\ :sup:`9`
Tera    T       10\ :sup:`12`
Peta    P       10\ :sup:`15`
Exa     E       10\ :sup:`18`
Zetta   Z       10\ :sup:`21`
Yotta   Y       10\ :sup:`24`
======= ======= ===================

Table 6 of WCS Paper I lists eleven units for which metric prefixes
are allowed.  However, in this implementation only prefixes greater
than unity are allowed for ``"a"`` (annum), ``"yr"`` (year), ``"pc"``
(parsec), ``"bit"``, and ``"byte"``, and only prefixes less than unity
are allowed for ``"mag"`` (stellar magnitude).

Metric prefix ``"P"`` (peta) is specifically forbidden for ``"a"``
(annum) to avoid confusion with ``"Pa"`` (Pascal, not peta-annum).
Note that metric prefixes are specifically disallowed for ``"h"``
(hour) and ``"d"`` (day) so that ``"ph"`` (photons) cannot be
interpreted as pico-hours, nor ``"cd"`` (candela) as centi-days.

Operators
---------

A compound unit is considered to be formed by a series of sub-strings
of component units & mathematical operations. Each of these
sub-strings must be separated by at least one space or a mathematical
operator (``*`` or ``/``).

Multiplication
^^^^^^^^^^^^^^

Multiplicative units can be specified either:

- by simply using one or more preceding spaces, e.g. ``str1 str2``
  (The recommended method).

- by the use of a single asterisk (``*``) with optional whitespace,
  e.g. ``str1 * str2``.

Division
^^^^^^^^

Units which form the denominator of a compound expression can be
specified either:

- by using a slash (``/``) with optional whitespace, e.g. ``str1 /
  str2``.  If such a syntax is used, it is recommended that no space
  is included between the slash and the unit string.

- by raising a multiplicative unit to a negative power (see below).

It should be stressed that the slash character only effects the
sub-string it immediately precedes. Thus, unless brackets are used,
subsequent sub-strings which also form part of the denominator of the
compound expression must also be preceded by a slash.  For example,
``str1 /str2 str3`` is equivalent to ``str1 str3 /str2`` whilst ``str1
/str2 /str3`` is equivalent to ``str1 /(str2 * str3)``.

Raising to Powers
^^^^^^^^^^^^^^^^^

A unit string raised to the power *y* is specified:

- by using two asterisks (``**``) followed by the index enclosed
  within round brackets and with no preceding or intervening spaces,
  e.g. ``str1**(y)`` or ``str1**(-y)``.

However, if *y* is positive, then the brackets need not be included,
but a following space is recommended if additional sub-strings follow.

Use of brackets
^^^^^^^^^^^^^^^

Any number of pairs of round brackets (``()``) may be used within the
string for a compound unit in order to prevent ambiguities. As
described within this section, a number of rules always/often require
their use. However, it is suggested that conservative use is made of
such pairs of brackets in order to minimize the total length of
compound strings.  (It should be remembered that a maximum of 68
characters are allowed in the card image of keywords.)

Avoidance of underflows & overflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The inclusion of numerical factors within the unit string should
generally be avoided (by the use of multiples and/or submultiples of
component basic units).

However, occasionally it may be preferable to include such factors on
the grounds of user-friendliness and/or to minimize the risk of
computer under- or overflows. In such cases, the numerical factor can
simply be considered a basic unit string.

The following additional guidelines are suggested:

- the numerical factor should precede any unit strings

- only powers of 10 are used as numerical factors

Mathematical Operations & Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A number of mathematical operations are supported.  It should be noted
that the (round) brackets are mandatory in all cases in which they are
included in the table.

=============== ==================================================
String          Meaning
=============== ==================================================
``str1*str2``   Multiplication
``str1 /str2``  Division
``str1**(y)``   Raised to the power *y*
``log(str1)``   Common Logarithm (to base 10)
``ln(str1)``    Natural Logarithm
``exp(str1)``   Exponential (exp\ :sup:`str1`\ )
``sqrt(str1)``  Square root
``sin(str1)``   Sine
``cos(str1)``   Cosine
``tan(str1)``   Tangent
``asin(str1)``  Arc Sine
``acos(str1)``  Arc Cosine
``atan(str1)``  Arc Tangent
``sinh(str1)``  Hyperbolic Sine
``cosh(str1)``  Hyperbolic Cosine
``tanh(str1)``  Hyperbolic Tangent
=============== ==================================================

Function types ``log()``, ``ln()`` and ``exp()`` may only occur at the
start of the units specification.
