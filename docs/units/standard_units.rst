Standard units
==============

Standard units are defined in the `astropy.units` package as object
instances.

All units are defined in term of basic 'irreducible' units. The
irreducible units include:

  - Length (meter)
  - Time (second)
  - Mass (kilogram)
  - Current (ampere)
  - Temperature (Kelvin)
  - Angular distance (radian)
  - Solid angle (steradian)
  - Luminous intensity (candela)
  - Stellar magnitude (mag)
  - Amount of substance (mole)
  - Photon count (photon)

(There are also some more obscure base units required by the FITS
standard that are no longer recommended for use.)

Units that involve combinations of fundamental units are instances of
`~astropy.units.core.CompositeUnit`. In most cases, one does not need
to worry about the various kinds of unit classes unless one wants to
design a more complex case.

There are many units already predefined in the module. One may use the
following function to list all the existing predefined units of a
given type::

  >>> u.equivalent_units(u.g)
  Primary name | Unit definition | Aliases
  g            | 1.00e-03 kg     | gram
  kg           | irreducible     | kilogram
  lb           | 4.54e-01 kg     | pound
  m_e          | 9.11e-31 kg     |
  m_p          | 1.67e-27 kg     |
  oz           | 2.83e-02 kg     | ounce
  solMass      | 1.99e+30 kg     |
  t            | 1.00e+03 kg     | tonne
  ton          | 9.07e+02 kg     |
  u            | 1.66e-27 kg     |
