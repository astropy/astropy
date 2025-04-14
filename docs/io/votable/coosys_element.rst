COOSYS element
--------------

The ``COOSYS`` element is allowed in the ``VOTABLE`` and ``RESOURCE`` elements. It
describes a coordinate system and has the following attributes:

- ``ID``: a label,
- ``system``: the reference frame of the coordinates. Since VOTable version 1.5, it
  follows the `IVOA reference frame vocabulary <http://www.ivoa.net/rdf/refframe>`_.
- ``equinox``: fixes the equatorial or ecliptic systems
- ``epoch``: the epoch of the positions, as Julian or Besselian years
- ``refposition``: it exists since VOTable version 1.5 and follows the `IVOA reference
  positions vocabulary <http://www.ivoa.net/rdf/refposition>`_.

They are all optional. The `~astropy.io.votable.tree.CooSys` class represents a
``COOSYS`` element. It can be created from scratch::

    >>> from astropy.io.votable.tree import CooSys
    >>> coosys = CooSys(system="FK4", equinox="B1950", epoch="B1980")

and converted into an astropy built-in frame::

    >>> coosys.to_astropy_frame()
    <FK4 Frame (equinox=B1950.000, obstime=B1950.000)>

Using the COOSYS element
------------------------

The ``COOSYS`` element is often used to describe the coordinates of a ``FIELD`` element.
In that case, the ``FIELD`` will have a ``ref`` attribute corresponding to the
``COOSYS``'s '``ID`` field. See for example this VOTable::

    >>> from astropy.io.votable import parse as parse_votable
    >>> from io import BytesIO
    >>> votable = parse_votable(BytesIO(str.encode("""<?xml version="1.0" encoding="utf-8"?>
    ...    <VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 http://www.ivoa.net/xml/VOTable/VOTable-1.3.xsd">
    ...    <RESOURCE type="results">
    ...      <COOSYS ID="coosys_c1356galcand" epoch="2016.0" system="ICRS"/>
    ...      <TABLE ID="result_S1744377631073" name="result_S1744377631073" nrows="3">
    ...      <FIELD ID="Source" datatype="long" name="Source" ucd="meta.id;meta.main">
    ...        <DESCRIPTION>Unique source identifier (unique within a particular Data Release) (source_id)</DESCRIPTION></FIELD>
    ...      <FIELD ID="RA_ICRS" datatype="double" name="RA_ICRS" ref="coosys_c1356galcand" ucd="pos.eq.ra;meta.main" unit="deg">
    ...        <DESCRIPTION>Right ascension (ICRS) at Ep=2016.0 (ra)</DESCRIPTION></FIELD>
    ...      <FIELD ID="DE_ICRS" datatype="double" name="DE_ICRS" ref="coosys_c1356galcand" ucd="pos.eq.dec;meta.main" unit="deg">
    ...        <DESCRIPTION>Declination (ICRS) at Ep=2016.0 (dec)</DESCRIPTION></FIELD>
    ...        <DATA><TABLEDATA>
    ...        <TR><TD>3680787975936</TD><TD>45.01912373733</TD><TD>0.16323692972</TD></TR>
    ...        <TR><TD>5570573422080</TD><TD>44.85334508786</TD><TD>0.12810861107</TD></TR>
    ...        <TR><TD>7941395337344</TD><TD>44.97411533912</TD><TD>0.22751845955</TD></TR>
    ...        </TABLEDATA></DATA>
    ...      </TABLE></RESOURCE></VOTABLE>""")))
    >>> table = votable.get_first_table()
    >>> for field in votable.iter_fields_and_params():
    ...     print(field)
    <FIELD ID="Source" datatype="long" name="Source" ucd="meta.id;meta.main"/>
    <FIELD ID="RA_ICRS" datatype="double" name="RA_ICRS" ref="coosys_c1356galcand" ucd="pos.eq.ra;meta.main" unit="deg"/>
    <FIELD ID="DE_ICRS" datatype="double" name="DE_ICRS" ref="coosys_c1356galcand" ucd="pos.eq.dec;meta.main" unit="deg"/>

We can now generate a ``SkyCoord`` object from the VOTable information::

   >>> from astropy.coordinates import SkyCoord
   >>> SkyCoord(table.array.data["RA_ICRS"],
   ...          table.array.data["DE_ICRS"],
   ...          frame=votable.get_coosys_by_id(
   ...                         table.get_field_by_id("RA_ICRS").ref
   ...                         ).to_astropy_frame(),
   ...          unit="deg")
   <SkyCoord (ICRS): (ra, dec) in deg
      [(45.01912374, 0.16323693), (44.85334509, 0.12810861),
       (44.97411534, 0.22751846)]>
