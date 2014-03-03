.. include:: references.txt

.. _relax:

Relax constants
===============

The ``relax`` keyword argument controls the handling of non-standard
FITS WCS keywords.

Note that the default value of ``relax`` is `True` for reading (to
accept all non standard keywords), and `False` for writing (to write
out only standard keywords), in accordance with `Postel's prescription
<http://catb.org/jargon/html/P/Postels-Prescription.html>`_:

    “Be liberal in what you accept, and conservative in what you send.”

.. _relaxread:

Header-reading relaxation constants
-----------------------------------

`~astropy.wcs.WCS`, `~astropy.wcs.Wcsprm` and
`~astropy.wcs.find_all_wcs` have a *relax* argument, which may be
either `True`, `False` or an `int`.

- If `True`, (default), all non-standard WCS extensions recognized by the parser
  will be handled.

- If `False`, none of the extensions (even those in the
  errata) will be handled.  Non-conformant keywords will be handled in
  the same way as non-WCS keywords in the header, i.e. by simply
  ignoring them.

- If an `int`, is is a bit field to provide fine-grained control over
  what non-standard WCS keywords to accept.  The flag bits are subject
  to change in future and should be set by using the constants
  beginning with ``WCSHDR_`` in the `astropy.wcs` module.

  For example, to accept ``CD00i00j`` and ``PC00i00j`` use::

      relax = astropy.wcs.WCSHDR_CD00i00j | astropy.wcs.WCSHDR_PC00i00j

  The parser always treats ``EPOCH`` as subordinate to ``EQUINOXa`` if
  both are present, and ``VSOURCEa`` is always subordinate to
  ``ZSOURCEa``.

  Likewise, ``VELREF`` is subordinate to the formalism of WCS Paper
  III.

The flag bits are:

- ``WCSHDR_none``: Don't accept any extensions (not even those in the
  errata).  Treat non-conformant keywords in the same way as non-WCS
  keywords in the header, i.e. simply ignore them.  (This is
  equivalent to passing `False`)

- ``WCSHDR_all``: Accept all extensions recognized by the parser.  (This
  is equivalent to the default behavior or passing `True`).

- ``WCSHDR_CROTAia``: Accept ``CROTAia``, ``iCROTna``, ``TCROTna``
- ``WCSHDR_EPOCHa``:  Accept ``EPOCHa``.
- ``WCSHDR_VELREFa``: Accept ``VELREFa``.

        The constructor always recognizes the AIPS-convention
        keywords, ``CROTAn``, ``EPOCH``, and ``VELREF`` for the
        primary representation ``(a = ' ')`` but alternates are
        non-standard.

        The constructor accepts ``EPOCHa`` and ``VELREFa`` only if
        ``WCSHDR_AUXIMG`` is also enabled.

- ``WCSHDR_CD00i00j``: Accept ``CD00i00j``.
- ``WCSHDR_PC00i00j``: Accept ``PC00i00j``.
- ``WCSHDR_PROJPn``: Accept ``PROJPn``.

        These appeared in early drafts of WCS Paper I+II (before they
        were split) and are equivalent to ``CDi_ja``, ``PCi_ja``, and
        ``PVi_ma`` for the primary representation ``(a = ' ')``.
        ``PROJPn`` is equivalent to ``PVi_ma`` with ``m`` = ``n`` <=
        9, and is associated exclusively with the latitude axis.

- ``WCSHDR_RADECSYS``: Accept ``RADECSYS``.  This appeared in early
  drafts of WCS Paper I+II and was subsequently replaced by
  ``RADESYSa``.  The construtor accepts ``RADECSYS`` only if
  ``WCSHDR_AUXIMG`` is also enabled.

- ``WCSHDR_VSOURCE``: Accept ``VSOURCEa`` or ``VSOUna``.  This appeared
  in early drafts of WCS Paper III and was subsequently dropped in
  favour of ``ZSOURCEa`` and ``ZSOUna``.  The constructor accepts
  ``VSOURCEa`` only if ``WCSHDR_AUXIMG`` is also enabled.

- ``WCSHDR_DOBSn``: Allow ``DOBSn``, the column-specific analogue of
  ``DATE-OBS``.  By an oversight this was never formally defined in
  the standard.

- ``WCSHDR_LONGKEY``: Accept long forms of the alternate binary table
  and pixel list WCS keywords, i.e. with "a" non- blank.
  Specifically::

        jCRPXna  TCRPXna  :  jCRPXn  jCRPna  TCRPXn  TCRPna  CRPIXja
           -     TPCn_ka  :    -     ijPCna    -     TPn_ka  PCi_ja
           -     TCDn_ka  :    -     ijCDna    -     TCn_ka  CDi_ja
        iCDLTna  TCDLTna  :  iCDLTn  iCDEna  TCDLTn  TCDEna  CDELTia
        iCUNIna  TCUNIna  :  iCUNIn  iCUNna  TCUNIn  TCUNna  CUNITia
        iCTYPna  TCTYPna  :  iCTYPn  iCTYna  TCTYPn  TCTYna  CTYPEia
        iCRVLna  TCRVLna  :  iCRVLn  iCRVna  TCRVLn  TCRVna  CRVALia
        iPVn_ma  TPVn_ma  :    -     iVn_ma    -     TVn_ma  PVi_ma
        iPSn_ma  TPSn_ma  :    -     iSn_ma    -     TSn_ma  PSi_ma

  where the primary and standard alternate forms together with the
  image-header equivalent are shown rightwards of the colon.

  The long form of these keywords could be described as quasi-
  standard.  ``TPCn_ka``, ``iPVn_ma``, and ``TPVn_ma`` appeared by
  mistake in the examples in WCS Paper II and subsequently these and
  also ``TCDn_ka``, ``iPSn_ma`` and ``TPSn_ma`` were legitimized by
  the errata to the WCS papers.

  Strictly speaking, the other long forms are non-standard and in fact
  have never appeared in any draft of the WCS papers nor in the
  errata.  However, as natural extensions of the primary form they are
  unlikely to be written with any other intention.  Thus it should be
  safe to accept them provided, of course, that the resulting keyword
  does not exceed the 8-character limit.

  If ``WCSHDR_CNAMn`` is enabled then also accept::

        iCNAMna  TCNAMna  :   ---   iCNAna    ---   TCNAna  CNAMEia
        iCRDEna  TCRDEna  :   ---   iCRDna    ---   TCRDna  CRDERia
        iCSYEna  TCSYEna  :   ---   iCSYna    ---   TCSYna  CSYERia

  Note that ``CNAMEia``, ``CRDERia``, ``CSYERia``, and their variants
  are not used by `astropy.wcs` but are stored as auxiliary information.

- ``WCSHDR_CNAMn``: Accept ``iCNAMn``, ``iCRDEn``, ``iCSYEn``,
  ``TCNAMn``, ``TCRDEn``, and ``TCSYEn``, i.e. with ``a`` blank.
  While non-standard, these are the obvious analogues of ``iCTYPn``,
  ``TCTYPn``, etc.

- ``WCSHDR_AUXIMG``: Allow the image-header form of an auxiliary WCS
  keyword with representation-wide scope to provide a default value
  for all images.  This default may be overridden by the
  column-specific form of the keyword.

  For example, a keyword like ``EQUINOXa`` would apply to all image
  arrays in a binary table, or all pixel list columns with alternate
  representation ``a`` unless overridden by ``EQUIna``.

  Specifically the keywords are::

        LATPOLEa  for LATPna
        LONPOLEa  for LONPna
        RESTFREQ  for RFRQna
        RESTFRQa  for RFRQna
        RESTWAVa  for RWAVna

  whose keyvalues are actually used by WCSLIB, and also keywords that
  provide auxiliary information that is simply stored in the wcsprm
  struct::

        EPOCH         -       ... (No column-specific form.)
        EPOCHa        -       ... Only if WCSHDR_EPOCHa is set.
        EQUINOXa  for EQUIna
        RADESYSa  for RADEna
        RADECSYS  for RADEna  ... Only if WCSHDR_RADECSYS is set.
        SPECSYSa  for SPECna
        SSYSOBSa  for SOBSna
        SSYSSRCa  for SSRCna
        VELOSYSa  for VSYSna
        VELANGLa  for VANGna
        VELREF        -       ... (No column-specific form.)
        VELREFa       -       ... Only if WCSHDR_VELREFa is set.
        VSOURCEa  for VSOUna  ... Only if WCSHDR_VSOURCE is set.
        WCSNAMEa  for WCSNna  ... Or TWCSna (see below).
        ZSOURCEa  for ZSOUna

        DATE-AVG  for DAVGn
        DATE-OBS  for DOBSn
        MJD-AVG   for MJDAn
        MJD-OBS   for MJDOBn
        OBSGEO-X  for OBSGXn
        OBSGEO-Y  for OBSGYn
        OBSGEO-Z  for OBSGZn

  where the image-header keywords on the left provide default values
  for the column specific keywords on the right.

  Keywords in the last group, such as ``MJD-OBS``, apply to all
  alternate representations, so ``MJD-OBS`` would provide a default
  value for all images in the header.

  This auxiliary inheritance mechanism applies to binary table image
  arrays and pixel lists alike.  Most of these keywords have no
  default value, the exceptions being ``LONPOLEa`` and ``LATPOLEa``,
  and also ``RADESYSa`` and ``EQUINOXa`` which provide defaults for
  each other.  Thus the only potential difficulty in using
  ``WCSHDR_AUXIMG`` is that of erroneously inheriting one of these four
  keywords.

  Unlike ``WCSHDR_ALLIMG``, the existence of one (or all) of these
  auxiliary WCS image header keywords will not by itself cause a
  `~astropy.wcs.Wcsprm` object to be created for alternate
  representation ``a``.  This is because they do not provide
  sufficient information to create a non-trivial coordinate
  representation when used in conjunction with the default values of
  those keywords, such as ``CTYPEia``, that are parameterized by axis
  number.

- ``WCSHDR_ALLIMG``: Allow the image-header form of *all* image header
  WCS keywords to provide a default value for all image arrays in a
  binary table (n.b. not pixel list).  This default may be overridden
  by the column-specific form of the keyword.

  For example, a keyword like ``CRPIXja`` would apply to all image
  arrays in a binary table with alternate representation ``a``
  unless overridden by ``jCRPna``.

  Specifically the keywords are those listed above for ``WCSHDR_AUXIMG``
  plus::

        WCSAXESa  for WCAXna

  which defines the coordinate dimensionality, and the following
  keywords which are parameterized by axis number::

        CRPIXja   for jCRPna
        PCi_ja    for ijPCna
        CDi_ja    for ijCDna
        CDELTia   for iCDEna
        CROTAi    for iCROTn
        CROTAia        -      ... Only if WCSHDR_CROTAia is set.
        CUNITia   for iCUNna
        CTYPEia   for iCTYna
        CRVALia   for iCRVna
        PVi_ma    for iVn_ma
        PSi_ma    for iSn_ma

        CNAMEia   for iCNAna
        CRDERia   for iCRDna
        CSYERia   for iCSYna

  where the image-header keywords on the left provide default values
  for the column specific keywords on the right.

  This full inheritance mechanism only applies to binary table image
  arrays, not pixel lists, because in the latter case there is no
  well-defined association between coordinate axis number and column
  number.

  Note that ``CNAMEia``, ``CRDERia``, ``CSYERia``, and their variants
  are not used by pywcs but are stored in the `~astropy.wcs.Wcsprm`
  object as auxiliary information.

  Note especially that at least one `~astropy.wcs.Wcsprm` object will
  be returned for each ``a`` found in one of the image header keywords
  listed above:

    - If the image header keywords for ``a`` **are not** inherited by
      a binary table, then the struct will not be associated with any
      particular table column number and it is up to the user to
      provide an association.

    - If the image header keywords for ``a`` **are** inherited by a
      binary table image array, then those keywords are considered to
      be "exhausted" and do not result in a separate
      `~astropy.wcs.Wcsprm` object.

.. _relaxwrite:

Header-writing relaxation constants
-----------------------------------

`~astropy.wcs.wcs.WCS.to_header` and `~astropy.wcs.wcs.WCS.to_header_string`
has a *relax* argument which may be either `True`, `False` or an
`int`.

- If `True`, write all recognized extensions.

- If `False` (default), write all extensions that are considered to be
  safe and recommended, equivalent to ``WCSHDO_safe`` (described below).

- If an `int`, is is a bit field to provide fine-grained control over
  what non-standard WCS keywords to accept.  The flag bits are subject
  to change in future and should be set by using the constants
  beginning with ``WCSHDO_`` in the `astropy.wcs` module.

The flag bits are:

- ``WCSHDO_none``: Don't use any extensions.

- ``WCSHDO_all``: Write all recognized extensions, equivalent to setting
  each flag bit.

- ``WCSHDO_safe``: Write all extensions that are considered to be safe
  and recommended.

- ``WCSHDO_DOBSn``: Write ``DOBSn``, the column-specific analogue of
  ``DATE-OBS`` for use in binary tables and pixel lists.  WCS Paper
  III introduced ``DATE-AVG`` and ``DAVGn`` but by an oversight
  ``DOBSn`` (the obvious analogy) was never formally defined by the
  standard.  The alternative to using ``DOBSn`` is to write
  ``DATE-OBS`` which applies to the whole table.  This usage is
  considered to be safe and is recommended.

- ``WCSHDO_TPCn_ka``: WCS Paper I defined

  - ``TPn_ka`` and ``TCn_ka`` for pixel lists

    but WCS Paper II uses ``TPCn_ka`` in one example and subsequently
    the errata for the WCS papers legitimized the use of

  - ``TPCn_ka`` and ``TCDn_ka`` for pixel lists

    provided that the keyword does not exceed eight characters.  This
    usage is considered to be safe and is recommended because of the
    non-mnemonic terseness of the shorter forms.

- ``WCSHDO_PVn_ma``: WCS Paper I defined

  - ``iVn_ma`` and ``iSn_ma`` for bintables and
  - ``TVn_ma`` and ``TSn_ma`` for pixel lists

    but WCS Paper II uses ``iPVn_ma`` and ``TPVn_ma`` in the examples
    and subsequently the errata for the WCS papers legitimized the use
    of

  - ``iPVn_ma`` and ``iPSn_ma`` for bintables and
  - ``TPVn_ma`` and ``TPSn_ma`` for pixel lists

    provided that the keyword does not exceed eight characters.  This
    usage is considered to be safe and is recommended because of the
    non-mnemonic terseness of the shorter forms.

- ``WCSHDO_CRPXna``: For historical reasons WCS Paper I defined

  - ``jCRPXn``, ``iCDLTn``, ``iCUNIn``, ``iCTYPn``, and ``iCRVLn`` for
    bintables and
  - ``TCRPXn``, ``TCDLTn``, ``TCUNIn``, ``TCTYPn``, and ``TCRVLn`` for
    pixel lists

    for use without an alternate version specifier.  However, because
    of the eight-character keyword constraint, in order to accommodate
    column numbers greater than 99 WCS Paper I also defined

  - ``jCRPna``, ``iCDEna``, ``iCUNna``, ``iCTYna`` and ``iCRVna`` for
    bintables and
  - ``TCRPna``, ``TCDEna``, ``TCUNna``, ``TCTYna`` and ``TCRVna`` for
    pixel lists

    for use with an alternate version specifier (the ``a``).  Like the
    ``PC``, ``CD``, ``PV``, and ``PS`` keywords there is an obvious
    tendency to confuse these two forms for column numbers up to 99.
    It is very unlikely that any parser would reject keywords in the
    first set with a non-blank alternate version specifier so this
    usage is considered to be safe and is recommended.

- ``WCSHDO_CNAMna``: WCS Papers I and III defined

  - ``iCNAna``,  ``iCRDna``,  and ``iCSYna``  for bintables and
  - ``TCNAna``,  ``TCRDna``,  and ``TCSYna``  for pixel lists

    By analogy with the above, the long forms would be

  - ``iCNAMna``, ``iCRDEna``, and ``iCSYEna`` for bintables and
  - ``TCNAMna``, ``TCRDEna``, and ``TCSYEna`` for pixel lists

    Note that these keywords provide auxiliary information only, none
    of them are needed to compute world coordinates.  This usage is
    potentially unsafe and is not recommended at this time.

- ``WCSHDO_WCSNna``: Write ``WCSNna`` instead of ``TWCSna`` for pixel
  lists.  While the constructor treats ``WCSNna`` and ``TWCSna`` as
  equivalent, other parsers may not.  Consequently, this usage is
  potentially unsafe and is not recommended at this time.

- ``WCSHDO_SIP``: Write out Simple Imaging Polynomial (SIP) keywords.
