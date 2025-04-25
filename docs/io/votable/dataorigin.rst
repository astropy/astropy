Introduction
------------

Extract basic provenance information from VOTable header. The information is described in
DataOrigin IVOA note: https://www.ivoa.net/documents/DataOrigin/.

DataOrigin includes both the query information (such as publisher, contact, versions, etc.)
and the Dataset origin (such as Creator, bibliographic links, URL, etc.)

This API retrieves Metadata from INFO in VOTable.


Getting Started
---------------

To extract DataOrigin from VOTable

Example: VizieR catalogue J/AJ/167/18

.. doctest-remote-data-all::

    >>> from astropy.io.votable import parse
    >>> from astropy.io.votable.dataorigin import extract_data_origin
    >>> votable = parse("https://vizier.cds.unistra.fr/viz-bin/conesearch/J/AJ/167/18/table4?RA=265.51&DEC=-22.71&SR=0.1")
    >>> data_origin = extract_data_origin(votable)
    >>> print(data_origin)  # doctest: +IGNORE_OUTPUT
    publisher: CDS
    server_software: 7.4.5
    service_protocol: ivo://ivoa.net/std/ConeSearch/v1.03
    request_date: 2025-03-05T14:18:05
    contact: cds-question@unistra.fr
    publisher: CDS

    ivoid: ivo://cds.vizier/j/aj/167/18
    citation: doi:10.26093/cds/vizier.51670018
    reference_url: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/167/18
    rights_uri: https://cds.unistra.fr/vizier-org/licences_vizier.html
    creator: Hong K.
    ...

Contents and metadata
---------------------

`astropy.io.votable.dataorigin.extract_data_origin` returns a `astropy.io.votable.dataorigin.DataOrigin` (class) container which is made of:

* a `astropy.io.votable.dataorigin.QueryOrigin` (class) container describing the request.
  ``QueryOrigin`` is considered to be unique for the whole VOTable.
  It includes metadata like  the publisher, the contact, date of execution, query, etc.

*  a list of `astropy.io.votable.dataorigin.DatasetOrigin` (class) container for each Element having DataOrigin information.
   ``DataSetOrigin`` is a basic provenance of the datasets queried. Each attribute is a list.
   It includes metadata like authors, ivoid, landing pages, ....

Examples
--------

Get the (Data Center) publisher and the Creator of the dataset

    >>> print(data_origin.query.publisher)
    CDS
    >>> print(data_origin.origin[0].creator)
    ['Hong K.']

Other capabilities
------------------

DataOrigin container includes VO Elements:

* Extract list of `astropy.io.votable.tree.Info`


    >>> # get DataOrigin with the description of each INFO
    >>> for dataset_origin in data_origin.origin:
    ...    for info in dataset_origin.infos:
    ...        print(f"{info.name}: {info.value} ({info.content})")
    ivoid: ivo://cds.vizier/j/aj/167/18 (IVOID of underlying data collection)
    creator: Hong K. (First author or institution)
    cites: bibcode:2024AJ....167...18H (Article or Data origin sources)
    editor: Astronomical Journal (AAS) (Editor name (article))
    original_date: 2024 (Year of the article publication)
    ...

* Extract tree node `astropy.io.votable.tree.Element`

The following example extracts the citation from the header (in APA style).

    >>> # get the Title retrieved in Element
    >>> origin = data_origin.origin[0]
    >>> vo_elt = origin.get_votable_element()
    >>> title = vo_elt.description if vo_elt else ""
    >>> print(f"APA: {','.join(origin.creator)} ({origin.publication_date[0]}). {title} [Dataset]. {data_origin.query.publisher}. {origin.citation[0]}")
    APA: Hong K. (2024-11-06). Period variations of 32 contact binaries (Hong+, 2024) [Dataset]. CDS. doi:10.26093/cds/vizier.51670018

* Add Data Origin INFO into VOTable:

.. doctest-skip::

    >>> votable = parse("votable.xml")
    >>> dataorigin.add_data_origin_info(votable, "query", "Data center name")
    >>> dataorigin.add_data_origin_info(votable.resources[0], "creator", "Author name")
