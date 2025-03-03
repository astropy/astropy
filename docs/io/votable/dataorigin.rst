.. _astropy-io-votable-dataorigin_intro:

********************************************
DataOrigin (`astropy.io.votable.dataorigin`)
********************************************

Introduction
============
Extract basic provenance information from VOTable header. The information is described in 
DataOrigin IVOA note: https://www.ivoa.net/documents/DataOrigin/.

DataOrigin includes both the query information (such as publisher, contact, versions, etc.) 
and the Dataset origin (such as Creator, bibliographic links, URL, etc.)

This API retrieves Metadata from INFO in VOTAble.


Getting Started
===============

To extract DataOrigin from VOTable

.. code-block:: python

    from astropy.io.votable import parse
    from astropy.io.votable.DataOrigin import extract_data_origin

    votable = parse("https://vizier.cds.unistra.fr/conesearch/II/246/out?RA=0&DEC=0&SR=0.1")
    data_origin = extract_data_origin(votable)
    print(data_origin)

Contents and metadata
=====================

``extract_data_origin`` returns a ``DataOrigin`` (class) container which is made of:

* a ``QueryOrigin`` (class) container describing the request.
      QueryOrigin is considered to be unique for the whole VOTable.
      It includes metadata like  the publisher, the contact, date of execution, query, etc.

*  a list of ``DatasetOrigin`` (class) container for each Element having DataOrigin information.
      DataSetOrigin is a basic provenance of the datasets queried. Each attribute is a list.
      It includes metadata like authors, ivoid, landing pages, ....

Example: Get the (Data Center) publisher and the Creator of the dataset

.. code-block:: python

    from astropy.io.votable import parse
    from astropy.io.votable.DataOrigin import extract_data_origin

    votable = parse("https://vizier.cds.unistra.fr/conesearch/II/246/out?RA=0&DEC=0&SR=0.1")
    data_origin = extract_data_origin(votable)

    uri_request = data_origin.query.publisher
    creators =  data_origin.origin[0].creator

Other capabilities
==================
DataOrigin container includes VO Elements:

* DataOrigin ``info`` : list of astropy.io.votable.tree.Info matched with DataOrigin

.. code-block:: python

    from astropy.io.votable import parse
    from astropy.io.votable.dataorigin import extract_data_origin

    votable = parse("https://vizier.cds.unistra.fr/viz-bin/conesearch/II/246/out?RA=0&DEC=0&SR=0.1")
    data_origin = extract_data_origin(votable)

    # get DataOrigin with the description of each INFO
    for dataset_origin in data_origin.origin:
        for info in dataset_origin.infos:
            print(f"{info.name}: {info.value} ({info.content})")

* DatasetOrigin ``get_vo_element()``: get the (tree) Element where are extracted DataOrigin

.. code-block:: python

    data_origin = extract_data_origin(votable)

    # get the Title retrieved in Element
    dataset_origin = data_origin.origin[0]
    vo_elt = dataset_origin.get_votable_element()
    if vo_elt:
        title = vo_elt.description

* Add INFO Data Origin into a VOTable

.. code-block:: python

        votable = parse("votable.xml")
        dataorigin.add_data_origin_info(votable, "publisher", "Data center name")
        dataorigin.add_data_origin_info(votable.resources[0], "creator", "Author name")

