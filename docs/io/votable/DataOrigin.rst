.. _astropy-io-votable:

*******************************************
DataOrigin (`astropy.io.votable`)
*******************************************

Introduction
============
DataOrigin is described in the IVOA note: https://www.ivoa.net/documents/DataOrigin/.

This API retrieves Metadata from INFO in VOTAble.


Getting Started
===============
To extract DataOrigin from VOTable
`~astropy.io.votable.DataOrigin`::

    from astropy.io.votable import parse
    from astropy.io.votable.DataOrigin import extract_data_origin

    votable = parse("votable.xml")
    data_origin = extract_data_origin(votable)
    print(data_origin)

Contents and metadata
=====================

```extract_data_origin``` returns a ```DataOrigin``` (class) container which is made of:

* a ```QueryOrigin``` (class) container describing the request.
      QueryOrigin is considered to be unique for the whole VOTable.
      It includes metadata like  the publisher, the contact, date of execution, query, etc.

*  a list of ```DatasetOrigin``` (class) container for each Element having DataOrigin information.
      DataSetOrigin is a basic provenance of the datasets queried. Each attribute is a list.
      It includes metadata like authors, ivoid, landing pages, ....

Example: Get the (Data Center) publisher and the Creator of the dataset

`~astropy.io.votable.DataOrigin`::

    from astropy.io.votable import parse
    from astropy.io.votable.DataOrigin import extract_data_origin

    votable = parse("votable.xml")
    data_origin = extract_data_origin(votable)

    uri_request = data_origin.query.publisher
    creators =  data_origin.origin[0].creator

Other capabilities
==================
DataOrigin container includes VO Elements:

* DataOrigin.info : list of astropy.io.votable.tree.Info matched with DataOrigin
`~astropy.io.votable.DataOrigin`::

    data_origin = extract_data_origin(votable)

    # get DataOrigin with the description of each INFO
    for dataset_origin in data_origin.origin:
      for info in dataset_origin.infos:
        print(f"{info.name}: {info.value} ({info.content})")

* DatasetOrigin.get_vo_element(): get the (tree) Element where are extracted DataOrigin
`~astropy.io.votable.DataOrigin`::

    data_origin = extract_data_origin(votable)

    # get the Title retrived in Element
    dataset_origin = data_origin.origin[0]
    vo_elt = dataset_origin.get_votable_element()
    if vo_elt:
        title = vo_elt.description

* Add INFO Data Origin into a VOTAble

`~astropy.io.votable.DataOrigin`::

        dataorigin.add_data_origin_info(votable, "publisher", "Data center name")
        dataorigin.add_data_origin_info(votable.resources[0], "creator", "Author name")

