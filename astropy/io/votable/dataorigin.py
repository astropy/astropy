# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Extract Data Origin in VOTable

References
----------
DataOrigin is a vocabulary described in the IVOA note: https://www.ivoa.net/documents/DataOrigin/

Notes
-----
This API retrieve Metadata from INFO in VOTable.
The information can be found at different level in a VOTable :

- global level
- resource level
- table level

Contents
--------

- Query information: Each element is considered to be unique in the VOTable.
  Information concerns publisher, date of execution, contact, request, etc.
- Dataset origin : basic provenance information.

Examples
--------
>>> from astropy.io.votable import parse
>>> from astropy.io.votable.dataorigin import extract_data_origin
>>> data_origin = extract_data_origin(parse("https://vizier.cds.unistra.fr/viz-bin/conesearch/II/246/out?RA=0&DEC=0&SR=0.1"))  # doctest: +REMOTE_DATA
>>> print(data_origin)  # doctest: +REMOTE_DATA +IGNORE_OUTPUT
publisher: CDS
server_software: 7.4.5
service_protocol: ivo://ivoa.net/std/ConeSearch/v1.03
request: https://vizier.cds.unistra.fr/viz-bin/conesearch/II/246/out?RA=0&DEC=0&SR=0.1
request_date: 2025-03-03T12:54:26
contact: cds-question@unistra.fr
...
>>> uri_request = data_origin.query.request  # doctest: +REMOTE_DATA
>>> creators =  data_origin.origin[0].creator  # doctest: +REMOTE_DATA
"""

import astropy.io.votable.tree

__all__ = [
    "DataOrigin",
    "DatasetOrigin",
    "QueryOrigin",
    "add_data_origin_info",
    "extract_data_origin",
]


DATAORIGIN_QUERY_INFO = (
    "ivoid_service",
    "publisher",
    "server_software",
    "service_protocol",
    "request",
    "query",
    "request_date",
    "contact",
)
DATAORIGIN_INFO = (
    "ivoid",
    "citation",
    "reference_url",
    "resource_version",
    "rights_uri",
    "rights",
    "creator",
    "editor",
    "article",
    "cites",
    "is_derived_from",
    "original_date",
    "publication_date",
    "last_update_date",
)


class QueryOrigin:
    """Data class storing query execution information that generated the VOTable.

    Notes
    -----
    The Query information should be unique in the whole VOTable.
    It includes reproducibility information to execute the query again.

    Attributes
    ----------
    ivoid_service : str
        IVOID of the service that produced the VOTable (default: None)

    publisher : str
        Data centre that produced the VOTable (default: None)

    server_software : str
        Software version (default: None)

    service_protocol : str
        IVOID of the protocol through which the data was retrieved (default: None)

    request : str
        Full request URL including a query string (default: None)

    query : str
        An input query in a formal language (e.g, ADQL) (default: None)

    request_date : str
        Query execution date (default: None)

    contact : str
        Email or URL to contact publisher (default: None)

    infos : list[astropy.io.votable.tree.Info]
        list of ``<INFO>`` used by DataOrigin (default: empty list)

    """

    def __init__(self):
        self.ivoid_service = None
        self.publisher = None
        self.server_software = None
        self.service_protocol = None
        self.request = None
        self.query = None
        self.request_date = None
        self.contact = None
        self.infos = []

    def __str__(self) -> str:
        s = []
        for info_name in DATAORIGIN_QUERY_INFO:
            info = getattr(self, info_name)
            if info:
                s.append(f"{info_name}: {info}")
        return "\n".join(s)


class DatasetOrigin:
    """Data class storing the basic provenance for a Dataset.

    Notes
    -----
    DatasetOrigin is dedicated to a specific Element in a VOTable.
    These ``<INFO>`` Elements describe a Resource, a TableElement or are Global.

    Attributes
    ----------
    ivoid : list
        IVOID of underlying data collection (default: None)

    citation : list
        Dataset identifier that can be used for citation (default: None)

    reference_url : list
        Dataset landing page (default: None)

    resource_version : list
        Dataset version (default: None)

    rights_uri : list
        Licence URI (default: None)

    rights : list
        Licence or Copyright text (default: None)

    creator : list
        The person(s) mainly involved in the creation of the resource (default: None)

    editor : list
        Editor name of the reference article (default: None)

    article : list
        Bibcode or DOI of a reference article (default: None)

    cites : list
        An Identifier (ivoid, DOI, bibcode) of second resource (default: None)

    is_derived_from : list
        An Identifier (ivoid, DOI, bibcode) of second resource (default: None)

    original_date : list
        Date of the original resource from which the present resource is derived (default: None)

    publication_date : list
        Date of first publication in the data centre (default: None)

    last_update_date : list
        Last data centre update (default: None)

    infos : list[astropy.io.votable.tree.Info]
        list of ``<INFO>`` used by DataOrigin (default: None)
    """

    def __init__(self, votable_element: astropy.io.votable.tree.Element = None):
        """
        Constructor

        Parameters
        ----------
        votable_element: astropy.io.votable.tree.Element, optional
                         indicates the VOTable element
        """
        self.ivoid = None
        self.citation = None
        self.reference_url = None
        self.resource_version = None
        self.rights_uri = None
        self.rights = None
        self.creator = None
        self.editor = None
        self.article = None
        self.cites = None
        self.is_derived_from = None
        self.original_date = None
        self.publication_date = None
        self.last_update_date = None
        self.__vo_elt = votable_element
        self.infos = []

    def get_votable_element(self) -> astropy.io.votable.tree.Element:
        """
        Get the VOTable element

        Returns
        -------
        astropy.io.votable.tree.Element
        """
        return self.__vo_elt

    def __str__(self) -> str:
        s = []
        for info_name in DATAORIGIN_INFO:
            info = getattr(self, info_name)
            if info:
                s.append(f"{info_name}: {','.join(info)}")
        return "\n".join(s)


class DataOrigin:
    """Data class storing both information about query execution
       and basic provenances of datasets used to generate the VOTable.

    Attributes
    ----------
    query : QueryOrigin
        request information (default: None)

    origin : list[DatasetOrigin]
        list of DatasetOrigin (default: empty)

    Notes
    -----
    The class includes an iterator on Attribute origin.
    """

    def __init__(self):
        self.query = QueryOrigin()
        self.origin = []
        self.__it = None

    def __str__(self) -> str:
        origin_list = []
        for origin in self.origin:
            origin_list.append(str(origin))
        return str(self.query) + "\n\n" + "\n\n".join(origin_list)

    def __iter__(self):
        self.__it = -1
        return self

    def __next__(self):
        self.__it += 1
        if self.__it >= len(self.origin):
            raise StopIteration
        return self.origin[self.__it]


def __empty_dataset_origin(o: DatasetOrigin) -> bool:
    """(internal) check if DataOrigin is filled"""
    for info in DATAORIGIN_INFO:
        v = getattr(o, info)
        if v is not None:
            return False
    return True


def __extract_generic_info(
    vo_element: astropy.io.votable.tree.Element, infos: list, data_origin: DataOrigin
):
    """(internal) extract info and populate DataOrigin

    Parameters
    ----------
    vo_element : astropy.io.votable.tree.Element
        VOTable element (votable, resource or table)

    infos : list[astropy.io.votable.tree.Info]
        list of ``<INFO>``

    data_origin : DataOrigin
        DataOrigin container to fill
    """
    if not infos:
        return

    dataset_origin = DatasetOrigin(vo_element)

    for info in infos:
        info_name = info.name.lower()
        for dataorigin_info in DATAORIGIN_INFO:
            if info_name == dataorigin_info:
                dataset_origin.infos.append(info)
                att = getattr(dataset_origin, dataorigin_info)
                if att is None:
                    setattr(dataset_origin, dataorigin_info, [info.value])
                else:
                    att.append(info.value)
                break

        for query_info in DATAORIGIN_QUERY_INFO:
            if info_name == query_info:
                data_origin.query.infos.append(info)
                setattr(data_origin.query, query_info, info.value)
                break

    if not __empty_dataset_origin(dataset_origin):
        data_origin.origin.append(dataset_origin)


def __extract_dali_info(infos: list, data_origin: DataOrigin):
    """(internal) append with DALI INFO

    Parameters
    ----------
    infos : list[astropy.io.votable.tree.Info]
        iterable info

    data_origin DataOrigin
        container to fill
    """
    if not data_origin.query.service_protocol:
        for info in infos:
            info_name = info.name.lower()
            if info_name == "standardid":
                if not data_origin.query.service_protocol:
                    if data_origin.info is None:
                        data_origin.infos = []
                    data_origin.quey.infos.append(info)
                    data_origin.query.service_protocol = info.value


def __extract_info_from_table(
    table: astropy.io.votable.tree.TableElement, data_origin: DataOrigin
):
    """(internal) extract and populate dataOrigin from astropy.io.votable.tree.TableElement

    Parameters
    ----------
    table : astropy.io.votable.tree.TableElement
        Table to explore

    data_origin : DataOrigin
        container to fill.
    """
    __extract_generic_info(table, table.infos, data_origin)


def __extract_info_from_resource(
    resource: astropy.io.votable.tree.Resource,
    data_origin: DataOrigin,
    recursive: bool = True,
):
    """(internal) extract and populate dataOrigin from astropy.io.votable.tree.Resource

    Parameters
    ----------
    param resource : astropy.io.votable.tree.Resource
        Resource to explore.

    data_origin : DataOrigin
        container to fill.

    recursive : bool, optional
        make a recursive search (default: True)
    """
    __extract_generic_info(resource, resource.infos, data_origin)
    __extract_dali_info(resource.infos, data_origin)
    if recursive:
        for table in resource.tables:
            __extract_info_from_table(table, data_origin)


def __extract_info_from_votable(
    votable: astropy.io.votable.tree.VOTableFile,
    data_origin: DataOrigin,
    recursive: bool = True,
):
    """(internal) extract and populate dataOrigin from astropy.io.votable.tree.VOTableFile

    Parameters
    ----------
    votable : astropy.io.votable.tree.VOTableFile
        VOTableFile to explore

    data_origin : DataOrigin
        container to fill

    recursive : bool, optional
        make a recursive search (default: True)
    """
    __extract_generic_info(votable, votable.infos, data_origin)
    if recursive:
        for resource in votable.resources:
            __extract_info_from_resource(resource, data_origin)


def extract_data_origin(vot_element: astropy.io.votable.tree.Element) -> DataOrigin:
    """Extract DataOrigin in a VO element

    Parameters
    ----------
    vot_element : astropy.io.votable.tree.Info
        VOTable Element to explore

    Returns
    -------
    DataOrigin

    Raises
    ------
    TypeError
        input ``vot_element`` type is not supported
    """
    data_origin = DataOrigin()
    if isinstance(vot_element, astropy.io.votable.tree.VOTableFile):
        __extract_info_from_votable(vot_element, data_origin)
    elif isinstance(vot_element, astropy.io.votable.tree.Resource):
        __extract_info_from_resource(vot_element, data_origin)
    elif isinstance(vot_element, astropy.io.votable.tree.TableElement):
        __extract_info_from_table(vot_element, data_origin)
    else:
        raise TypeError("input vot_element type is not supported.")

    return data_origin


def add_data_origin_info(
    vot_element: astropy.io.votable.tree.Element,
    info_name: str,
    info_value: str,
    content: str | None = None,
) -> None:
    """Update VOTable element with information compatible
       with DataOrigin vocabulary.

    Notes
    -----
    The function checks information name and adds the
    VOTable element with a new ``<INFO>``.

    Parameters
    ----------
    vot_element : astropy.io.votable.tree.Element
        VOTable element where to add the information

    info_name : str
        Attribute name (see DATAORIGIN_INFO, DATAORIGIN_QUERY_INFO)

    info_value : str
        value

    content : str, optional
        Content in ``<INFO>`` (default: None)

    Raises
    ------
    TypeError
        input type not managed or information name not recognized
    ValueError
        ``info_name`` already exists in ``vot_element``
    ValueError
        ``info_name`` is an unknown DataOrigin name.
    """
    if info_name in DATAORIGIN_INFO:
        if not isinstance(
            vot_element,
            (
                astropy.io.votable.tree.VOTableFile,
                astropy.io.votable.tree.Resource,
                astropy.io.votable.tree.TableElement,
            ),
        ):
            raise TypeError("Unsupported vot_element type.")

        vot_element.infos.extend(
            [astropy.io.votable.tree.Info(name=info_name, value=info_value)]
        )
        return

    elif info_name in DATAORIGIN_QUERY_INFO:
        if not isinstance(vot_element, astropy.io.votable.tree.VOTableFile):
            raise TypeError(
                "Bad type of vot_element: this information needs VOTableFile."
            )

        for info in vot_element.get_infos_by_name(info_name):
            raise ValueError(f"QueryOrigin {info_name} already exists")
        new_info = astropy.io.votable.tree.Info(name=info_name, value=info_value)
        if content:
            new_info.content = content
        vot_element.infos.extend([new_info])
        return

    raise ValueError("Unknown DataOrigin info name.")
