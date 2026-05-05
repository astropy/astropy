# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Extract Data Origin in VOTable.

References
----------
DataOrigin is a vocabulary described in the IVOA note: https://www.ivoa.net/documents/DataOrigin/

Notes
-----
This API retrieve Metadata from INFO in VOTable.
The information can be found at different level in a VOTable:

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
For more information, please see :ref:`DataOrigin documentation <astropy-io-votable-dataorigin>`.
"""

import astropy.io.votable.tree
from astropy.utils.decorators import deprecated_attribute

__all__ = [
    "DataOrigin",
    "DatasetOrigin",
    "QueryOrigin",
    "add_data_origin_info",
    "extract_data_origin",
]


DATAORIGIN_QUERY_INFO = (
    "service_ivoid",
    "publisher",
    "server_software",
    "service_protocol",
    "request",
    "query",
    "request_date",
    "contact",
)


DATAORIGIN_INFO = (
    "data_ivoid",
    "citation",
    "reference_url",
    "resource_version",
    "rights_uri",
    "rights",
    "creator",
    "journal",
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
    service_ivoid : str
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

    _INFO_MAPPING = ("standardid",)  # DALI INFO

    def __init__(self, votable_element: astropy.io.votable.tree.Element = None):
        self.service_ivoid = None
        self.publisher = None
        self.server_software = None
        self.service_protocol = None
        self.request = None
        self.query = None
        self.request_date = None
        self.contact = None
        self.infos = []

    @property
    def standardID(self) -> list:
        """Compatibility with previous version"""
        return self.service_protocol

    @standardID.setter
    def standardID(self, value: list):
        """Compatibility with previous version"""
        if not self.service_protocol:
            self.service_protocol = value

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
    data_ivoid : list
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

    journal : list
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

    _INFO_MAPPING = ("editor", "ivoid")  # obsolete INFO

    def __init__(self, votable_element: astropy.io.votable.tree.Element = None):
        """
        Constructor

        Parameters
        ----------
        votable_element: astropy.io.votable.tree.Element, optional
                         indicates the VOTable element
        """
        self.data_ivoid = None
        self.citation = None
        self.reference_url = None
        self.resource_version = None
        self.rights_uri = None
        self.rights = None
        self.creator = None
        self.journal = None
        self.article = None
        self.cites = None
        self.is_derived_from = None
        self.original_date = None
        self.publication_date = None
        self.last_update_date = None
        self.__vo_elt = votable_element
        self.infos = []
        self.ivoid = deprecated_attribute(
            name="ivoid", alternative="data_ivoid", since="8.0"
        )
        self.editor = deprecated_attribute(
            name="editor", alternative="journal", since="8.0"
        )

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

    def is_empty(self) -> bool:
        """check if DataOrigin is filled

        Returns
        -------
        bool
        """
        for info in DATAORIGIN_INFO:
            v = getattr(self, info)
            if v is not None:
                return False
        return True


class DataOrigin:
    """Class parsing a VOTable and storing both information about query execution
    and basic provenances.
    """

    def __init__(self, vot_element=None):
        self.query = QueryOrigin()
        self.origin = []
        self.__it = None

        self.__vot_element = vot_element
        if vot_element:
            self.parse()

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

    def __extract_generic_info(
        self, vo_fragment: astropy.io.votable.tree.Element, infos: list
    ):
        """(internal) extract info and populate DataOrigin

        Parameters
        ----------
        vo_fragment : astropy.io.votable.tree.Element
            VOTable element (votable, resource or table)

        infos : list[astropy.io.votable.tree.Info]
            list of ``<INFO>``

        """
        if not infos:
            return

        dataset_origin = DatasetOrigin(vo_fragment)

        dataset_origin_terms = DATAORIGIN_INFO + DatasetOrigin._INFO_MAPPING
        query_origin_terms = DATAORIGIN_QUERY_INFO + QueryOrigin._INFO_MAPPING

        for info in infos:
            info_name = info.name.lower()
            for dataset_info in dataset_origin_terms:
                if info_name == dataset_info:
                    dataset_origin.infos.append(info)
                    att = getattr(dataset_origin, dataset_info)
                    if att is None or isinstance(att, property):
                        setattr(dataset_origin, dataset_info, [info.value])
                    else:
                        att.append(info.value)
                    break

            for query_info in query_origin_terms:
                if info_name == query_info:
                    self.query.infos.append(info)
                    setattr(self.query, query_info, info.value)
                    break

        if not dataset_origin.is_empty():
            self.origin.append(dataset_origin)

    def __extract_info_from_table(self, table: astropy.io.votable.tree.TableElement):
        """(internal) extract and populate dataOrigin from astropy.io.votable.tree.TableElement

        Parameters
        ----------
        table : astropy.io.votable.tree.TableElement
            Table to explore.
        """
        self.__extract_generic_info(table, table.infos)

    def __extract_info_from_resource(
        self,
        resource: astropy.io.votable.tree.Resource,
        recursive: bool = True,
    ):
        """(internal) extract and populate dataOrigin from astropy.io.votable.tree.Resource

        Parameters
        ----------
        resource : astropy.io.votable.tree.Resource
            Resource to explore.

        recursive : bool, optional
            make a recursive search (default: True)
        """
        self.__extract_generic_info(resource, resource.infos)
        if recursive:
            for table in resource.tables:
                self.__extract_info_from_table(table)

    def __extract_info_from_votable(
        self,
        votable: astropy.io.votable.tree.VOTableFile,
        recursive: bool = True,
    ):
        """(internal) extract and populate dataOrigin from astropy.io.votable.tree.VOTableFile

        Parameters
        ----------
        votable : astropy.io.votable.tree.VOTableFile
            VOTableFile to explore.

        recursive : bool, optional
            make a recursive search (default: True)
        """
        self.__extract_generic_info(votable, votable.infos)
        if recursive:
            for resource in votable.resources:
                self.__extract_info_from_resource(resource)

    def parse(self) -> None:
        """Extract DataOrigin in a VO element

        Raises
        ------
        TypeError
            input ``vot_element`` type is not supported
        """
        if isinstance(self.__vot_element, astropy.io.votable.tree.VOTableFile):
            self.__extract_info_from_votable(self.__vot_element)
        elif isinstance(self.__vot_element, astropy.io.votable.tree.Resource):
            self.__extract_info_from_resource(self.__vot_element)
        elif isinstance(self.__vot_element, astropy.io.votable.tree.TableElement):
            self.__extract_info_from_table(self.__vot_element)
        else:
            raise TypeError("input vot_element type is not supported.")

    @staticmethod
    def __clean_votable_info(vot_element: astropy.io.votable.tree.Element) -> None:
        """(internal) Clean existing DataOrigin INFO in the VOTable Element

        Parameters
        ----------
        vot_element : astropy.io.votable.tree.Element
            VOTable Element where to remove the INFO

        """
        for info in vot_element.infos[0:]:
            if (
                info.name in DATAORIGIN_QUERY_INFO
                or info.name in DATAORIGIN_INFO
                or info.name in QueryOrigin._INFO_MAPPING
                or info.name in DatasetOrigin._INFO_MAPPING
            ):
                vot_element.infos.remove(info)

        if isinstance(vot_element, astropy.io.votable.tree.Resource):
            for table in vot_element.resources:
                DataOrigin.__clean_votable_info(table)
        elif isinstance(vot_element, astropy.io.votable.tree.VOTableFile):
            for resource in vot_element.resources:
                DataOrigin.__clean_votable_info(resource)

    @staticmethod
    def __append_votable_info(
        vot_element: astropy.io.votable.tree.Element,
        name: str,
        value: str | list,
        content: str | None = None,
        unique: bool = False,
    ) -> None:
        """(internal) add new DATAOrigin info in the VOTable Element

        Parameters
        ----------
        vot_element : astropy.io.votable.tree.Element
            VOTable Element where to add a new INFO

        name : str
            INFO name

        value: str | list
            the INFO value

        content: str, optional
            INFO description (default: None)

        unique: bool, optional
            the INFO element is unique (default: False)

        """
        if not isinstance(
            vot_element,
            (
                astropy.io.votable.tree.VOTableFile,
                astropy.io.votable.tree.Resource,
                astropy.io.votable.tree.TableElement,
            ),
        ):
            raise TypeError("input vot_element type is not supported.")

        for info in vot_element.infos:
            if info.name == name:
                if unique:
                    return
                if info.value == value:
                    return

        values = [value] if not isinstance(value, list) else value
        for val in values:
            new_info = astropy.io.votable.tree.Info(name=name, value=val)
            if content:
                new_info.content = content
            vot_element.infos.extend([new_info])

    def update_votable(self):
        """Update the VOTable fragment with DataOrigin <INFO>

        Returns
        -------
        astropy.io.votable.tree.Element
        """
        if not self.__vot_element:
            raise ValueError("VOTable not parsed yet (please call parse method first)")

        # clean existing DataOrigin info
        DataOrigin.__clean_votable_info(self.__vot_element)

        for item in DATAORIGIN_QUERY_INFO:
            att = getattr(self.query, item)
            if not att:
                continue

            DataOrigin.__append_votable_info(
                self.__vot_element, name=item, value=att, unique=True
            )

        for origin_info in self.origin:
            for item in DATAORIGIN_INFO:
                att = getattr(origin_info, item)
                if not att:
                    continue

                vot_fragment = origin_info.get_votable_element()
                if not vot_fragment:
                    vot_fragment = self.__vo_elt
                DataOrigin.__append_votable_info(vot_fragment, name=item, value=att)

        return self.__vot_element


def extract_data_origin(vot_element: astropy.io.votable.tree.Element) -> DataOrigin:
    """Extract DataOrigin in a VO element
      (keep compatibility with previous version)

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
    return DataOrigin(vot_element)


def add_data_origin_info(
    vot_element: astropy.io.votable.tree.VOTableFile
    | astropy.io.votable.tree.Resource
    | astropy.io.votable.tree.TableElement,
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
    vot_element : astropy.io.votable.tree.VOTableFile | astropy.io.votable.tree.Resource | astropy.io.votable.tree.TableElement
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

    elif info_name in DATAORIGIN_QUERY_INFO:
        if not isinstance(vot_element, astropy.io.votable.tree.VOTableFile):
            raise TypeError(
                "Bad type of vot_element: this information needs VOTableFile."
            )

        for info in vot_element.get_infos_by_name(info_name):
            raise ValueError(f"QueryOrigin {info_name} already exists")

    else:
        raise ValueError("Unknown DataOrigin info name.")

    new_info = astropy.io.votable.tree.Info(name=info_name, value=info_value)
    if content:
        new_info.content = content
    vot_element.infos.extend([new_info])
