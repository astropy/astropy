"""
This file defines the nodes that make up the VOTABLE XML tree.

.. _BINARY: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC27
.. _COOSYS: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC19
.. _DESCRIPTION: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC19
.. _FIELD: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC24
.. _FIELDref: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _FITS: http://fits.gsfc.nasa.gov/fits_documentation.html
.. _GROUP: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _ID: http://www.w3.org/TR/REC-xml/#id
.. _INFO: http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC19
.. _LINK: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC22
.. _multidimensional arrays: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC12
.. _numerical accuracy: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC26
.. _PARAM: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC24
.. _PARAMref: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _RESOURCE: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC21
.. _TABLE: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC23
.. _TABLEDATA: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC25
.. _unified content descriptor: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC28
.. _unique type: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC29
.. _units: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC27
.. _VALUES: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC30
.. _VOTABLE: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC9
"""

# TODO: Test FITS parsing

from __future__ import division, absolute_import

from .util import IS_PY3K

# STDLIB
import base64
import codecs
import io
from math import ceil
from operator import attrgetter
import re
import sys
import urllib2
if IS_PY3K:
    basestring = (str, bytes)

# THIRD-PARTY
import numpy as np
try:
    import pyfits
    _has_pyfits = True
except ImportError:
    _has_pyfits = False

# LOCAL
from . import converters
from .voexceptions import warn_or_raise, vo_warn, vo_raise, \
    vo_reraise, warn_unknown_attrs, UnimplementedWarning, \
    VOTableChangeWarning, W06, W07, W08, W09, W10, W11, W12, W13, \
    W15, W17, W18, W19, W20, W21, W22, W26, W27, W28, W29, W32, W33, \
    W35, W36, W37, W38, W40, W41, W42, W43, W44, W45, W48, E08, E09, \
    E10, E11, E12, E13, E14, E15, E16, E17, E18, E19, E20, E21
from . import ucd as ucd_mod
from . import util
from . import xmlutil
from astropy import __version__ as astropy_version

try:
    from . import iterparser
    _has_c_tabledata_writer = True
except ImportError:
    _has_c_tabledata_writer = False

# The default number of rows to read in each chunk before converting
# to an array.
DEFAULT_CHUNK_SIZE = 256
RESIZE_AMOUNT = 1.5

######################################################################
# FACTORY FUNCTIONS


def _lookup_by_id_factory(iterator, element_name, doc):
    """
    Creates a function useful for looking up an element by ID.

    Parameters
    ----------
    iterator : generator
        A generator that iterates over some arbitrary set of elements

    element_name : str
        The XML element name of the elements being iterated over (used
        for error messages only).

    doc : str
        A docstring to apply to the generated function.

    Returns
    -------
    factory : function
        A function that looks up an element by ID
    """
    def lookup_by_id(self, ref, before=None):
        """
        Given an XML id *ref*, finds the first element in the iterator
        with the attribute ID == *ref*.  If *before* is provided, will
        stop searching at the object *before*.  This is important,
        since "forward references" are not allowed in the VOTABLE
        format.
        """
        for element in getattr(self, iterator)():
            if element is before:
                if element.ID == ref:
                    vo_raise(
                        "%s references itself" % element_name,
                        element._config, element._pos, KeyError)
                break
            if element.ID == ref:
                return element
        raise KeyError(
            "No %s with ID '%s' found before the referencing %s" %
            (element_name, ref, element_name))

    lookup_by_id.__doc__ = doc
    return lookup_by_id


def _lookup_by_id_or_name_factory(iterator, element_name, doc):
    """
    Like `_lookup_by_id_factory`, but also looks in the "name" attribute.
    """
    def lookup_by_id(self, ref, before=None):
        """
        Given an key *ref*, finds the first element in the iterator
        with the attribute ID == *ref* or name == *ref*.  If *before*
        is provided, will stop searching at the object *before*.  This
        is important, since "forward references" are not allowed in
        the VOTABLE format.
        """
        for element in getattr(self, iterator)():
            if element is before:
                if ref in (element.ID, element.name):
                    vo_raise(
                        "%s references itself" % element_name,
                        element._config, element._pos, KeyError)
                break
            if ref in (element.ID, element.name):
                return element
        raise KeyError(
            "No %s with ID or name '%s' found before the referencing %s" %
            (element_name, ref, element_name))

    lookup_by_id.__doc__ = doc
    return lookup_by_id


######################################################################
# ATTRIBUTE CHECKERS
def check_astroyear(year, field, config={}, pos=None):
    """
    Raises a `~astropy.io.vo.exceptions.VOTableSpecError` if *year* is
    not a valid astronomical year as defined by the VOTABLE standard.

    Parameters
    ----------
    year : str
        An astronomical year string

    field : str
        The name of the field this year was found in (used for error
        message)

    config, pos : optional
        Information about the source of the value
    """
    if (year is not None and
        re.match(r"^[JB]?[0-9]+([.][0-9]*)?$", year) is None):
        warn_or_raise(W07, W07, (field, year), config, pos)
        return False
    return True


def check_string(string, attr_name, config={}, pos=None):
    """
    Raises a `~astropy.io.vo.voexceptions.VOTableSpecError` if
    *string* is not a string or Unicode string.

    Parameters
    ----------
    string : str
        An astronomical year string

    field : str
        The name of the field this year was found in (used for error
        message)

    config, pos : optional
        Information about the source of the value
    """
    if string is not None and not isinstance(string, basestring):
        warn_or_raise(W08, W08, attr_name, config, pos)
        return False
    return True


def resolve_id(ID, id, config={}, pos=None):
    if ID is None and id is not None:
        warn_or_raise(W09, W09, (), config, pos)
        return id
    return ID


def check_ucd(ucd, config={}, pos=None):
    """
    Warns or raises a `~astropy.io.vo.voexceptions.VOTableSpecError`
    if *ucd* is not a valid `unified content descriptor`_ string as
    defined by the VOTABLE standard.

    Parameters
    ----------
    ucd : str
        A UCD string.

    config, pos : optional
        Information about the source of the value
    """
    if config.get('version_1_1_or_later'):
        try:
            ucd_mod.parse_ucd(
                ucd,
                check_controlled_vocabulary=config.get(
                    'version_1_2_or_later', False),
                has_colon=config.get('version_1_2_or_later', False))
        except ValueError as e:
            # This weird construction is for Python 3 compatibility
            if config.get('pedantic'):
                vo_raise(W06, (ucd, unicode(e)), config, pos)
            else:
                vo_warn(W06, (ucd, unicode(e)), config, pos)
                return False
    return True


######################################################################
# ELEMENT CLASSES
class Element(object):
    """
    A base class for all classes that represent XML elements in the
    VOTABLE file.
    """
    def _add_unknown_tag(self, iterator, tag, data, config, pos):
        warn_or_raise(W10, W10, tag, config, pos)

    def _ignore_add(self, iterator, tag, data, config, pos):
        warn_unknown_attrs(tag, data.iterkeys(), config, pos)

    def _add_definitions(self, iterator, tag, data, config, pos):
        if config.get('version_1_1_or_later'):
            warn_or_raise(W22, W22, (), config, pos)
        warn_unknown_attrs(tag, data.iterkeys(), config, pos)


class SimpleElement(Element):
    """
    A base class for simple elements, such as FIELD, PARAM and INFO
    that don't require any special parsing or outputting machinery.
    """
    def __init__(self):
        Element.__init__(self)

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start and tag != self._element_name:
                self._add_unknown_tag(iterator, tag, data, config, pos)
            elif tag == self._element_name:
                break

        return self

    def to_xml(self, w, **kwargs):
        w.element(self._element_name,
                  attrib=xmlutil.object_attrs(self, self._attr_list))


class SimpleElementWithContent(SimpleElement):
    """
    A base class for simple elements, such as FIELD, PARAM and INFO
    that don't require any special parsing or outputting machinery.
    """
    def __init__(self):
        SimpleElement.__init__(self)

        self._content = None

    def to_xml(self, w, **kwargs):
        w.element(self._element_name, self._content,
                  attrib=xmlutil.object_attrs(self, self._attr_list))

    def _set_content(self, content):
        check_string(content, 'content', self._config, self._pos)
        self._content = content
    def _del_content(self):
        self._content = None
    content = property(
        attrgetter('_content'), _set_content, _del_content,
        """
        The content of the element.
        """)


class Link(SimpleElement):
    """
    A class for storing LINK_ elements, which are used to reference
    external documents and servers through a URI.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _attr_list = ['ID', 'content_role', 'content_type', 'title', 'value',
                  'href', 'action']
    _element_name = 'LINK'

    def __init__(self, ID=None, title=None, value=None, href=None, action=None,
                 id=None, config={}, pos=None, **kwargs):
        self._config = config
        self._pos = pos

        SimpleElement.__init__(self)

        content_role = kwargs.get('content-role') or kwargs.get('content_role')
        content_type = kwargs.get('content-type') or kwargs.get('content_type')

        if 'gref' in kwargs:
            warn_or_raise(W11, W11, (), config, pos)

        self.ID           = resolve_id(ID, id, config, pos)
        self.content_role = content_role
        self.content_type = content_type
        self.title        = title
        self.value        = value
        self.href         = href
        self.action       = action

        warn_unknown_attrs(
            'LINK', kwargs.iterkeys(), config, pos,
            ['content-role', 'content_role', 'content-type', 'content_type',
             'gref'])

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID_ of the LINK_ element.
        """)

    def _set_content_role(self, content_role):
        if content_role not in (None, 'query', 'hints', 'doc', 'location'):
            vo_warn(W45, (content_role,), self._config, self._pos)
        self._content_role = content_role
    def _del_content_role(self):
        self._content_role = None
    content_role = property(
        attrgetter('_content_role'), _set_content_role, _del_content_role,
        """
        Defines the MIME role of the referenced object.  Must be one of:

          None, 'query', 'hints', 'doc' or 'location'
        """)

    def _set_content_type(self, content_type):
        xmlutil.check_mime_content_type(content_type, self._config, self._pos)
        self._content_type = content_type
    def _del_content_type(self):
        self._content_type = None
    content_type = property(
        attrgetter('_content_type'), _set_content_type, _del_content_type,
        """
        Defines the MIME content type of the referenced object.
        """)

    def _set_href(self, href):
        xmlutil.check_anyuri(href, self._config, self._pos)
        self._href = href
    def _del_href(self):
        self._href = None
    href = property(
        attrgetter('_href'), _set_href, _del_href,
        """
        A URI to an arbitrary protocol.  The vo package only
        supports http and anonymous ftp.
        """)


class Info(SimpleElementWithContent):
    """
    A class for storing INFO elements, which contain arbitrary
    key-value pairs for extensions to the standard.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _element_name = 'INFO'
    _attr_list_11 = ['ID', 'name', 'value']
    _attr_list_12 = _attr_list_11 + ['xtype', 'ref', 'unit', 'ucd', 'utype']

    def __init__(self, ID=None, name=None, value=None, id=None, xtype=None,
                 ref=None, unit=None, ucd=None, utype=None,
                 config={}, pos=None, **extra):
        self._config = config
        self._pos = pos

        SimpleElementWithContent.__init__(self)

        self.ID      = (resolve_id(ID, id, config, pos) or
                        xmlutil.fix_id(name, config, pos))
        self.name    = name
        self.value   = value
        self.xtype   = xtype
        self.ref     = ref
        self.unit    = unit
        self.ucd     = ucd
        self.utype   = utype

        if config.get('version_1_2_or_later'):
            self._attr_list = self._attr_list_12
        else:
            self._attr_list = self._attr_list_11
            if xtype is not None:
                warn_unknown_attrs('INFO', ['xtype'], config, pos)
            if ref is not None:
                warn_unknown_attrs('INFO', ['ref'], config, pos)
            if unit is not None:
                warn_unknown_attrs('INFO', ['unit'], config, pos)
            if ucd is not None:
                warn_unknown_attrs('INFO', ['ucd'], config, pos)
            if utype is not None:
                warn_unknown_attrs('INFO', ['utype'], config, pos)

        warn_unknown_attrs('INFO', extra.iterkeys(), config, pos)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID_ of the INFO element.  Used for cross-referencing.
        """)

    def _set_name(self, name):
        if name is None:
            warn_or_raise(W35, W35, ('name'), self._config, self._pos)
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name
    name = property(
        attrgetter('_name'), _set_name, None,
        """
        [*required*] The key of the key-value pair.
        """)

    def _set_value(self, value):
        if value is None:
            warn_or_raise(W35, W35, ('value'), self._config, self._pos)
        check_string(value, 'value', self._config, self._pos)
        self._value = value
    value = property(
        attrgetter('_value'), _set_value, None,
        """
        The value of the key-value pair.  (Always stored as a
        string or unicode string).
        """)

    def _set_content(self, content):
        check_string(content, 'content', self._config, self._pos)
        self._content = content
    def _del_content(self):
        self._content = None
    content = property(
        attrgetter('_content'), _set_content, _del_content,
        """
        The content inside the INFO element.
        """)

    def _set_xtype(self, xtype):
        if xtype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('xtype', 'INFO', '1.2'), config, pos)
        check_string(xtype, 'xtype', self._config, self._pos)
        self._xtype = xtype
    def _del_xtype(self):
        self._xtype = None
    xtype = property(
        attrgetter('_xtype'), _set_xtype, _del_xtype,
        """
        Extended data type information.
        """)

    def _set_ref(self, ref):
        if ref is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('ref', 'INFO', '1.2'), config, pos)
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        # TODO: actually apply the reference
        # if ref is not None:
        #     try:
        #         other = self._votable.get_values_by_id(ref, before=self)
        #     except KeyError:
        #         vo_raise(
        #             "VALUES ref='%s', which has not already been defined." %
        #             self.ref, self._config, self._pos, KeyError)
        #     self.null = other.null
        #     self.type = other.type
        #     self.min = other.min
        #     self.min_inclusive = other.min_inclusive
        #     self.max = other.max
        #     self.max_inclusive = other.max_inclusive
        #     self._options[:] = other.options
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        Refer to another INFO_ element by ID_, defined previously in
        the document.
        """)

    def _set_unit(self, unit):
        # TODO: Validate unit more accurately
        if unit is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('unit', 'INFO', '1.2'), config, pos)
        xmlutil.check_token(unit, 'unit', self._config, self._pos)
        self._unit = unit
    def _del_unit(self):
        self._unit = None
    unit = property(
        attrgetter('_unit'), _set_unit, _del_unit,
        """
        A string specifying the units_ for the INFO_.
        """)

    def _set_utype(self, utype):
        if utype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('utype', 'INFO', '1.2'), config, pos)
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the INFO_.
        """
        )


class Values(Element):
    """
    A class to represent the VALUES_ element, used within FIELD_ and
    PARAM_ elements to define the domain of values.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, votable, field, ID=None, null=None, ref=None,
                 type="legal", id=None, config={}, pos=None, **extras):
        self._config  = config
        self._pos = pos

        Element.__init__(self)

        self._votable = votable
        self._field   = field
        self.ID       = resolve_id(ID, id, config, pos)
        self.null     = null
        self._ref     = ref
        self.type     = type

        self.min           = None
        self.max           = None
        self.min_inclusive = True
        self.max_inclusive = True
        self._options      = []

        warn_unknown_attrs('VALUES', extras.iterkeys(), config, pos)

    def _get_null(self):
        return self._null
    def _set_null(self, null):
        if null is not None and isinstance(null, basestring):
            try:
                null_val = self._field.converter.parse_scalar(
                    null, self._config, self._pos)[0]
            except:
                warn_or_raise(W36, W36, null, self._config, self._pos)
                null_val = self._field.converter.parse_scalar(
                    '0', self._config, self._pos)[0]
        else:
            null_val = null
        self._null = null_val
    def _del_null(self):
        self._null = None
    null = property(
        _get_null, _set_null, _del_null,
        """
        For integral datatypes, *null* is used to define the value
        used for missing values.
        """)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the VALUES_ element, used for cross-referencing.
        May be ``None`` or a string conforming to XML ID_ syntax.
        """)

    def _set_type(self, type):
        if type not in ('legal', 'actual'):
            vo_raise(E08, type, self._config, self._pos)
        self._type = type
    type = property(
        attrgetter('_type'), _set_type, None,
        """
        [*required*] Defines the applicability of the domain defined
        by this VALUES_ element.  Must be one of the following
        strings:

          - 'legal': The domain of this column applies in general to
            this datatype. (default)

          - 'actual': The domain of this column applies only to the
            data enclosed in the parent table.
        """)

    def _set_ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        if ref is not None:
            try:
                other = self._votable.get_values_by_id(ref, before=self)
            except KeyError:
                warn_or_raise(W43, W43, ('VALUES', self.ref), self._config,
                              self._pos)
                ref = None
            else:
                self.null = other.null
                self.type = other.type
                self.min = other.min
                self.min_inclusive = other.min_inclusive
                self.max = other.max
                self.max_inclusive = other.max_inclusive
                self._options[:] = other.options
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        Refer to another VALUES_ element by ID_, defined previously in
        the document, for MIN/MAX/OPTION information.
        """)

    def _set_min(self, min):
        if hasattr(self._field, 'converter') and min is not None:
            self._min = self._field.converter.parse(min)[0]
        else:
            self._min = min
    def _del_min(self):
        self._min = None
    min = property(
        attrgetter('_min'), _set_min, _del_min,
        """
        The minimum value of the domain.  See :attr:`min_inclusive`.
        """)

    def _set_min_inclusive(self, inclusive):
        if inclusive == 'yes':
            self._min_inclusive = True
        elif inclusive == 'no':
            self._min_inclusive = False
        else:
            self._min_inclusive = bool(inclusive)
    def _del_min_inclusive(self):
        self._min_inclusive = True
    min_inclusive = property(
        attrgetter('_min_inclusive'), _set_min_inclusive, _del_min_inclusive,
        """
        When True, the domain includes the minimum value.
        """)

    def _set_max(self, max):
        if hasattr(self._field, 'converter') and max is not None:
            self._max = self._field.converter.parse(max)[0]
        else:
            self._max = max
    def _del_max(self):
        self._max = None
    max = property(
        attrgetter('_max'), _set_max, _del_max,
        """
        The maximum value of the domain.  See :attr:`max_inclusive`.
        """)

    def _set_max_inclusive(self, inclusive):
        if inclusive == 'yes':
            self._max_inclusive = True
        elif inclusive == 'no':
            self._max_inclusive = False
        else:
            self._max_inclusive = bool(inclusive)
    def _del_max_inclusive(self):
        self._max_inclusive = True
    max_inclusive = property(
        attrgetter('_max_inclusive'), _set_max_inclusive, _del_max_inclusive,
        """
        When True, the domain includes the maximum value.
        """)

    options = property(
        attrgetter('_options'), None, None,
        """
        A list of string key-value tuples defining other OPTION
        elements for the domain.  All options are ignored -- they are
        stored for round-tripping purposes only.
        """)

    def parse(self, iterator, config):
        if self.ref is not None:
            for start, tag, data, pos in iterator:
                if start:
                    warn_or_raise(W44, W44, tag, config, pos)
                else:
                    if tag != 'VALUES':
                        warn_or_raise(W44, W44, tag, config, pos)
                    break
        else:
            for start, tag, data, pos in iterator:
                if start:
                    if tag == 'MIN':
                        if 'value' not in data:
                            vo_raise(E09, 'MIN', config, pos)
                        self.min = data['value']
                        self.min_inclusive = data.get('inclusive', 'yes')
                        warn_unknown_attrs(
                            'MIN', data.iterkeys(), config, pos,
                            ['value', 'inclusive'])
                    elif tag == 'MAX':
                        if 'value' not in data:
                            vo_raise(E09, 'MAX', config, pos)
                        self.max = data['value']
                        self.max_inclusive = data.get('inclusive', 'yes')
                        warn_unknown_attrs(
                            'MAX', data.iterkeys(), config, pos,
                            ['value', 'inclusive'])
                    elif tag == 'OPTION':
                        if 'value' not in data:
                            vo_raise(E09, 'OPTION', config, pos)
                        xmlutil.check_token(
                            data.get('name'), 'name', config, pos)
                        self.options.append(
                            (data.get('name'), data.get('value')))
                        warn_unknown_attrs(
                            'OPTION', data.iterkeys(), config, pos,
                            ['data', 'name'])
                elif tag == 'VALUES':
                    break

        return self

    def is_defaults(self):
        # If there's nothing meaningful or non-default to write,
        # don't write anything.
        return (self.ref is None and self.null is None and self.ID is None and
                self.max is None and self.min is None and self.options == [])

    def to_xml(self, w, **kwargs):
        def yes_no(value):
            if value:
                return 'yes'
            return 'no'

        if self.is_defaults():
            return

        if self.ref is not None:
            w.element('VALUES', attrib=xmlutil.object_attrs(self, ['ref']))
        else:
            with w.tag('VALUES',
                       attrib=xmlutil.object_attrs(
                           self, ['ID', 'null', 'ref'])):
                if self.min is not None:
                    w.element(
                        'MIN',
                        value=self._field.converter.output(self.min, False),
                        inclusive=yes_no(self.min_inclusive))
                if self.max is not None:
                    w.element(
                        'MAX',
                        value=self._field.converter.output(self.max, False),
                        inclusive=yes_no(self.max_inclusive))
                for name, value in self.options:
                    w.element(
                        'OPTION',
                        name=name,
                        value=value)


class Field(SimpleElement):
    """
    A class that represents the FIELD_ element, which describes the
    datatype of a particular column of data.

    The keyword arguments correspond to setting members of the same
    name, documented below.

    If *ID* is provided, it is used for the column name in the
    resulting recarray of the table.  If no *ID* is provided, *name*
    is used instead.  If neither is provided, an exception will be
    raised.
    """
    _attr_list_11 = ['ID', 'name', 'datatype', 'arraysize', 'ucd',
                     'unit', 'width', 'precision', 'utype', 'ref']
    _attr_list_12 = _attr_list_11 + ['xtype']
    _element_name = 'FIELD'

    def __init__(self, votable, ID=None, name=None, datatype=None,
                 arraysize=None, ucd=None, unit=None, width=None,
                 precision=None, utype=None, ref=None, type=None, id=None,
                 xtype=None,
                 config={}, pos=None, **extra):
        self._config = config
        self._pos = pos

        SimpleElement.__init__(self)

        if config.get('version_1_2_or_later'):
            self._attr_list = self._attr_list_12
        else:
            self._attr_list = self._attr_list_11
            if xtype is not None:
                warn_unknown_attrs(self._element_name, ['xtype'], config, pos)

        # TODO: REMOVE ME ----------------------------------------
        # This is a terrible hack to support Simple Image Access
        # Protocol results from archive.noao.edu.  It creates a field
        # for the coordinate projection type of type "double", which
        # actually contains character data.  We have to hack the field
        # to store character data, or we can't read it in.  A warning
        # will be raised when this happens.
        if (not config.get('pedantic') and name == 'cprojection' and
            ID == 'cprojection' and ucd == 'VOX:WCS_CoordProjection' and
            datatype == 'double'):
            datatype = 'char'
            arraysize = '3'
            vo_warn(W40, (), config, pos)
        # ----------------------------------------

        self.description = None
        self._votable = votable

        self.ID = (resolve_id(ID, id, config, pos) or
                   xmlutil.fix_id(name, config, pos))
        self.name = name
        if name is None:
            if (self._element_name == 'PARAM' and
                not config.get('version_1_1_or_later')):
                pass
            else:
                warn_or_raise(W15, W15, self._element_name, config, pos)
            self.name = self.ID

        if self._ID is None and name is None:
            vo_raise(W12, self._element_name, config, pos)

        datatype_mapping = {
            'string'        : 'char',
            'unicodeString' : 'unicodeChar',
            'int16'         : 'short',
            'int32'         : 'int',
            'int64'         : 'long',
            'float32'       : 'float',
            'float64'       : 'double'}

        if datatype in datatype_mapping:
            warn_or_raise(W13, W13, (datatype, datatype_mapping[datatype]),
                          config, pos)
            datatype = datatype_mapping[datatype]

        self.ref        = ref
        self.datatype   = datatype
        self.arraysize  = arraysize
        self.ucd        = ucd
        self.unit       = unit
        self.width      = width
        self.precision  = precision
        self.utype      = utype
        self.type       = type
        self._links     = util.HomogeneousList(Link)
        self.title      = self.name
        self.values     = Values(self._votable, self)
        self.xtype      = xtype

        self._setup(config, pos)

        warn_unknown_attrs(self._element_name, extra.iterkeys(), config, pos)

    @classmethod
    def uniqify_names(cls, fields):
        """
        Make sure that all names and titles in a list of fields are
        unique, by appending numbers if necessary.
        """
        unique = {}
        for field in fields:
            i = 2
            new_id = field.ID
            while new_id in unique:
                new_id = field.ID + "_%d" % i
                i += 1
            if new_id != field.ID:
                vo_warn(W32, (field.ID, new_id), field._config, field._pos)
            field.ID = new_id
            unique[new_id] = field.ID

        for field in fields:
            i = 2
            if field.name is None or field.ID == xmlutil.fix_id(
                field.name, field._config, field._pos):
                new_name = '_%s' % field.ID
                implicit = True
            else:
                new_name = field.name
                implicit = False
            while new_name in unique:
                new_name = field.name + " %d" % i
                i += 1
            if (field.name is not None and
                not implicit and
                new_name != field.name):
                vo_warn(W33, (field.name, new_name), field._config, field._pos)
            field._unique_name = new_name
            unique[new_name] = field.name

    def _setup(self, config, pos):
        if self.values._ref is not None:
            self.values.ref = self.values._ref
        self.converter = converters.get_converter(self, config, pos)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the FIELD_ element, used for cross-referencing.
        May be ``None`` or a string conforming to XML ID_ syntax.
        """)

    def _set_name(self, name):
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name
    def _del_name(self):
        self._name = None
    name = property(
        attrgetter('_name'), _set_name, _del_name,
        """
        An optional name for the FIELD_.
        """)

    def _set_datatype(self, datatype):
        if datatype is None:
            if self._config.get('version_1_1_or_later'):
                vo_raise(E10, self._element_name, self._config, self._pos)
            else:
                datatype = 'char'
        if datatype not in converters.converter_mapping:
            vo_raise(E06, (datatype, self.ID), self._config, self._pos)
        self._datatype = datatype
    datatype = property(
        attrgetter('_datatype'), _set_datatype, None,
        """
        [*required*] The datatype of the column.  Valid values (as
        defined by the spec) are:

          'boolean', 'bit', 'unsignedByte', 'short', 'int', 'long',
          'char', 'unicodeChar', 'float', 'double', 'floatComplex', or
          'doubleComplex'

        Many VOTABLE files in the wild use 'string' instead of 'char',
        so that is also a valid option, though 'string' will always be
        converted to 'char' when writing the file back out.
        """)

    def _set_precision(self, precision):
        if precision is not None and not re.match("^[FE]?[0-9]+$", precision):
            vo_raise(E11, precision, self._config, self._pos)
        self._precision = precision
    def _del_precision(self):
        self._precision = None
    precision = property(
        attrgetter('_precision'), _set_precision, _del_precision,
        """
        Along with :attr:`width`, defines the `numerical accuracy`_
        associated with the data.  These values are used to limit the
        precision when writing floating point values back to the XML
        file.  Otherwise, it is purely informational -- the Numpy
        recarray containing the data itself does not use this
        information.
        """)

    def _set_width(self, width):
        if width is not None:
            width = int(width)
            if width <= 0:
                vo_raise(E12, width, self._config, self._pos)
        self._width = width
    def _del_width(self):
        self._width = None
    width = property(
        attrgetter('_width'), _set_width, _del_width,
        """
        Along with :attr:`precision`, defines the `numerical
        accuracy`_ associated with the data.  These values are used to
        limit the precision when writing floating point values back to
        the XML file.  Otherwise, it is purely informational -- the
        Numpy recarray containing the data itself does not use this
        information.
        """)

    # ref on FIELD and PARAM behave differently than elsewhere -- here
    # they're just informational, such as to refer to a coordinate
    # system.
    def _set_ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        On FIELD_ elements, ref is used only for informational
        purposes, for example to refer to a COOSYS_ element.
        """
        )

    def _set_ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            check_ucd(ucd, self._config, self._pos)
        self._ucd = ucd
    def _del_ucd(self):
        self._ucd = None
    ucd = property(
        attrgetter('_ucd'), _set_ucd, _del_ucd,
        """
        The `unified content descriptor`_ for the FIELD_.
        """
        )

    def _set_utype(self, utype):
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the FIELD_.
        """
        )

    def _set_unit(self, unit):
        # TODO: Validate unit more accurately
        xmlutil.check_token(unit, 'unit', self._config, self._pos)
        self._unit = unit
    def _del_unit(self):
        self._unit = None
    unit = property(
        attrgetter('_unit'), _set_unit, _del_unit,
        """
        A string specifying the units_ for the FIELD_.
        """)

    def _set_arraysize(self, arraysize):
        if (arraysize is not None and
            not re.match("^([0-9]+x)*[0-9]*[*]?(s\W)?$", arraysize)):
            vo_raise(E13, arraysize, self._config, self._pos)
        self._arraysize = arraysize
    def _del_arraysize(self):
        self._arraysize = None
    arraysize = property(
        attrgetter('_arraysize'), _set_arraysize, _del_arraysize,
        """
        Specifies the size of the multidimensional array if this
        FIELD_ contains more than a single value.

        See `multidimensional arrays`_.
        """)

    def _set_type(self, type):
        self._type = type
    def _del_type(self):
        self._type = None
    type = property(
        attrgetter('_type'), _set_type, _del_type,
        """
        The type attribute on FIELD_ elements is reserved for future
        extensions.
        """)

    def _set_xtype(self, xtype):
        if xtype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('xtype', 'FIELD', '1.2'), config, pos)
        check_string(xtype, 'xtype', self._config, self._pos)
        self._xtype = xtype
    def _del_xtype(self):
        self._xtype = None
    xtype = property(
        attrgetter('_xtype'), _set_xtype, _del_xtype,
        """
        Extended data type information.
        """)

    def _set_values(self, values):
        assert values is None or isinstance(values, Values)
        self._values = values
    def _del_values(self):
        self._values = None
    values = property(
        attrgetter('_values'), _set_values, _del_values,
        """
        A :class:`Values` instance (or ``None``) defining the domain of the
        column.
        """)

    links = property(
        attrgetter('_links'), None, None,
        """
        A list of :class:`Link` instances used to reference more details
        about the meaning of the FIELD_.  This is purely informational
        and is not used by the vo package.
        """)

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start:
                if tag == 'VALUES':
                    self.values.__init__(
                        self._votable, self, config=config, pos=pos, **data)
                    self.values.parse(iterator, config)
                elif tag == 'LINK':
                    link = Link(config=config, pos=pos, **data)
                    self.links.append(link)
                    link.parse(iterator, config)
                elif tag == 'DESCRIPTION':
                    warn_unknown_attrs(
                        'DESCRIPTION', data.iterkeys(), config, pos)
                elif tag != self._element_name:
                    self._add_unknown_tag(iterator, tag, data, config, pos)
            else:
                if tag == 'DESCRIPTION':
                    if self.description is not None:
                        warn_or_raise(
                            W17, W17, self._element_name, config, pos)
                    self.description = data or None
                elif tag == self._element_name:
                    break

        if self.description is not None:
            self.title = " ".join(x.strip() for x in
                                  self.description.split("\n"))
        else:
            self.title = self.name

        self._setup(config, pos)

        return self

    def to_xml(self, w, **kwargs):
        with w.tag(self._element_name,
                   attrib=xmlutil.object_attrs(self, self._attr_list)):
            if self.description is not None:
                w.element("DESCRIPTION", self.description, wrap=True)
            if not self.values.is_defaults():
                self.values.to_xml(w, **kwargs)
            for link in self.links:
                link.to_xml(w, **kwargs)


class Param(Field):
    """
    A class to represent the PARAM_ element, which are constant-valued
    columns in the data.

    :class:`Param` objects are a subclass of :class:`Field`, and have
    all of its methods and members.  Additionally, it defines :attr:`value`.
    """
    _attr_list_11 = Field._attr_list_11 + ['value']
    _attr_list_12 = Field._attr_list_12 + ['value']
    _element_name = 'PARAM'

    def __init__(self, votable, ID=None, name=None, value=None, datatype=None,
                 arraysize=None, ucd=None, unit=None, width=None,
                 precision=None, utype=None, type=None, id=None, config={},
                 pos=None, **extra):
        self._value = value
        Field.__init__(self, votable, ID=ID, name=name, datatype=datatype,
                       arraysize=arraysize, ucd=ucd, unit=unit,
                       precision=precision, utype=utype, type=type,
                       id=id, config=config, pos=pos, **extra)

    def _set_value(self, value):
        if value is None:
            vo_raise(E14, (), self._config, self._pos)
        if isinstance(value, basestring):
            self._value = self.converter.parse(
                value, self._config, self._pos)[0]
        else:
            self._value = value
    value = property(
        attrgetter('_value'), _set_value, None,
        """
        [*required*] The constant value of the parameter.  Its type is
        determined by the :attr:`~Field.datatype` member.
        """)

    def _setup(self, config, pos):
        Field._setup(self, config, pos)
        self.value = self._value

    def to_xml(self, w, **kwargs):
        tmp_value = self._value
        self._value = self.converter.output(tmp_value, False)
        # We must always have a value
        if self._value in (None, ''):
            self._value = " "
        Field.to_xml(self, w, **kwargs)
        self._value = tmp_value


class CooSys(SimpleElement):
    """
    A class representing the COOSYS_ element, which defines a
    coordinate system.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _attr_list = ['ID', 'equinox', 'epoch', 'system']
    _element_name = 'COOSYS'

    def __init__(self, ID=None, equinox=None, epoch=None, system=None, id=None,
                 config={}, pos=None, **extra):
        self._config = config
        self._pos = pos

        if config.get('version_1_2_or_later'):
            warn_or_raise(W27, W27, (), config, pos)

        SimpleElement.__init__(self)

        self.ID      = resolve_id(ID, id, config, pos)
        self.equinox = equinox
        self.epoch   = epoch
        self.system  = system

        warn_unknown_attrs('COOSYS', extra.iterkeys(), config, pos)

    def _set_ID(self, ID):
        if self._config.get('version_1_1_or_later'):
            if ID is None:
                vo_raise(E15, (), self._config, self._pos)
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    ID = property(
        attrgetter('_ID'), _set_ID, None,
        """
        [*required*] The XML ID of the COOSYS_ element, used for
        cross-referencing.  May be ``None`` or a string conforming to
        XML ID_ syntax.
        """)

    def _set_system(self, system):
        if system not in ('eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5',
                          'galactic', 'supergalactic', 'xy', 'barycentric',
                          'geo_app'):
            warn_or_raise(E16, E16, system, self._config, self._pos)
        self._system = system
    def _del_system(self):
        self._system = None
    system = property(
        attrgetter('_system'), _set_system, _del_system,
        """
        Specifies the type of coordinate system.  Valid choices are:

          'eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5', 'galactic',
          'supergalactic', 'xy', 'barycentric', or 'geo_app'
        """)

    def _set_equinox(self, equinox):
        check_astroyear(equinox, 'equinox', self._config, self._pos)
        self._equinox = equinox
    def _del_equinox(self):
        self._equinox = None
    equinox = property(
        attrgetter('_equinox'), _set_equinox, _del_equinox,
        """
        A parameter required to fix the equatorial or ecliptic systems
        (as e.g. "J2000" as the default "eq_FK5" or "B1950" as the
        default "eq_FK4").
        """)

    def _set_epoch(self, epoch):
        check_astroyear(epoch, 'epoch', self._config, self._pos)
        self._epoch = epoch
    def _del_epoch(self):
        self._epoch = None
    epoch = property(
        attrgetter('_epoch'), _set_epoch, _del_epoch,
        """
        Specifies the epoch of the positions.  It must be a string
        specifying an astronomical year.
        """)


class FieldRef(SimpleElement):
    """
    A class representing the FIELDref_ element, which is used inside
    of GROUP_ elements to refer to FIELD_ elements defined elsewhere.
    """
    _attr_list_11 = ['ref']
    _attr_list_12 = _attr_list_11 + ['ucd', 'utype']
    _element_name = "FIELDref"

    def __init__(self, table, ref, ucd=None, utype=None, config={}, pos=None,
                 **extra):
        """
        *table* is the :class:`Table` object that this :class:`FieldRef`
        is a member of.

        *ref* is the ID to reference a :class:`Field` object defined
        elsewhere.
        """
        self._config = config
        self._pos = pos

        SimpleElement.__init__(self)
        self._table = table
        self.ref    = ref
        self.ucd    = ucd
        self.utype  = utype

        if config.get('version_1_2_or_later'):
            self._attr_list = self._attr_list_12
        else:
            self._attr_list = self._attr_list_11
            if ucd is not None:
                warn_unknown_attrs(self._element_name, ['ucd'], config, pos)
            if utype is not None:
                warn_unknown_attrs(self._element_name, ['utype'], config, pos)

    def _set_ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        The ID_ of the FIELD_ that this FIELDref_ references.
        """)

    def _set_ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            if not self._config.get('version_1_2_or_later'):
                warn_or_raise(
                    W28, W28, ('ucd', 'FIELDref', '1.2'), config, pos)
            check_ucd(ucd, self._config, self._pos)
        self._ucd = ucd
    def _del_ucd(self):
        self._ucd = None
    ucd = property(
        attrgetter('_ucd'), _set_ucd, _del_ucd,
        """
        The `unified content descriptor`_ for the FIELDref_.
        """
        )

    def _set_utype(self, utype):
        if utype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('utype', 'FIELDref', '1.2'), config, pos)
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the FIELDref_.
        """
        )

    def get_ref(self):
        """
        Lookup the :class:`Field` instance that this :class:`FieldRef`
        references.
        """
        for field in self._table._votable.iter_fields_and_params():
            if isinstance(field, Field) and field.ID == self.ref:
                return field
        vo_raise(
            "No field named '%s'" % self.ref,
            self._config, self._pos, KeyError)


class ParamRef(SimpleElement):
    """
    A class representing the PARAMref_ element, which is used inside
    of GROUP_ elements to refer to PARAM_ elements defined elsewhere.

    The keyword arguments correspond to setting members of the same
    name, documented below.

    It contains the following publicly-accessible members:

      *ref*: An XML ID refering to a <PARAM> element.
    """
    _attr_list_11 = ['ref']
    _attr_list_12 = _attr_list_11 + ['ucd', 'utype']
    _element_name = "PARAMref"

    def __init__(self, table, ref, ucd=None, utype=None, config={}, pos=None):
        self._config = config
        self._pos = pos

        Element.__init__(self)
        self._table = table
        self.ref    = ref
        self.ucd    = ucd
        self.utype  = utype

        if config.get('version_1_2_or_later'):
            self._attr_list = self._attr_list_12
        else:
            self._attr_list = self._attr_list_11
            if ucd is not None:
                warn_unknown_attrs(self._element_name, ['ucd'], config, pos)
            if utype is not None:
                warn_unknown_attrs(self._element_name, ['utype'], config, pos)

    def _set_ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        The ID_ of the PARAM_ that this PARAMref_ references.
        """
        )

    def _set_ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            if not self._config.get('version_1_2_or_later'):
                warn_or_raise(
                    W28, W28, ('ucd', 'PARAMref', '1.2'), config, pos)
            check_ucd(ucd, self._config, self._pos)
        self._ucd = ucd
    def _del_ucd(self):
        self._ucd = None
    ucd = property(
        attrgetter('_ucd'), _set_ucd, _del_ucd,
        """
        The `unified content descriptor`_ for the PARAMref_.
        """
        )

    def _set_utype(self, utype):
        if utype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('utype', 'PARAMref', '1.2'), config, pos)
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the PARAMref_.
        """
        )

    def get_ref(self):
        """
        Lookup the :class:`Param` instance that this :class:`PARAMref`
        references.
        """
        for param in self._table._votable.iter_fields_and_params():
            if isinstance(param, Param) and param.ID == self.ref:
                return param
        vo_raise(
            "No params named '%s'" % self.ref,
            self._config, self._pos, KeyError)


class Group(Element):
    """
    Stores information about the grouping of FIELD_ and PARAM_
    elements.

    This information is currently ignored by the vo package---that is
    the columns in the recarray are always flat---but the grouping
    information is stored so that it can be written out again to the
    XML file.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """

    def __init__(self, table, ID=None, name=None, ref=None, ucd=None,
                 utype=None, id=None, config={}, pos=None, **extra):
        self._config     = config
        self._pos        = pos

        Element.__init__(self)
        self._table = table

        self.ID          = (resolve_id(ID, id, config, pos)
                            or xmlutil.fix_id(name, config, pos))
        self.name        = name
        self.ref         = ref
        self.ucd         = ucd
        self.utype       = utype
        self.description = None

        self._entries = util.HomogeneousList(
            (FieldRef, ParamRef, Group, Param))

        warn_unknown_attrs('GROUP', extra.iterkeys(), config, pos)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the GROUP_ element.  May be ``None`` or a string
        conforming to XML ID_ syntax.
        """)

    def _set_name(self, name):
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name
    def _del_name(self):
        self._name = None
    name = property(
        attrgetter('_name'), _set_name, _del_name,
        """
        An optional name for the grouping.
        """)

    def _set_ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        Currently ignored, as it's not clear from the spec how this is
        meant to work.
        """)

    def _set_ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            check_ucd(ucd, self._config)
        self._ucd = ucd
    def _del_ucd(self):
        self._ucd = None
    ucd = property(
        attrgetter('_ucd'), _set_ucd, _del_ucd,
        """
        The `unified content descriptor`_ for the GROUP_.
        """)

    def _set_utype(self, utype):
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the GROUP_.
        """)

    def _set_description(self, description):
        self._description = description
    def _del_description(self):
        self._description = None
    description = property(
        attrgetter('_description'), _set_description, _del_description,
        """
        An optional string describing the GROUP_.  Corresponds to the
        DESCRIPTION_ element.
        """)

    entries = property(
        attrgetter('_entries'), None, None,
        """
        A list of members of the GROUP_.  This list may only contain
        objects of type :class:`Param`, :class:`Group`,
        :class:`ParamRef` and :class:`FieldRef`.
        """)

    def _add_fieldref(self, iterator, tag, data, config, pos):
        fieldref = FieldRef(self._table, config=config, pos=pos, **data)
        self.entries.append(fieldref)

    def _add_paramref(self, iterator, tag, data, config, pos):
        paramref = ParamRef(self._table, config=config, pos=pos, **data)
        self.entries.append(paramref)

    def _add_param(self, iterator, tag, data, config, pos):
        if isinstance(self._table, VOTableFile):
            votable = self._table
        else:
            votable = self._table._votable
        param = Param(votable, config=config, pos=pos, **data)
        self.entries.append(param)
        param.parse(iterator, config)

    def _add_group(self, iterator, tag, data, config, pos):
        group = Group(self._table, config=config, pos=pos, **data)
        self.entries.append(group)
        group.parse(iterator, config)

    def parse(self, iterator, config):
        tag_mapping = {
            'FIELDref'    : self._add_fieldref,
            'PARAMref'    : self._add_paramref,
            'PARAM'       : self._add_param,
            'GROUP'       : self._add_group,
            'DESCRIPTION' : self._ignore_add}

        for start, tag, data, pos in iterator:
            if start:
                tag_mapping.get(tag, self._add_unknown_tag)(
                    iterator, tag, data, config, pos)
            else:
                if tag == 'DESCRIPTION':
                    if self.description is not None:
                        warn_or_raise(W17, W17, 'GROUP', config, pos)
                    self.description = data or None
                elif tag == 'GROUP':
                    break
        return self

    def to_xml(self, w, **kwargs):
        with w.tag(
            'GROUP',
            attrib=xmlutil.object_attrs(
                self, ['ID', 'name', 'ref', 'ucd', 'utype'])):
            if self.description is not None:
                w.element("DESCRIPTION", self.description, wrap=True)
            for entry in self.entries:
                entry.to_xml(w, **kwargs)

    def iter_fields_and_params(self):
        """
        Recursively iterate over all :class:`Param` elements in this
        :class:`Group`.
        """
        for entry in self.entries:
            if isinstance(entry, Param):
                yield entry
            elif isinstance(entry, Group):
                for field in entry.iter_fields_and_params():
                    yield field

    def iter_groups(self):
        """
        Recursively iterate over all sub-:class:`Group` instances in
        this :class:`Group`.
        """
        for entry in self.entries:
            if isinstance(entry, Group):
                yield entry
                for group in entry.iter_groups():
                    yield group


class Table(Element):
    """
    A class to store a TABLE_ element, which optionally contains data.

    It contains the following publicly-accessible members, all of
    which are mutable:

        *array*: A Numpy recarray of the data itself, where each row
        is a row of votable data, and columns are named and typed
        based on the <FIELD> elements of the table.

        *mask*: A Numpy recarray of only boolean values, set to *True*
        wherever a value is undefined.

    If the Table contains no data, (for example, its enclosing
    :class:`Resource` has :attr:`~Resource.type` == 'meta') *array*
    and *mask* will be zero-length arrays.

    .. note::
        In a future version of the vo package, the *array* and *mask*
        elements will likely be combined into a single Numpy masked
        record array.  However, there are a number of deficiencies the
        current implementation of Numpy that prevent this.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, votable, ID=None, name=None, ref=None, ucd=None,
                 utype=None, nrows=None, id=None, config={}, pos=None,
                 **extra):
        self._config = config
        self._pos = pos
        self._empty = False

        Element.__init__(self)
        self._votable = votable

        self.ID = (resolve_id(ID, id, config, pos)
                   or xmlutil.fix_id(name, config, pos))
        self.name = name
        xmlutil.check_id(ref, 'ref', config, pos)
        self._ref = ref
        self.ucd = ucd
        self.utype = utype
        if nrows is not None:
            nrows = int(nrows)
            assert nrows >= 0
        self._nrows = nrows
        self.description = None
        self.format = 'tabledata'

        self._fields = util.HomogeneousList(Field)
        self._params = util.HomogeneousList(Param)
        self._groups = util.HomogeneousList(Group)
        self._links  = util.HomogeneousList(Link)
        self._infos  = util.HomogeneousList(Info)

        self.array = np.array([])
        self.mask  = np.array([])

        warn_unknown_attrs('TABLE', extra.iterkeys(), config, pos)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the TABLE_ element, used for cross-referencing.
        May be ``None`` or a string conforming to XML ID_ syntax.
        """)

    def _set_name(self, name):
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name
    def _del_name(self):
        self._name = None
    name = property(
        attrgetter('_name'), _set_name, _del_name,
        """
        An optional name for the table.
        """)

    def _set_ref(self, ref):
        # When the ref changes, we want to verify that it will work
        # by actually going and looking for the referenced table.
        # If found, set a bunch of properties in this table based
        # on the other one.
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        if ref is not None:
            try:
                table = self._votable.get_table_by_id(ref, before=self)
            except KeyError:
                warn_or_raise(
                    W43, W43, ('TABLE', self.ref), self._config, self._pos)
                ref = None
            else:
                self._fields = table.fields
                self._params = table.params
                self._groups = table.groups
                self._links  = table.links
        else:
            del self._fields[:]
            del self._params[:]
            del self._groups[:]
            del self._links[:]
        self._ref = ref
    def _del_ref(self):
        self._ref = None
    ref = property(
        attrgetter('_ref'), _set_ref, _del_ref,
        """
        Refer to another TABLE, previously defined, by the *ref* ID_
        for all metadata (FIELD_, PARAM_ etc.) information.
        """)

    def _set_ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            check_ucd(ucd, self._config)
        self._ucd = ucd
    def _del_ucd(self):
        self._ucd = None
    ucd = property(
        attrgetter('_ucd'), _set_ucd, _del_ucd,
        """
        The `unified content descriptor`_ for the TABLE_.
        """
        )

    def _set_format(self, format):
        format = format.lower()
        if format == 'fits':
            vo_raise("fits format can not be written out, only read.",
                     self._config, self._pos, NotImplementedError)
        if format not in ('tabledata', 'binary'):
            vo_raise("Invalid format '%s'" % format,
                     self._config, self._pos)
        self._format = format
    format = property(
        attrgetter('_format'), _set_format, None,
        """
        [*required*] The serialization format of the table.  Must be
        one of:

          'tabledata' (TABLEDATA_), 'binary' (BINARY_), 'fits' (FITS_).

        Note that the 'fits' format, since it requires an external
        file, can not be written out.  Any file read in with 'fits'
        format will be read out, by default, in 'tabledata' format.
        """)

    nrows = property(
        attrgetter('_nrows'), None, None,
        """
        [*immutable*] The number of rows in the table, as specified in
        the XML file.
        """)

    def _set_description(self, description):
        self._description = description
    def _del_description(self):
        self._description = None
    description = property(
        attrgetter('_description'), _set_description, _del_description,
        """
        An optional string describing the TABLE_.  Corresponds to the
        DESCRIPTION_ element.
        """)

    fields = property(
        attrgetter('_fields'), None, None,
        """
        A list of :class:`Field` objects describing the types of each
        of the data columns.
        """)

    params = property(
        attrgetter('_params'), None, None,
        """
        A list of parameters (constant-valued columns) for the
        table.  Must contain only :class:`Param` objects.
        """)

    groups = property(
        attrgetter('_groups'), None, None,
        """
        A list of :class:`Group` objects describing how the columns
        and parameters are grouped.  Currently this information is
        only kept around for round-tripping and informational
        purposes.
        """)

    links = property(
        attrgetter('_links'), None, None,
        """
        A list of :class:`Link` objects (pointers to other documents
        or servers through a URI) for the table.
        """)

    infos = property(
        attrgetter('_infos'), None, None,
        """
        A list of :class:`Info` objects for the table.  Allows for
        post-operational diagnostics.
        """)

    def is_empty(self):
        """
        Returns True if this table doesn't contain any real data
        because it was skipped over by the parser (through use of the
        `table_number` kwarg).
        """
        return self._empty

    def create_arrays(self, nrows=0, config={}):
        """
        Create new arrays to hold the data based on the current set of
        fields, and store them in the *array* and *mask* member
        variables.  Any data in existing arrays will be lost.

        *nrows*, if provided, is the number of rows to allocate.
        """
        if nrows is None:
            nrows = 0

        fields = self.fields

        if len(fields) == 0:
            array = np.recarray((nrows,), dtype='O')
            mask = np.zeros((nrows,), dtype='b')
        else:
            # for field in fields: field._setup(config)
            Field.uniqify_names(fields)
            ids = [x.ID for x in fields]
            names = [x._unique_name for x in fields]
            formats = [x.converter.format for x in fields]

            descr = np.format_parser(formats, ids, names).dtype
            array = np.recarray((nrows,), dtype=descr)
            descr_mask = []
            for d in descr.descr:
                new_type = (d[1][1] == 'O' and 'O') or 'bool'
                if len(d) == 2:
                    descr_mask.append((d[0], new_type))
                elif len(d) == 3:
                    descr_mask.append((d[0], new_type, d[2]))
            mask = np.zeros((nrows,), dtype=descr_mask)

        self.array = array
        self.mask = mask

    def _resize_strategy(self, size):
        """
        Return a new (larger) size based on size, used for
        reallocating an array when it fills up.  This is in its own
        function so the resizing strategy can be easily replaced.
        """
        # Once we go beyond 0, make a big step -- after that use a
        # factor of 1.5 to help keep memory usage compact
        if size == 0:
            return 512
        return int(ceil(size * RESIZE_AMOUNT))

    def _add_field(self, iterator, tag, data, config, pos):
        field = Field(self._votable, config=config, pos=pos, **data)
        self.fields.append(field)
        field.parse(iterator, config)

    def _add_param(self, iterator, tag, data, config, pos):
        param = Param(self._votable, config=config, pos=pos, **data)
        self.params.append(param)
        param.parse(iterator, config)

    def _add_group(self, iterator, tag, data, config, pos):
        group = Group(self, config=config, pos=pos, **data)
        self.groups.append(group)
        group.parse(iterator, config)

    def _add_link(self, iterator, tag, data, config, pos):
        link = Link(config=config, pos=pos, **data)
        self.links.append(link)
        link.parse(iterator, config)

    def parse(self, iterator, config):
        columns = config.get('columns')

        # If we've requested to read in only a specific table, skip
        # all others
        table_number = config.get('table_number')
        current_table_number = config.get('_current_table_number')
        skip_table = False
        if current_table_number is not None:
            config['_current_table_number'] += 1
            if (table_number is not None and
                table_number != current_table_number):
                skip_table = True
                self._empty = True

        if self.ref is not None:
            # This table doesn't have its own datatype descriptors, it
            # just references those from another table.

            # This is to call the property setter to go and get the
            # referenced information
            self.ref = self.ref

            for start, tag, data, pos in iterator:
                if start:
                    if tag == 'DATA':
                        warn_unknown_attrs(
                            'DATA', data.iterkeys(), config, pos)
                        break
                else:
                    if tag == 'TABLE':
                        return self
                    elif tag == 'DESCRIPTION':
                        if self.description is not None:
                            warn_or_raise(W17, W17, 'RESOURCE', config, pos)
                        self.description = data or None
        else:
            tag_mapping = {
                'FIELD'       : self._add_field,
                'PARAM'       : self._add_param,
                'GROUP'       : self._add_group,
                'LINK'        : self._add_link,
                'DESCRIPTION' : self._ignore_add}

            for start, tag, data, pos in iterator:
                if start:
                    if tag == 'DATA':
                        warn_unknown_attrs(
                            'DATA', data.iterkeys(), config, pos)
                        break

                    tag_mapping.get(tag, self._add_unknown_tag)(
                        iterator, tag, data, config, pos)
                else:
                    if tag == 'DESCRIPTION':
                        if self.description is not None:
                            warn_or_raise(W17, W17, 'RESOURCE', config, pos)
                        self.description = data or None
                    elif tag == 'TABLE':
                        # For error checking purposes
                        Field.uniqify_names(self.fields)
                        return self

        self.create_arrays(nrows=self._nrows, config=config)
        fields = self.fields
        names = [x.ID for x in fields]
        # Deal with a subset of the columns, if requested.
        if not columns:
            colnumbers = range(len(fields))
        else:
            if isinstance(columns, basestring):
                columns = [columns]
            columns = np.asarray(columns)
            if issubclass(columns.dtype.type, np.integer):
                if np.any(columns < 0) or np.any(columns > len(fields)):
                    raise ValueError(
                        "Some specified column numbers out of range")
                colnumbers = columns
            elif issubclass(columns.dtype.type, np.character):
                try:
                    colnumbers = [names.index(x) for x in columns]
                except ValueError:
                    raise ValueError(
                        "Columns '%s' not found in fields list" % columns)
            else:
                raise TypeError("Invalid columns list")

        if not skip_table:
            for start, tag, data, pos in iterator:
                if start:
                    if tag == 'TABLEDATA':
                        warn_unknown_attrs(
                            'TABLEDATA', data.iterkeys(), config, pos)
                        self.array, self.mask = self._parse_tabledata(
                            iterator, colnumbers, fields, config)
                        break
                    elif tag == 'BINARY':
                        warn_unknown_attrs(
                            'BINARY', data.iterkeys(), config, pos)
                        self.array, self.mask = self._parse_binary(
                            iterator, colnumbers, fields, config)
                        break
                    elif tag == 'FITS':
                        warn_unknown_attrs(
                            'FITS', data.iterkeys(), config, pos, ['extnum'])
                        try:
                            extnum = int(data.get('extnum', 0))
                            if extnum < 0:
                                raise ValueError()
                        except ValueError:
                            vo_raise(E17, (), config, pos)
                        self.array, self.mask = self._parse_fits(
                            iterator, extnum, config)
                        break
                    else:
                        warn_or_raise(W37, W37, tag, config, pos)
                        break

        for start, tag, data, pos in iterator:
            if not start and tag == 'DATA':
                break

        for start, tag, data, pos in iterator:
            if tag == 'INFO':
                if start:
                    if not config.get('version_1_2_or_later'):
                        warn_or_raise(
                            W26, W26, ('INFO', 'TABLE', '1.2'), config, pos)
                    info = Info(config=config, pos=pos, **data)
                    self.infos.append(info)
                    info.parse(iterator, config)
                else:
                    info.content = data
            elif not start and tag == 'TABLE':
                break

        return self

    def _parse_tabledata(self, iterator, colnumbers, fields, config):
        # Since we don't know the number of rows up front, we'll
        # reallocate the record array to make room as we go.  This
        # prevents the need to scan through the XML twice.  The
        # allocation is by factors of 1.5.
        invalid = config.get('invalid', 'exception')

        array = self.array
        mask = self.mask
        # Need to have only one reference so that we can resize the
        # array
        del self.array
        del self.mask

        parsers = [field.converter.parse for field in fields]
        binparsers = [field.converter.binparse for field in fields]

        numrows = 0
        alloc_rows = len(array)
        colnumbers_bits = [i in colnumbers for i in range(len(fields))]
        row_default = [x.converter.default for x in fields]
        mask_default = [True] * len(fields)
        array_chunk = []
        mask_chunk = []
        chunk_size = config.get('chunk_size', DEFAULT_CHUNK_SIZE)
        for start, tag, data, pos in iterator:
            if tag == 'TR':
                # Now parse one row
                row = row_default[:]
                row_mask = mask_default[:]
                i = 0
                for start, tag, data, pos in iterator:
                    if start:
                        binary = (data.get('encoding', None) == 'base64')
                        warn_unknown_attrs(
                            tag, data.iterkeys(), config, pos, ['encoding'])
                    else:
                        if tag == 'TD':
                            if i >= len(fields):
                                vo_raise(E20, len(fields), config, pos)

                            if colnumbers_bits[i]:
                                try:
                                    if binary:
                                        if IS_PY3K:
                                            rawdata = base64.b64decode(
                                                data.encode('ascii'))
                                        else:
                                            rawdata = data.decode('base64')
                                        buf = io.BytesIO(rawdata)
                                        buf.seek(0)
                                        try:
                                            value, mask_value = binparsers[i](
                                                buf.read)
                                        except Exception as e:
                                            vo_reraise(e, config, pos,
                                                       "(in row %d, col '%s')" %
                                                       (len(array_chunk),
                                                        fields[i].ID))
                                    else:
                                        try:
                                            value, mask_value = parsers[i](
                                                data, config, pos)
                                        except Exception as e:
                                            vo_reraise(e, config, pos,
                                                       "(in row %d, col '%s')" %
                                                       (len(array_chunk),
                                                        fields[i].ID))
                                except Exception as e:
                                    if invalid == 'exception':
                                        vo_reraise(e, config, pos)
                                else:
                                    row[i] = value
                                    row_mask[i] = mask_value
                        elif tag == 'TR':
                            break
                        else:
                            self._add_unknown_tag(
                                iterator, tag, data, config, pos)
                        i += 1

                if i < len(fields):
                    vo_raise(E21, (i, len(fields)), config, pos)

                array_chunk.append(tuple(row))
                mask_chunk.append(tuple(row_mask))

                if len(array_chunk) == chunk_size:
                    while numrows + chunk_size > alloc_rows:
                        alloc_rows = self._resize_strategy(alloc_rows)
                    if alloc_rows != len(array):
                        array.resize((alloc_rows,))
                        mask.resize((alloc_rows,))
                    array[numrows:numrows + chunk_size] = array_chunk
                    mask[numrows:numrows + chunk_size] = mask_chunk
                    numrows += chunk_size
                    array_chunk = []
                    mask_chunk = []

            elif not start and tag == 'TABLEDATA':
                break

        # Now, resize the array to the exact number of rows we need and
        # put the last chunk values in there.
        if len(array_chunk):
            alloc_rows = numrows + len(array_chunk)
            array.resize((alloc_rows,))
            mask.resize((alloc_rows,))
            array[numrows:] = array_chunk
            mask[numrows:] = mask_chunk
            numrows += len(array_chunk)

        if (self.nrows is not None and
            self.nrows >= 0 and
            self.nrows != numrows):
            warn_or_raise(W18, W18, (self.nrows, numrows), config, pos)
        self._nrows = numrows

        return array, mask

    def _parse_binary(self, iterator, colnumbers, fields, config):
        fields = self.fields

        have_local_stream = False
        for start, tag, data, pos in iterator:
            if tag == 'STREAM':
                if start:
                    warn_unknown_attrs(
                        'STREAM', data.iterkeys(), config, pos,
                        ['type', 'href', 'actuate', 'encoding', 'expires',
                         'rights'])
                    if 'href' not in data:
                        have_local_stream = True
                        if data.get('encoding', None) != 'base64':
                            warn_or_raise(
                                W38, W38, data.get('encoding', None),
                                config, pos)
                    else:
                        href = data['href']
                        xmlutil.check_anyuri(href, config, pos)
                        encoding = data.get('encoding', None)
                else:
                    buffer = data
                    break

        if have_local_stream:
            buffer = base64.b64decode(buffer.encode('ascii'))
            string_io = io.BytesIO(buffer)
            string_io.seek(0)
            read = string_io.read
        else:
            if not (href.startswith('http') or
                    href.startswith('ftp') or
                    href.startswith('file')):
                vo_raise(
                    "The vo package only supports remote data through http, " +
                    "ftp or file",
                    self._config, self._pos, NotImplementedError)
            fd = urllib2.urlopen(href)
            if encoding is not None:
                if encoding == 'gzip':
                    import gzip
                    fd = gzip.GzipFile(href, 'rb', fileobj=fd)
                elif encoding == 'base64':
                    fd = codecs.EncodedFile(fd, 'base64')
                else:
                    vo_raise(
                        "Unknown encoding type '%s'" % encoding,
                        self._config, self._pos, NotImplementedError)
            read = fd.read

        def careful_read(length):
            result = read(length)
            if len(result) != length:
                raise EOFError
            return result

        array = self.array
        mask = self.mask
        # Need to have only one reference so that we can resize the
        # array
        del self.array
        del self.mask

        binparsers = [field.converter.binparse for field in fields]

        numrows = 0
        alloc_rows = len(array)
        while True:
            # Resize result arrays if necessary
            if numrows >= alloc_rows:
                alloc_rows = self._resize_strategy(alloc_rows)
                array.resize((alloc_rows,))
                mask.resize((alloc_rows,))

            row_data = []
            row_mask_data = []
            try:
                for i, binparse in enumerate(binparsers):
                    try:
                        value, value_mask = binparse(careful_read)
                    except EOFError:
                        raise
                    except Exception as e:
                        vo_reraise(e, config, pos,
                                   "(in row %d, col '%s')" %
                                   (numrows, fields[i].ID))
                    row_data.append(value)
                    row_mask_data.append(value_mask)
            except EOFError:
                break

            row = [x.converter.default for x in fields]
            row_mask = [False] * len(fields)
            for i in colnumbers:
                row[i] = row_data[i]
                row_mask[i] = row_mask_data[i]

            array[numrows] = tuple(row)
            mask[numrows] = tuple(row_mask)
            numrows += 1

        array = np.resize(array, (numrows,))
        mask = np.resize(mask, (numrows,))

        return array, mask

    def _parse_fits(self, iterator, extnum, config):
        if not _has_pyfits:
            vo_raise(
                "Input file contains FITS data, but pyfits is not installed.",
                config, None, ImportError)

        for start, tag, data, pos in iterator:
            if tag == 'STREAM':
                if start:
                    warn_unknown_attrs(
                        'STREAM', data.iterkeys(), config, pos,
                        ['type', 'href', 'actuate', 'encoding', 'expires',
                         'rights'])
                    href = data['href']
                    encoding = data.get('encoding', None)
                else:
                    break

        if not (href.startswith('http') or
                href.startswith('ftp') or
                href.startswith('file')):
            vo_raise(
                "The vo package only supports remote data through http, "
                "ftp or file",
                self._config, self._pos, NotImplementedError)

        fd = urllib2.urlopen(href)
        if encoding is not None:
            if encoding == 'gzip':
                import gzip
                fd = gzip.GzipFile(href, 'r', fileobj=fd)
            elif encoding == 'base64':
                fd = codecs.EncodedFile(fd, 'base64')
            else:
                vo_raise(
                    "Unknown encoding type '%s'" % encoding,
                    self._config, self._pos, NotImplementedError)

        fits = pyfits.open(fd)

        array = fits[int(extnum)].data
        if array.dtype != self.array.dtype:
            warn_or_raise(W19, W19, (), self._config, self._pos)

        return array, self.mask

    def to_xml(self, w, **kwargs):
        with w.tag(
            'TABLE',
             attrib=xmlutil.object_attrs(
                       self,
                       ('ID', 'name', 'ref', 'ucd', 'utype', 'nrows'))):

            if self.description is not None:
                w.element("DESCRIPTION", self.description, wrap=True)

            for element_set in (self.fields, self.params):
                for element in element_set:
                    element._setup({}, None)

            if self.ref is None:
                for element_set in (self.fields, self.params, self.groups,
                                    self.links):
                    for element in element_set:
                        element.to_xml(w, **kwargs)

            if len(self.array):
                with w.tag('DATA'):
                    if self.format == 'fits':
                        self.format = 'tabledata'

                    if self.format == 'tabledata':
                        self._write_tabledata(w, **kwargs)
                    elif self.format == 'binary':
                        self._write_binary(w, **kwargs)

            if self.ref is None and kwargs['version_1_2_or_later']:
                for element in self._infos:
                    element.to_xml(w, **kwargs)

    def _write_tabledata(self, w, **kwargs):
        fields = self.fields
        array = self.array
        mask = self.mask

        write_null_values = kwargs.get('write_null_values', False)
        with w.tag('TABLEDATA'):
            w._flush()
            if (_has_c_tabledata_writer and
                not kwargs.get('_debug_python_based_parser')):
                fields = [field.converter.output for field in fields]
                indent = len(w._tags) - 1
                iterparser.write_tabledata(w.write, array, mask, fields,
                                           write_null_values, indent, 1 << 8)
            else:
                write = w.write
                indent_spaces = w.get_indentation_spaces()
                if not IS_PY3K:
                    indent_spaces = indent_spaces.encode('ascii')
                tr_start = indent_spaces + "<TR>\n"
                tr_end = indent_spaces + "</TR>\n"
                td = indent_spaces + " <TD>%s</TD>\n"
                td_empty = indent_spaces + " <TD/>\n"
                fields = [(i, field.converter.output)
                          for i, field in enumerate(fields)]
                for row in xrange(len(array)):
                    write(tr_start)
                    array_row = array[row]
                    mask_row = mask[row]
                    for i, output in fields:
                        masked = mask_row[i]
                        if not np.all(masked) or write_null_values:
                            try:
                                val = output(array_row[i], masked)
                            except Exception as e:
                                vo_reraise(e, config, pos,
                                           "(in row %d, col '%s')" %
                                           (row, fields[i].ID))
                            write(td % val)
                        else:
                            write(td_empty)
                    write(tr_end)

    def _write_binary(self, w, **kwargs):
        fields = self.fields
        array = self.array
        mask = self.mask

        with w.tag('BINARY'):
            with w.tag('STREAM', encoding='base64'):
                fields = [(i, field.converter.binoutput)
                          for (i, field) in enumerate(fields)]

                data = io.BytesIO()
                for row in xrange(len(array)):
                    array_row = array[row]
                    array_mask = mask[row]
                    for i, converter in fields:
                        try:
                            chunk = converter(array_row[i], array_mask[i])
                        except Exception as e:
                            vo_reraise(e, config, pos,
                                       "(in row %d, col '%s')" %
                                       (row, fields[i].ID))
                        data.write(chunk)

                w._flush()
                if IS_PY3K:
                    w.write(base64.b64encode(data.getvalue()).decode('ascii'))
                else:
                    w.write(data.getvalue().encode('base64'))

    def iter_fields_and_params(self):
        """
        Recursively iterate over all FIELD and PARAM elements in the
        TABLE.
        """
        for param in self.params:
            yield param
        for field in self.fields:
            yield field
        for group in self.groups:
            for field in group.iter_fields_and_params():
                yield field

    get_field_by_id = _lookup_by_id_factory(
        'iter_fields_and_params', 'FIELD or PARAM',
        """
        Looks up a FIELD or PARAM element by the given ID.
        """)

    get_field_by_id_or_name = _lookup_by_id_or_name_factory(
        'iter_fields_and_params', 'FIELD or PARAM',
        """
        Looks up a FIELD or PARAM element by the given ID or name.
        """)

    def iter_groups(self):
        """
        Recursively iterate over all GROUP elements in the TABLE.
        """
        for group in self.groups:
            yield group
            for g in group.iter_groups():
                yield g

    get_group_by_id = _lookup_by_id_factory(
        'iter_groups', 'GROUP',
        """
        Looks up a GROUP element by the given ID.  Used by the group's
        "ref" attribute
        """)


class Resource(Element):
    """
    A class to store the information in a RESOURCE_ element.  Each
    resource may contain zero-or-more TABLE_ elements and zero-or-more
    nested RESOURCE_ elements.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, name=None, ID=None, utype=None, type='results',
                 id=None, config={}, pos=None, **kwargs):
        self._config           = config
        self._pos              = pos

        Element.__init__(self)
        self.name              = name
        self.ID                = resolve_id(ID, id, config, pos)
        self.utype             = utype
        self.type              = type
        self._extra_attributes = kwargs
        self.description       = None

        self._coordinate_systems = util.HomogeneousList(CooSys)
        self._params             = util.HomogeneousList(Param)
        self._infos              = util.HomogeneousList(Info)
        self._links              = util.HomogeneousList(Link)
        self._tables             = util.HomogeneousList(Table)
        self._resources          = util.HomogeneousList(Resource)

        warn_unknown_attrs('RESOURCE', kwargs.iterkeys(), config, pos)

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the RESOURCE_ element, used for
        cross-referencing.  May be ``None`` or a string conforming to
        XML ID_ syntax.
        """)

    def _set_name(self, name):
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name
    def _del_name(self):
        self._name = None
    name = property(
        attrgetter('_name'), _set_name, _del_name,
        """
        An optional name for the RESOURCE_.
        """)

    def _set_type(self, type):
        if type not in ('results', 'meta'):
            vo_raise(E18, type, self._config, self._pos)
        self._type = type
    type = property(
        attrgetter('_type'), _set_type, None,
        """
        [*required*] The type of the resource.  Must be either:

          - 'results': This resource contains actual result values
            (default)

          - 'meta': This resource contains only datatype descriptions
            (FIELD_ elements), but no actual data.
        """)

    def _set_utype(self, utype):
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype
    def _del_utype(self):
        self._utype = None
    utype = property(
        attrgetter('_utype'), _set_utype, _del_utype,
        """
        The usage-specific or `unique type`_ of the FIELD_.
        """
        )

    extra_attributes = property(
        attrgetter('_extra_attributes'), None, None,
        """
        A dictionary of string keys to string values containing any
        extra attributes of the RESOURCE_ element that are not defined
        in the specification.  (The specification explicitly allows
        for extra attributes here, but nowhere else.)
        """)

    def _set_description(self, description):
        self._description = description
    def _del_description(self):
        self._description = None
    description = property(
        attrgetter('_description'), _set_description, _del_description,
        """
        An optional string describing the RESOURCE_.  Corresponds to the
        DESCRIPTION_ element.
        """
        )

    coordinate_systems = property(
        attrgetter('_coordinate_systems'), None, None,
        """
        A list of coordinate system definitions (COOSYS_ elements) for
        the RESOURCE_.  Must contain only :class:`CooSys` objects.
        """)

    infos = property(
        attrgetter('_infos'), None, None,
        """
        A list of informational parameters (key-value pairs) for the
        resource.  Must only contain :class:`Info` objects.
        """)

    params = property(
        attrgetter('_params'), None, None,
        """
        A list of parameters (constant-valued columns) for the
        resource.  Must contain only :class:`Param` objects.
        """)

    links = property(
        attrgetter('_links'), None, None,
        """
        A list of links (pointers to other documents or servers
        through a URI) for the resource.  Must contain only
        :class:`Link` objects.
        """)

    tables = property(
        attrgetter('_tables'), None, None,
        """
        A list of tables in the resource.  Must contain only
        :class:`Table` objects.
        """)

    resources = property(
        attrgetter('_resources'), None, None,
        """
        A list of nested resources inside this resource.  Must contain
        only :class:`Resource` objects.
        """)

    def _add_table(self, iterator, tag, data, config, pos):
        table = Table(self._votable, config=config, pos=pos, **data)
        self.tables.append(table)
        table.parse(iterator, config)

    def _add_info(self, iterator, tag, data, config, pos):
        info = Info(config=config, pos=pos, **data)
        self.infos.append(info)
        info.parse(iterator, config)

    def _add_param(self, iterator, tag, data, config, pos):
        param = Param(self._votable, config=config, pos=pos, **data)
        self.params.append(param)
        param.parse(iterator, config)

    def _add_coosys(self, iterator, tag, data, config, pos):
        coosys = CooSys(config=config, pos=pos, **data)
        self.coordinate_systems.append(coosys)
        coosys.parse(iterator, config)

    def _add_resource(self, iterator, tag, data, config, pos):
        resource = Resource(config=config, pos=pos, **data)
        self.resources.append(resource)
        resource.parse(self._votable, iterator, config)

    def _add_link(self, iterator, tag, data, config, pos):
        link = Link(config=config, pos=pos, **data)
        self.links.append(link)
        link.parse(iterator, config)

    def parse(self, votable, iterator, config):
        self._votable = votable

        tag_mapping = {
            'TABLE'       : self._add_table,
            'INFO'        : self._add_info,
            'PARAM'       : self._add_param,
            'COOSYS'      : self._add_coosys,
            'RESOURCE'    : self._add_resource,
            'LINK'        : self._add_link,
            'DESCRIPTION' : self._ignore_add
            }

        for start, tag, data, pos in iterator:
            if start:
                tag_mapping.get(tag, self._add_unknown_tag)(
                    iterator, tag, data, config, pos)
            elif tag == 'DESCRIPTION':
                if self.description is not None:
                    warn_or_raise(W17, W17, 'RESOURCE', config, pos)
                self.description = data or None
            elif tag == 'INFO':
                self.infos[-1].content = data

        del self._votable

        return self

    def to_xml(self, w, **kwargs):
        attrs = xmlutil.object_attrs(self, ('ID', 'type', 'utype'))
        attrs.update(self.extra_attributes)
        with w.tag('RESOURCE', attrib=attrs):
            if self.description is not None:
                w.element("DESCRIPTION", self.description, wrap=True)
            for element_set in (self.coordinate_systems, self.params,
                                self.infos, self.links, self.tables,
                                self.resources):
                for element in element_set:
                    element.to_xml(w, **kwargs)

    def iter_tables(self):
        """
        Recursively iterates over all tables in the resource and
        nested resources.
        """
        for table in self.tables:
            yield table
        for resource in self.resources:
            for table in resource.iter_tables():
                yield table

    def iter_fields_and_params(self):
        """
        Recursively iterates over all FIELD_ and PARAM_ elements in
        the resource, its tables and nested resources.
        """
        for param in self.params:
            yield param
        for table in self.tables:
            for param in table.iter_fields_and_params():
                yield param
        for resource in self.resources:
            for param in resource.iter_fields_and_params():
                yield param

    def iter_coosys(self):
        """
        Recursively iterates over all the COOSYS_ elements in the
        resource and nested resources.
        """
        for coosys in self.coordinate_systems:
            yield coosys
        for resource in self.resources:
            for coosys in resource.iter_coosys():
                yield coosys


class VOTableFile(Element):
    """
    A class to represent the top-level VOTABLE_ element.

    The keyword arguments correspond to setting members of the same
    name, documented below.

    *version* is settable at construction time only, since conformance
    tests for building the rest of the structure depend on it.
    """

    def __init__(self, ID=None, id=None, config={}, pos=None, version="1.2"):
        self._config             = config
        self._pos                = pos

        Element.__init__(self)
        self.ID                  = resolve_id(ID, id, config, pos)
        self.description         = None
        self._coordinate_systems = util.HomogeneousList(CooSys)
        self._params             = util.HomogeneousList(Param)
        self._infos              = util.HomogeneousList(Info)
        self._resources          = util.HomogeneousList(Resource)
        self._groups             = util.HomogeneousList(Group)
        version = str(version)
        assert version in ("1.0", "1.1", "1.2")
        self._version            = version

    def _set_ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID
    def _del_ID(self):
        self._ID = None
    ID = property(
        attrgetter('_ID'), _set_ID, _del_ID,
        """
        The XML ID of the VOTABLE_ element, used for
        cross-referencing.  May be ``None`` or a string conforming to
        XML ID_ syntax.
        """)

    def _get_version(self):
        return self._version
    version = property(
        _get_version, None, None,
        """
        The version of the VOTable specification that the file uses.
        """)

    def _set_description(self, description):
        self._description = description
    def _del_description(self):
        self._description = None
    description = property(
        attrgetter('_description'), _set_description, _del_description,
        """
        An optional string describing the VOTABLE_.  Corresponds to
        the DESCRIPTION_ element.
        """)

    coordinate_systems = property(
        attrgetter('_coordinate_systems'), None, None,
        """
        A list of coordinate system descriptions for the file.  Must
        contain only :class:`CooSys` objects.
        """)

    params = property(
        attrgetter('_params'), None, None,
        """
        A list of parameters (constant-valued columns) that apply to
        the entire file.  Must contain only :class:`Param` objects.
        """)

    infos = property(
        attrgetter('_infos'), None, None,
        """
        A list of informational parameters (key-value pairs) for the
        entire file.  Must only contain :class:`Info` objects.
        """)

    resources = property(
        attrgetter('_resources'), None, None,
        """
        A list of resources, in the order they appear in the file.
        Must only contain :class:`Resource` objects.
        """)

    groups = property(
        attrgetter('_groups'), None, None,
        """
        A list of groups, in the order they appear in the file.  Only
        supported as a child of the VOTABLE element in VOTable 1.2 or
        later.
        """)

    def _add_param(self, iterator, tag, data, config, pos):
        param = Param(self, config=config, pos=pos, **data)
        self.params.append(param)
        param.parse(iterator, config)

    def _add_resource(self, iterator, tag, data, config, pos):
        resource = Resource(config=config, pos=pos, **data)
        self.resources.append(resource)
        resource.parse(self, iterator, config)

    def _add_coosys(self, iterator, tag, data, config, pos):
        coosys = CooSys(config=config, pos=pos, **data)
        self.coordinate_systems.append(coosys)
        coosys.parse(iterator, config)

    def _add_info(self, iterator, tag, data, config, pos):
        info = Info(config=config, pos=pos, **data)
        self.infos.append(info)
        info.parse(iterator, config)

    def _add_group(self, iterator, tag, data, config, pos):
        if not config.get('version_1_2_or_later'):
            warn_or_raise(W26, W26, ('GROUP', 'VOTABLE', '1.2'), config, pos)
        group = Group(self, config=config, pos=pos, **data)
        self.groups.append(group)
        group.parse(iterator, config)

    def parse(self, iterator, config):
        config['_current_table_number'] = 0

        for start, tag, data, pos in iterator:
            if start:
                if tag == 'VOTABLE':
                    if 'version' not in data:
                        warn_or_raise(W20, W20, self.version, config, pos)
                        config['version'] = self.version
                    else:
                        config['version'] = self._version = data['version']
                        if config['version'].lower().startswith('v'):
                            warn_or_raise(
                                W29, W29, config['version'], config, pos)
                            self._version = config['version'] = \
                                            config['version'][1:]
                        if config['version'] not in ('1.1', '1.2'):
                            vo_warn(W21, config['version'], config, pos)

                    if 'xmlns' in data:
                        correct_ns = ('http://www.ivoa.net/xml/VOTable/v%s' %
                                      config['version'])
                        if data['xmlns'] != correct_ns:
                            vo_warn(
                                W41, (correct_ns, data['xmlns']), config, pos)
                    else:
                        vo_warn(W42, (), config, pos)

                    break
                else:
                    vo_raise(E19, (), config, pos)
        config['version_1_1_or_later'] = \
            util.version_compare(config['version'], '1.1') >= 0
        config['version_1_2_or_later'] = \
            util.version_compare(config['version'], '1.2') >= 0

        tag_mapping = {
            'PARAM'       : self._add_param,
            'RESOURCE'    : self._add_resource,
            'COOSYS'      : self._add_coosys,
            'INFO'        : self._add_info,
            'DEFINITIONS' : self._add_definitions,
            'DESCRIPTION' : self._ignore_add,
            'GROUP'       : self._add_group}

        for start, tag, data, pos in iterator:
            if start:
                tag_mapping.get(tag, self._add_unknown_tag)(
                    iterator, tag, data, config, pos)
            elif tag == 'DESCRIPTION':
                if self.description is not None:
                    warn_or_raise(W17, W17, 'VOTABLE', config, pos)
                self.description = data or None
            elif tag == 'INFO':
                self.infos[-1].content = data
        return self

    def to_xml(self, fd, write_null_values=False,
               _debug_python_based_parser=False,
               _astropy_version=None):
        """
        Write to an XML file.

        Parameters
        ----------
        fd : str path or writable file-like object
           Where to write the file.

        write_null_values : bool
           When `True`, write the 'null' value (specified in the null
           attribute of the VALUES element for each FIELD) for empty
           values.  When False (default), simply write no value.
        """
        kwargs = {
            'write_null_values': write_null_values,
            'version': self.version,
            'version_1_1_or_later':
                util.version_compare(self.version, '1.1') >= 0,
            'version_1_2_or_later':
                util.version_compare(self.version, '1.2') >= 0,
            '_debug_python_based_parser': _debug_python_based_parser}

        fd = util.convert_to_writable_filelike(fd)
        w = xmlutil.XMLWriter(fd)
        version = self.version
        if _astropy_version is None:
            lib_version = astropy_version
        else:
            lib_version = _astropy_version

        xml_header = """
<?xml version="1.0" encoding="utf-8"?>
<!-- Produced with astropy.io.vo version %(lib_version)s
     http://www.astropy.org/ -->\n"""
        w.write(xml_header.lstrip() % locals())

        with w.tag('VOTABLE',
                   {'version': version,
                    'xmlns:xsi': "http://www.w3.org/2001/XMLSchema-instance",
                    'xsi:noNamespaceSchemaLocation':
                        "http://www.ivoa.net/xml/VOTable/v%s" % version,
                    'xmlns':
                        "http://www.ivoa.net/xml/VOTable/v%s" % version}):
            if self.description is not None:
                w.element("DESCRIPTION", self.description, wrap=True)
            element_sets = [self.coordinate_systems, self.params,
                            self.infos, self.resources]
            if kwargs['version_1_2_or_later']:
                element_sets[0] = self.groups
            for element_set in element_sets:
                for element in element_set:
                    element.to_xml(w, **kwargs)

    def iter_tables(self):
        """
        Iterates over all tables in the VOTable file in a "flat" way,
        ignoring the nesting of resources etc.
        """
        for resource in self.resources:
            for table in resource.iter_tables():
                yield table

    def get_first_table(self):
        """
        Often, you know there is only one table in the file, and
        that's all you need.  This method returns that first table.
        """
        for table in self.iter_tables():
            if not table.is_empty():
                return table
        raise IndexError("No table found in VOTABLE file.")

    get_table_by_id = _lookup_by_id_factory(
        'iter_tables', 'TABLE',
        """
        Looks up a TABLE element by the given ID.  Used by the table
        "ref" attribute.
        """)

    def get_table_by_index(self, idx):
        """
        Get a table by its ordinal position in the file.
        """
        for i, table in enumerate(self.iter_tables()):
            if i == idx:
                return table
        raise IndexError("No table at index %d found in VOTABLE file." % idx)

    def iter_fields_and_params(self):
        """
        Recursively iterate over all FIELD and PARAM elements in the
        VOTABLE file.
        """
        for resource in self.resources:
            for field in resource.iter_fields_and_params():
                yield field

    get_field_by_id = _lookup_by_id_factory(
        'iter_fields_and_params', 'FIELD',
        """
        Looks up a FIELD element by the given ID.  Used by the field's
        "ref" attribute.
        """)

    get_field_by_id_or_name = _lookup_by_id_or_name_factory(
        'iter_fields_and_params', 'FIELD',
        """
        Looks up a FIELD element by the given ID or name.
        """)

    def iter_values(self):
        """
        Recursively iterate over all VALUES_ elements in the VOTABLE
        file.
        """
        for field in self.iter_fields_and_params():
            yield field.values

    get_values_by_id = _lookup_by_id_factory(
        'iter_values', 'VALUES',
        """
        Looks up a VALUES element by the given ID.  Used by the values
        "ref" attribute.
        """)

    def iter_groups(self):
        """
        Recursively iterate over all GROUP elements in the VOTABLE
        file.
        """
        for table in self.iter_tables():
            for group in table.iter_groups():
                yield group

    get_group_by_id = _lookup_by_id_factory(
        'iter_groups', 'GROUP',
        """
        Looks up a GROUP element by the given ID.  Used by the group's
        "ref" attribute
        """)

    def iter_coosys(self):
        """
        Recursively iterate over all COOSYS elements in the VOTABLE
        file.
        """
        for coosys in self.coordinate_systems:
            yield coosys
        for resource in self.resources:
            for coosys in resource.iter_coosys():
                yield coosys

    get_coosys_by_id = _lookup_by_id_factory(
        'iter_coosys', 'COOSYS',
        """Looks up a COOSYS element by the given ID.""")

    def set_all_tables_format(self, format):
        """
        Set the output storage format of all tables in the file.
        """
        for table in self.iter_tables():
            table.format = format
