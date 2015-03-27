# Licensed under a 3-clause BSD style license - see LICENSE.rst
# TODO: Test FITS parsing

from __future__ import absolute_import, division, print_function, unicode_literals
from ...extern import six
from ...extern.six.moves import xrange, urllib

# STDLIB
import base64
import codecs
import io
import re
import sys
import warnings

# THIRD-PARTY
import numpy as np
from numpy import ma

# LOCAL
from .. import fits
from ... import __version__ as astropy_version
from ...utils.collections import HomogeneousList
from ...utils.xml.writer import XMLWriter
from ...utils.exceptions import AstropyDeprecationWarning
from ...utils.misc import InheritDocstrings

from . import converters
from .exceptions import (warn_or_raise, vo_warn, vo_raise, vo_reraise,
    warn_unknown_attrs,
    W06, W07, W08, W09, W10, W11, W12, W13, W15, W17, W18, W19, W20,
    W21, W22, W26, W27, W28, W29, W32, W33, W35, W36, W37, W38, W40,
    W41, W42, W43, W44, W45, W50, W52, W53, E06, E08, E09, E10, E11,
    E12, E13, E14, E15, E16, E17, E18, E19, E20, E21)
from . import ucd as ucd_mod
from . import util
from . import xmlutil

try:
    from . import tablewriter
    _has_c_tabledata_writer = True
except ImportError:
    _has_c_tabledata_writer = False


__all__ = [
    'Link', 'Info', 'Values', 'Field', 'Param', 'CooSys',
    'FieldRef', 'ParamRef', 'Group', 'Table', 'Resource',
    'VOTableFile'
    ]


# The default number of rows to read in each chunk before converting
# to an array.
DEFAULT_CHUNK_SIZE = 256
RESIZE_AMOUNT = 1.5

######################################################################
# FACTORY FUNCTIONS


def _resize(masked, new_size):
    """
    Masked arrays can not be resized inplace, and `np.resize` and
    `ma.resize` are both incompatible with structured arrays.
    Therefore, we do all this.
    """
    new_array = ma.zeros((new_size,), dtype=masked.dtype)
    length = min(len(masked), new_size)
    try:
        # Pre Numpy 1.10 way
        new_array.data[:length] = masked.data[:length]
    except TypeError:
        # Numpy 1.10 and later
        new_array[:length] = masked[:length]
    else:
        if length != 0:
            new_array.mask[:length] = masked.mask[:length]
    return new_array


def _lookup_by_attr_factory(attr, unique, iterator, element_name, doc):
    """
    Creates a function useful for looking up an element by a given
    attribute.

    Parameters
    ----------
    attr : str
        The attribute name

    unique : bool
        Should be `True` if the attribute is unique and therefore this
        should return only one value.  Otherwise, returns a list of
        values.

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
        A function that looks up an element by the given attribute.
    """
    def lookup_by_attr(self, ref, before=None):
        """
        Given a string *ref*, finds the first element in the iterator
        where the given attribute == *ref*.  If *before* is provided,
        will stop searching at the object *before*.  This is
        important, since "forward references" are not allowed in the
        VOTABLE format.
        """
        for element in getattr(self, iterator)():
            if element is before:
                if getattr(element, attr, None) == ref:
                    vo_raise(
                        "%s references itself" % element_name,
                        element._config, element._pos, KeyError)
                break
            if getattr(element, attr, None) == ref:
                yield element

    def lookup_by_attr_unique(self, ref, before=None):
        for element in lookup_by_attr(self, ref, before=before):
            return element
        raise KeyError(
            "No %s with %s '%s' found before the referencing %s" %
            (element_name, attr, ref, element_name))

    if unique:
        lookup_by_attr_unique.__doc__ = doc
        return lookup_by_attr_unique
    else:
        lookup_by_attr.__doc__ = doc
        return lookup_by_attr


def _lookup_by_id_or_name_factory(iterator, element_name, doc):
    """
    Like `_lookup_by_attr_factory`, but looks in both the "ID" and
    "name" attributes.
    """
    def lookup_by_id_or_name(self, ref, before=None):
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

    lookup_by_id_or_name.__doc__ = doc
    return lookup_by_id_or_name


def _get_default_unit_format(config):
    """
    Get the default unit format as specified in the VOTable spec.
    """
    # In the future, this should take into account the VOTable
    # version.
    return 'cds'


def _get_unit_format(config):
    """
    Get the unit format based on the configuration.
    """
    if config.get('unit_format') is None:
        format = _get_default_unit_format(config)
    else:
        format = config['unit_format']
    return format



######################################################################
# ATTRIBUTE CHECKERS
def check_astroyear(year, field, config=None, pos=None):
    """
    Raises a `~astropy.io.votable.exceptions.VOTableSpecError` if
    *year* is not a valid astronomical year as defined by the VOTABLE
    standard.

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


def check_string(string, attr_name, config=None, pos=None):
    """
    Raises a `~astropy.io.votable.exceptions.VOTableSpecError` if
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
    if string is not None and not isinstance(string, six.string_types):
        warn_or_raise(W08, W08, attr_name, config, pos)
        return False
    return True


def resolve_id(ID, id, config=None, pos=None):
    if ID is None and id is not None:
        warn_or_raise(W09, W09, (), config, pos)
        return id
    return ID


def check_ucd(ucd, config=None, pos=None):
    """
    Warns or raises a
    `~astropy.io.votable.exceptions.VOTableSpecError` if *ucd* is not
    a valid `unified content descriptor`_ string as defined by the
    VOTABLE standard.

    Parameters
    ----------
    ucd : str
        A UCD string.

    config, pos : optional
        Information about the source of the value
    """
    if config is None:
        config = {}
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
                vo_raise(W06, (ucd, six.text_type(e)), config, pos)
            else:
                vo_warn(W06, (ucd, six.text_type(e)), config, pos)
                return False
    return True


######################################################################
# PROPERTY MIXINS
class _IDProperty(object):
    @property
    def ID(self):
        """
        The XML ID_ of the element.  May be `None` or a string
        conforming to XML ID_ syntax.
        """
        return self._ID

    @ID.setter
    def ID(self, ID):
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID

    @ID.deleter
    def ID(self):
        self._ID = None


class _NameProperty(object):
    @property
    def name(self):
        """An optional name for the element."""
        return self._name

    @name.setter
    def name(self, name):
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name

    @name.deleter
    def name(self):
        self._name = None


class _XtypeProperty(object):
    @property
    def xtype(self):
        """Extended data type information."""
        return self._xtype

    @xtype.setter
    def xtype(self, xtype):
        if xtype is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(
                W28, W28, ('xtype', self._element_name, '1.2'),
                self._config, self._pos)
        check_string(xtype, 'xtype', self._config, self._pos)
        self._xtype = xtype

    @xtype.deleter
    def xtype(self):
        self._xtype = None


class _UtypeProperty(object):
    _utype_in_v1_2 = False

    @property
    def utype(self):
        """The usage-specific or `unique type`_ of the element."""
        return self._utype

    @utype.setter
    def utype(self, utype):
        if (self._utype_in_v1_2 and
            utype is not None and
            not self._config.get('version_1_2_or_later')):
            warn_or_raise(
                W28, W28, ('utype', self._element_name, '1.2'),
                self._config, self._pos)
        check_string(utype, 'utype', self._config, self._pos)
        self._utype = utype

    @utype.deleter
    def utype(self):
        self._utype = None


class _UcdProperty(object):
    _ucd_in_v1_2 = False

    @property
    def ucd(self):
        """The `unified content descriptor`_ for the element."""
        return self._ucd

    @ucd.setter
    def ucd(self, ucd):
        if ucd is not None and ucd.strip() == '':
            ucd = None
        if ucd is not None:
            if (self._ucd_in_v1_2 and
                not self._config.get('version_1_2_or_later')):
                warn_or_raise(
                    W28, W28, ('ucd', self._element_name, '1.2'),
                    self._config, self._pos)
            check_ucd(ucd, self._config, self._pos)
        self._ucd = ucd

    @ucd.deleter
    def ucd(self):
        self._ucd = None


class _DescriptionProperty(object):
    @property
    def description(self):
        """
        An optional string describing the element.  Corresponds to the
        DESCRIPTION_ element.
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @description.deleter
    def description(self):
        self._description = None


######################################################################
# ELEMENT CLASSES
@six.add_metaclass(InheritDocstrings)
class Element(object):
    """
    A base class for all classes that represent XML elements in the
    VOTABLE file.
    """
    def _add_unknown_tag(self, iterator, tag, data, config, pos):
        warn_or_raise(W10, W10, tag, config, pos)

    def _ignore_add(self, iterator, tag, data, config, pos):
        warn_unknown_attrs(tag, six.iterkeys(data), config, pos)

    def _add_definitions(self, iterator, tag, data, config, pos):
        if config.get('version_1_1_or_later'):
            warn_or_raise(W22, W22, (), config, pos)
        warn_unknown_attrs(tag, six.iterkeys(data), config, pos)

    def parse(self, iterator, config):
        """
        For internal use. Parse the XML content of the children of the
        element.

        Parameters
        ----------
        iterator : xml iterator
            An iterator over XML elements as returned by
            `~astropy.utils.xml.iterparser.get_xml_iterator`.

        config : dict
            The configuration dictionary that affects how certain
            elements are read.

        Returns
        -------
        self : Element
            Returns self as a convenience.
        """
        raise NotImplementedError()

    def to_xml(self, w, **kwargs):
        """
        For internal use. Output the element to XML.

        Parameters
        ----------
        w : astropy.utils.xml.writer.XMLWriter object
            An XML writer to write to.

        kwargs : dict
            Any configuration parameters to control the output.
        """
        raise NotImplementedError()


class SimpleElement(Element):
    """
    A base class for simple elements, such as FIELD, PARAM and INFO
    that don't require any special parsing or outputting machinery.
    """
    def __init__(self):
        Element.__init__(self)

    def __repr__(self):
        buff = io.StringIO()
        SimpleElement.to_xml(self, XMLWriter(buff))
        return buff.getvalue()

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start and tag != self._element_name:
                self._add_unknown_tag(iterator, tag, data, config, pos)
            elif tag == self._element_name:
                break

        return self

    def to_xml(self, w, **kwargs):
        w.element(self._element_name,
                  attrib=w.object_attrs(self, self._attr_list))


class SimpleElementWithContent(SimpleElement):
    """
    A base class for simple elements, such as FIELD, PARAM and INFO
    that don't require any special parsing or outputting machinery.
    """
    def __init__(self):
        SimpleElement.__init__(self)

        self._content = None

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start and tag != self._element_name:
                self._add_unknown_tag(iterator, tag, data, config, pos)
            elif tag == self._element_name:
                if data:
                    self.content = data
                break

        return self

    def to_xml(self, w, **kwargs):
        w.element(self._element_name, self._content,
                  attrib=w.object_attrs(self, self._attr_list))

    @property
    def content(self):
        """The content of the element."""
        return self._content

    @content.setter
    def content(self, content):
        check_string(content, 'content', self._config, self._pos)
        self._content = content

    @content.deleter
    def content(self):
        self._content = None


class Link(SimpleElement, _IDProperty):
    """
    LINK_ elements: used to reference external documents and servers through a URI.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _attr_list = ['ID', 'content_role', 'content_type', 'title', 'value',
                  'href', 'action']
    _element_name = 'LINK'

    def __init__(self, ID=None, title=None, value=None, href=None, action=None,
                 id=None, config=None, pos=None, **kwargs):
        if config is None:
            config = {}
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
            'LINK', six.iterkeys(kwargs), config, pos,
            ['content-role', 'content_role', 'content-type', 'content_type',
             'gref'])

    @property
    def content_role(self):
        """
        Defines the MIME role of the referenced object.  Must be one of:

          None, 'query', 'hints', 'doc', 'location' or 'type'
        """
        return self._content_role

    @content_role.setter
    def content_role(self, content_role):
        if ((content_role == 'type' and
             not self._config['version_1_3_or_later']) or
             content_role not in
             (None, 'query', 'hints', 'doc', 'location')):
            vo_warn(W45, (content_role,), self._config, self._pos)
        self._content_role = content_role

    @content_role.deleter
    def content_role(self):
        self._content_role = None

    @property
    def content_type(self):
        """Defines the MIME content type of the referenced object."""
        return self._content_type

    @content_type.setter
    def content_type(self, content_type):
        xmlutil.check_mime_content_type(content_type, self._config, self._pos)
        self._content_type = content_type

    @content_type.deleter
    def content_type(self):
        self._content_type = None

    @property
    def href(self):
        """
        A URI to an arbitrary protocol.  The vo package only supports
        http and anonymous ftp.
        """
        return self._href

    @href.setter
    def href(self, href):
        xmlutil.check_anyuri(href, self._config, self._pos)
        self._href = href

    @href.deleter
    def href(self):
        self._href = None

    def to_table_column(self, column):
        meta = {}
        for key in self._attr_list:
            val = getattr(self, key, None)
            if val is not None:
                meta[key] = val

        column.meta.setdefault('links', [])
        column.meta['links'].append(meta)

    @classmethod
    def from_table_column(cls, d):
        return cls(**d)


class Info(SimpleElementWithContent, _IDProperty, _XtypeProperty,
           _UtypeProperty):
    """
    INFO_ elements: arbitrary key-value pairs for extensions to the standard.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _element_name = 'INFO'
    _attr_list_11 = ['ID', 'name', 'value']
    _attr_list_12 = _attr_list_11 + ['xtype', 'ref', 'unit', 'ucd', 'utype']
    _utype_in_v1_2 = True

    def __init__(self, ID=None, name=None, value=None, id=None, xtype=None,
                 ref=None, unit=None, ucd=None, utype=None,
                 config=None, pos=None, **extra):
        if config is None:
            config = {}
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

        warn_unknown_attrs('INFO', six.iterkeys(extra), config, pos)

    @property
    def name(self):
        """[*required*] The key of the key-value pair."""
        return self._name

    @name.setter
    def name(self, name):
        if name is None:
            warn_or_raise(W35, W35, ('name'), self._config, self._pos)
        xmlutil.check_token(name, 'name', self._config, self._pos)
        self._name = name

    @property
    def value(self):
        """
        [*required*] The value of the key-value pair.  (Always stored
        as a string or unicode string).
        """
        return self._value

    @value.setter
    def value(self, value):
        if value is None:
            warn_or_raise(W35, W35, ('value'), self._config, self._pos)
        check_string(value, 'value', self._config, self._pos)
        self._value = value

    @property
    def content(self):
        """The content inside the INFO element."""
        return self._content

    @content.setter
    def content(self, content):
        check_string(content, 'content', self._config, self._pos)
        self._content = content

    @content.deleter
    def content(self):
        self._content = None

    @property
    def ref(self):
        """
        Refer to another INFO_ element by ID_, defined previously in
        the document.
        """
        return self._ref

    @ref.setter
    def ref(self, ref):
        if ref is not None and not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('ref', 'INFO', '1.2'),
                          self._config, self._pos)
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

    @ref.deleter
    def ref(self):
        self._ref = None

    @property
    def unit(self):
        """A string specifying the units_ for the INFO_."""
        return self._unit

    @unit.setter
    def unit(self, unit):
        if unit is None:
            self._unit = None
            return

        from ... import units as u

        if not self._config.get('version_1_2_or_later'):
            warn_or_raise(W28, W28, ('unit', 'INFO', '1.2'),
                          self._config, self._pos)

        # First, parse the unit in the default way, so that we can
        # still emit a warning if the unit is not to spec.
        default_format = _get_default_unit_format(self._config)
        unit_obj = u.Unit(
            unit, format=default_format, parse_strict='silent')
        if isinstance(unit_obj, u.UnrecognizedUnit):
            warn_or_raise(W50, W50, (unit,),
                          self._config, self._pos)

        format = _get_unit_format(self._config)
        if format != default_format:
            unit_obj = u.Unit(
                unit, format=format, parse_strict='silent')

        self._unit = unit_obj

    @unit.deleter
    def unit(self):
        self._unit = None

    def to_xml(self, w, **kwargs):
        attrib = w.object_attrs(self, self._attr_list)
        if 'unit' in attrib:
            attrib['unit'] = self.unit.to_string('cds')
        w.element(self._element_name, self._content,
                  attrib=attrib)


class Values(Element, _IDProperty):
    """
    VALUES_ element: used within FIELD_ and PARAM_ elements to define the domain of values.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, votable, field, ID=None, null=None, ref=None,
                 type="legal", id=None, config=None, pos=None, **extras):
        if config is None:
            config = {}
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

        warn_unknown_attrs('VALUES', six.iterkeys(extras), config, pos)

    def __repr__(self):
        buff = io.StringIO()
        self.to_xml(XMLWriter(buff))
        return buff.getvalue()

    @property
    def null(self):
        """
        For integral datatypes, *null* is used to define the value
        used for missing values.
        """
        return self._null

    @null.setter
    def null(self, null):
        if null is not None and isinstance(null, six.string_types):
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

    @null.deleter
    def null(self):
        self._null = None

    @property
    def type(self):
        """
        [*required*] Defines the applicability of the domain defined
        by this VALUES_ element.  Must be one of the following
        strings:

          - 'legal': The domain of this column applies in general to
            this datatype. (default)

          - 'actual': The domain of this column applies only to the
            data enclosed in the parent table.
        """
        return self._type

    @type.setter
    def type(self, type):
        if type not in ('legal', 'actual'):
            vo_raise(E08, type, self._config, self._pos)
        self._type = type

    @property
    def ref(self):
        """
        Refer to another VALUES_ element by ID_, defined previously in
        the document, for MIN/MAX/OPTION information.
        """
        return self._ref

    @ref.setter
    def ref(self, ref):
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

    @ref.deleter
    def ref(self):
        self._ref = None

    @property
    def min(self):
        """
        The minimum value of the domain.  See :attr:`min_inclusive`.
        """
        return self._min

    @min.setter
    def min(self, min):
        if hasattr(self._field, 'converter') and min is not None:
            self._min = self._field.converter.parse(min)[0]
        else:
            self._min = min

    @min.deleter
    def min(self):
        self._min = None

    @property
    def min_inclusive(self):
        """When `True`, the domain includes the minimum value."""
        return self._min_inclusive

    @min_inclusive.setter
    def min_inclusive(self, inclusive):
        if inclusive == 'yes':
            self._min_inclusive = True
        elif inclusive == 'no':
            self._min_inclusive = False
        else:
            self._min_inclusive = bool(inclusive)

    @min_inclusive.deleter
    def min_inclusive(self):
        self._min_inclusive = True

    @property
    def max(self):
        """
        The maximum value of the domain.  See :attr:`max_inclusive`.
        """
        return self._max

    @max.setter
    def max(self, max):
        if hasattr(self._field, 'converter') and max is not None:
            self._max = self._field.converter.parse(max)[0]
        else:
            self._max = max

    @max.deleter
    def max(self):
        self._max = None

    @property
    def max_inclusive(self):
        """When `True`, the domain includes the maximum value."""
        return self._max_inclusive

    @max_inclusive.setter
    def max_inclusive(self, inclusive):
        if inclusive == 'yes':
            self._max_inclusive = True
        elif inclusive == 'no':
            self._max_inclusive = False
        else:
            self._max_inclusive = bool(inclusive)

    @max_inclusive.deleter
    def max_inclusive(self):
        self._max_inclusive = True

    @property
    def options(self):
        """
        A list of string key-value tuples defining other OPTION
        elements for the domain.  All options are ignored -- they are
        stored for round-tripping purposes only.
        """
        return self._options

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
                            'MIN', six.iterkeys(data), config, pos,
                            ['value', 'inclusive'])
                    elif tag == 'MAX':
                        if 'value' not in data:
                            vo_raise(E09, 'MAX', config, pos)
                        self.max = data['value']
                        self.max_inclusive = data.get('inclusive', 'yes')
                        warn_unknown_attrs(
                            'MAX', six.iterkeys(data), config, pos,
                            ['value', 'inclusive'])
                    elif tag == 'OPTION':
                        if 'value' not in data:
                            vo_raise(E09, 'OPTION', config, pos)
                        xmlutil.check_token(
                            data.get('name'), 'name', config, pos)
                        self.options.append(
                            (data.get('name'), data.get('value')))
                        warn_unknown_attrs(
                            'OPTION', six.iterkeys(data), config, pos,
                            ['data', 'name'])
                elif tag == 'VALUES':
                    break

        return self

    def is_defaults(self):
        """
        Are the settings on this ``VALUE`` element all the same as the
        XML defaults?
        """
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
            w.element('VALUES', attrib=w.object_attrs(self, ['ref']))
        else:
            with w.tag('VALUES',
                       attrib=w.object_attrs(
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

    def to_table_column(self, column):
        # Have the ref filled in here
        meta = {}
        for key in ['ID', 'null']:
            val = getattr(self, key, None)
            if val is not None:
                meta[key] = val
        if self.min is not None:
            meta['min'] = {
                'value': self.min,
                'inclusive': self.min_inclusive}
        if self.max is not None:
            meta['max'] = {
                'value': self.max,
                'inclusive': self.max_inclusive}
        if len(self.options):
            meta['options'] = dict(self.options)

        column.meta['values'] = meta

    def from_table_column(self, column):
        if not 'values' in column.meta:
            return

        meta = column.meta['values']
        for key in ['ID', 'null']:
            val = meta.get(key, None)
            if val is not None:
                setattr(self, key, val)
        if 'min' in meta:
            self.min = meta['min']['value']
            self.min_inclusive = meta['min']['inclusive']
        if 'max' in meta:
            self.max = meta['max']['value']
            self.max_inclusive = meta['max']['inclusive']
        if 'options' in meta:
            self._options = list(meta['options'].items())


class Field(SimpleElement, _IDProperty, _NameProperty, _XtypeProperty,
            _UtypeProperty, _UcdProperty):
    """
    FIELD_ element: describes the datatype of a particular column of data.

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
                 config=None, pos=None, **extra):
        if config is None:
            config = {}
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
            'float64'       : 'double',
            # The following appear in some Vizier tables
            'unsignedInt'   : 'long',
            'unsignedShort' : 'int'
        }

        datatype_mapping.update(config.get('datatype_mapping', {}))

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
        self._links     = HomogeneousList(Link)
        self.title      = self.name
        self.values     = Values(self._votable, self)
        self.xtype      = xtype

        self._setup(config, pos)

        warn_unknown_attrs(self._element_name, six.iterkeys(extra), config, pos)

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
            if field.name is None:
                new_name = field.ID
                implicit = True
            else:
                new_name = field.name
                implicit = False
            if new_name != field.ID:
                while new_name in unique:
                    new_name = field.name + " %d" % i
                    i += 1

            if (not implicit and
                new_name != field.name):
                vo_warn(W33, (field.name, new_name), field._config, field._pos)
            field._unique_name = new_name
            unique[new_name] = field.name

    def _setup(self, config, pos):
        if self.values._ref is not None:
            self.values.ref = self.values._ref
        self.converter = converters.get_converter(self, config, pos)

    @property
    def datatype(self):
        """
        [*required*] The datatype of the column.  Valid values (as
        defined by the spec) are:

          'boolean', 'bit', 'unsignedByte', 'short', 'int', 'long',
          'char', 'unicodeChar', 'float', 'double', 'floatComplex', or
          'doubleComplex'

        Many VOTABLE files in the wild use 'string' instead of 'char',
        so that is also a valid option, though 'string' will always be
        converted to 'char' when writing the file back out.
        """
        return self._datatype

    @datatype.setter
    def datatype(self, datatype):
        if datatype is None:
            if self._config.get('version_1_1_or_later'):
                warn_or_raise(E10, E10, self._element_name, self._config,
                              self._pos)
            datatype = 'char'
        if datatype not in converters.converter_mapping:
            vo_raise(E06, (datatype, self.ID), self._config, self._pos)
        self._datatype = datatype

    @property
    def precision(self):
        """
        Along with :attr:`width`, defines the `numerical accuracy`_
        associated with the data.  These values are used to limit the
        precision when writing floating point values back to the XML
        file.  Otherwise, it is purely informational -- the Numpy
        recarray containing the data itself does not use this
        information.
        """
        return self._precision

    @precision.setter
    def precision(self, precision):
        if precision is not None and not re.match(r"^[FE]?[0-9]+$", precision):
            vo_raise(E11, precision, self._config, self._pos)
        self._precision = precision

    @precision.deleter
    def precision(self):
        self._precision = None

    @property
    def width(self):
        """
        Along with :attr:`precision`, defines the `numerical
        accuracy`_ associated with the data.  These values are used to
        limit the precision when writing floating point values back to
        the XML file.  Otherwise, it is purely informational -- the
        Numpy recarray containing the data itself does not use this
        information.
        """
        return self._width

    @width.setter
    def width(self, width):
        if width is not None:
            width = int(width)
            if width <= 0:
                vo_raise(E12, width, self._config, self._pos)
        self._width = width

    @width.deleter
    def width(self):
        self._width = None

    # ref on FIELD and PARAM behave differently than elsewhere -- here
    # they're just informational, such as to refer to a coordinate
    # system.
    @property
    def ref(self):
        """
        On FIELD_ elements, ref is used only for informational
        purposes, for example to refer to a COOSYS_ element.
        """
        return self._ref

    @ref.setter
    def ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref

    @ref.deleter
    def ref(self):
        self._ref = None

    @property
    def unit(self):
        """A string specifying the units_ for the FIELD_."""
        return self._unit

    @unit.setter
    def unit(self, unit):
        if unit is None:
            self._unit = None
            return

        from ... import units as u

        # First, parse the unit in the default way, so that we can
        # still emit a warning if the unit is not to spec.
        default_format = _get_default_unit_format(self._config)
        unit_obj = u.Unit(
            unit, format=default_format, parse_strict='silent')
        if isinstance(unit_obj, u.UnrecognizedUnit):
            warn_or_raise(W50, W50, (unit,),
                          self._config, self._pos)

        format = _get_unit_format(self._config)
        if format != default_format:
            unit_obj = u.Unit(
                unit, format=format, parse_strict='silent')

        self._unit = unit_obj

    @unit.deleter
    def unit(self):
        self._unit = None

    @property
    def arraysize(self):
        """
        Specifies the size of the multidimensional array if this
        FIELD_ contains more than a single value.

        See `multidimensional arrays`_.
        """
        return self._arraysize

    @arraysize.setter
    def arraysize(self, arraysize):
        if (arraysize is not None and
            not re.match(r"^([0-9]+x)*[0-9]*[*]?(s\W)?$", arraysize)):
            vo_raise(E13, arraysize, self._config, self._pos)
        self._arraysize = arraysize

    @arraysize.deleter
    def arraysize(self):
        self._arraysize = None

    @property
    def type(self):
        """
        The type attribute on FIELD_ elements is reserved for future
        extensions.
        """
        return self._type

    @type.setter
    def type(self, type):
        self._type = type

    @type.deleter
    def type(self):
        self._type = None

    @property
    def values(self):
        """
        A :class:`Values` instance (or `None`) defining the domain
        of the column.
        """
        return self._values

    @values.setter
    def values(self, values):
        assert values is None or isinstance(values, Values)
        self._values = values

    @values.deleter
    def values(self):
        self._values = None

    @property
    def links(self):
        """
        A list of :class:`Link` instances used to reference more
        details about the meaning of the FIELD_.  This is purely
        informational and is not used by the `astropy.io.votable`
        package.
        """
        return self._links

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
                        'DESCRIPTION', six.iterkeys(data), config, pos)
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
                                  self.description.splitlines())
        else:
            self.title = self.name

        self._setup(config, pos)

        return self

    def to_xml(self, w, **kwargs):
        attrib = w.object_attrs(self, self._attr_list)
        if 'unit' in attrib:
            attrib['unit'] = self.unit.to_string('cds')
        with w.tag(self._element_name, attrib=attrib):
            if self.description is not None:
                w.element('DESCRIPTION', self.description, wrap=True)
            if not self.values.is_defaults():
                self.values.to_xml(w, **kwargs)
            for link in self.links:
                link.to_xml(w, **kwargs)

    def to_table_column(self, column):
        """
        Sets the attributes of a given `astropy.table.Column` instance
        to match the information in this `Field`.
        """
        for key in ['ucd', 'width', 'precision', 'utype', 'xtype']:
            val = getattr(self, key, None)
            if val is not None:
                column.meta[key] = val
        if not self.values.is_defaults():
            self.values.to_table_column(column)
        for link in self.links:
            link.to_table_column(column)
        if self.description is not None:
            column.description = self.description
        if self.unit is not None:
            # TODO: Use units framework when it's available
            column.unit = self.unit
        if isinstance(self.converter, converters.FloatingPoint):
            column.format = self.converter.output_format

    @classmethod
    def from_table_column(cls, votable, column):
        """
        Restores a `Field` instance from a given
        `astropy.table.Column` instance.
        """
        kwargs = {}
        for key in ['ucd', 'width', 'precision', 'utype', 'xtype']:
            val = column.meta.get(key, None)
            if val is not None:
                kwargs[key] = val
        # TODO: Use the unit framework when available
        if column.unit is not None:
            kwargs['unit'] = column.unit
        kwargs['name'] = column.name
        result = converters.table_column_to_votable_datatype(column)
        kwargs.update(result)

        field = cls(votable, **kwargs)

        if column.description is not None:
            field.description = column.description
        field.values.from_table_column(column)
        if 'links' in column.meta:
            for link in column.meta['links']:
                field.links.append(Link.from_table_column(link))

        # TODO: Parse format into precision and width
        return field


class Param(Field):
    """
    PARAM_ element: constant-valued columns in the data.

    :class:`Param` objects are a subclass of :class:`Field`, and have
    all of its methods and members.  Additionally, it defines :attr:`value`.
    """
    _attr_list_11 = Field._attr_list_11 + ['value']
    _attr_list_12 = Field._attr_list_12 + ['value']
    _element_name = 'PARAM'

    def __init__(self, votable, ID=None, name=None, value=None, datatype=None,
                 arraysize=None, ucd=None, unit=None, width=None,
                 precision=None, utype=None, type=None, id=None, config=None,
                 pos=None, **extra):
        self._value = value
        Field.__init__(self, votable, ID=ID, name=name, datatype=datatype,
                       arraysize=arraysize, ucd=ucd, unit=unit,
                       precision=precision, utype=utype, type=type,
                       id=id, config=config, pos=pos, **extra)

    @property
    def value(self):
        """
        [*required*] The constant value of the parameter.  Its type is
        determined by the :attr:`~Field.datatype` member.
        """
        return self._value

    @value.setter
    def value(self, value):
        if value is None:
            value = ""
        if ((six.PY3 and isinstance(value, six.text_type)) or
            (not six.PY3 and isinstance(value, six.string_types))):
            self._value = self.converter.parse(
                value, self._config, self._pos)[0]
        else:
            self._value = value

    def _setup(self, config, pos):
        Field._setup(self, config, pos)
        self.value = self._value

    def to_xml(self, w, **kwargs):
        tmp_value = self._value
        self._value = self.converter.output(tmp_value, False)
        # We must always have a value
        if self._value is None:
            self._value = ""
        Field.to_xml(self, w, **kwargs)
        self._value = tmp_value


class CooSys(SimpleElement):
    """
    COOSYS_ element: defines a coordinate system.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    _attr_list = ['ID', 'equinox', 'epoch', 'system']
    _element_name = 'COOSYS'

    def __init__(self, ID=None, equinox=None, epoch=None, system=None, id=None,
                 config=None, pos=None, **extra):
        if config is None:
            config = {}
        self._config = config
        self._pos = pos

        if config.get('version_1_2_or_later'):
            warn_or_raise(W27, W27, (), config, pos)

        SimpleElement.__init__(self)

        self.ID      = resolve_id(ID, id, config, pos)
        self.equinox = equinox
        self.epoch   = epoch
        self.system  = system

        warn_unknown_attrs('COOSYS', six.iterkeys(extra), config, pos)

    @property
    def ID(self):
        """
        [*required*] The XML ID of the COOSYS_ element, used for
        cross-referencing.  May be `None` or a string conforming to
        XML ID_ syntax.
        """
        return self._ID

    @ID.setter
    def ID(self, ID):
        if self._config.get('version_1_1_or_later'):
            if ID is None:
                vo_raise(E15, (), self._config, self._pos)
        xmlutil.check_id(ID, 'ID', self._config, self._pos)
        self._ID = ID

    @property
    def system(self):
        """
        Specifies the type of coordinate system.  Valid choices are:

          'eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5', 'galactic',
          'supergalactic', 'xy', 'barycentric', or 'geo_app'
        """
        return self._system

    @system.setter
    def system(self, system):
        if system not in ('eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5',
                          'galactic', 'supergalactic', 'xy', 'barycentric',
                          'geo_app'):
            warn_or_raise(E16, E16, system, self._config, self._pos)
        self._system = system

    @system.deleter
    def system(self):
        self._system = None

    @property
    def equinox(self):
        """
        A parameter required to fix the equatorial or ecliptic systems
        (as e.g. "J2000" as the default "eq_FK5" or "B1950" as the
        default "eq_FK4").
        """
        return self._equinox

    @equinox.setter
    def equinox(self, equinox):
        check_astroyear(equinox, 'equinox', self._config, self._pos)
        self._equinox = equinox

    @equinox.deleter
    def equinox(self):
        self._equinox = None

    @property
    def epoch(self):
        """
        Specifies the epoch of the positions.  It must be a string
        specifying an astronomical year.
        """
        return self._epoch

    @epoch.setter
    def epoch(self, epoch):
        check_astroyear(epoch, 'epoch', self._config, self._pos)
        self._epoch = epoch

    @epoch.deleter
    def epoch(self):
        self._epoch = None


class FieldRef(SimpleElement, _UtypeProperty, _UcdProperty):
    """
    FIELDref_ element: used inside of GROUP_ elements to refer to remote FIELD_ elements.
    """
    _attr_list_11 = ['ref']
    _attr_list_12 = _attr_list_11 + ['ucd', 'utype']
    _element_name = "FIELDref"
    _utype_in_v1_2 = True
    _ucd_in_v1_2 = True

    def __init__(self, table, ref, ucd=None, utype=None, config=None, pos=None,
                 **extra):
        """
        *table* is the :class:`Table` object that this :class:`FieldRef`
        is a member of.

        *ref* is the ID to reference a :class:`Field` object defined
        elsewhere.
        """
        if config is None:
            config = {}
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

    @property
    def ref(self):
        """The ID_ of the FIELD_ that this FIELDref_ references."""
        return self._ref

    @ref.setter
    def ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref

    @ref.deleter
    def ref(self):
        self._ref = None

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


class ParamRef(SimpleElement, _UtypeProperty, _UcdProperty):
    """
    PARAMref_ element: used inside of GROUP_ elements to refer to remote PARAM_ elements.

    The keyword arguments correspond to setting members of the same
    name, documented below.

    It contains the following publicly-accessible members:

      *ref*: An XML ID referring to a <PARAM> element.
    """
    _attr_list_11 = ['ref']
    _attr_list_12 = _attr_list_11 + ['ucd', 'utype']
    _element_name = "PARAMref"
    _utype_in_v1_2 = True
    _ucd_in_v1_2 = True

    def __init__(self, table, ref, ucd=None, utype=None, config=None, pos=None):
        if config is None:
            config = {}

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

    @property
    def ref(self):
        """The ID_ of the PARAM_ that this PARAMref_ references."""
        return self._ref

    @ref.setter
    def ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref

    @ref.deleter
    def ref(self):
        self._ref = None

    def get_ref(self):
        """
        Lookup the :class:`Param` instance that this :class:``PARAMref``
        references.
        """
        for param in self._table._votable.iter_fields_and_params():
            if isinstance(param, Param) and param.ID == self.ref:
                return param
        vo_raise(
            "No params named '%s'" % self.ref,
            self._config, self._pos, KeyError)


class Group(Element, _IDProperty, _NameProperty, _UtypeProperty,
            _UcdProperty, _DescriptionProperty):
    """
    GROUP_ element: groups FIELD_ and PARAM_ elements.

    This information is currently ignored by the vo package---that is
    the columns in the recarray are always flat---but the grouping
    information is stored so that it can be written out again to the
    XML file.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """

    def __init__(self, table, ID=None, name=None, ref=None, ucd=None,
                 utype=None, id=None, config=None, pos=None, **extra):
        if config is None:
            config = {}
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

        self._entries = HomogeneousList(
            (FieldRef, ParamRef, Group, Param))

        warn_unknown_attrs('GROUP', six.iterkeys(extra), config, pos)

    def __repr__(self):
        return '<GROUP>... {0} entries ...</GROUP>'.format(len(self._entries))

    @property
    def ref(self):
        """
        Currently ignored, as it's not clear from the spec how this is
        meant to work.
        """
        return self._ref

    @ref.setter
    def ref(self, ref):
        xmlutil.check_id(ref, 'ref', self._config, self._pos)
        self._ref = ref

    @ref.deleter
    def ref(self):
        self._ref = None

    @property
    def entries(self):
        """
        [read-only] A list of members of the GROUP_.  This list may
        only contain objects of type :class:`Param`, :class:`Group`,
        :class:`ParamRef` and :class:`FieldRef`.
        """
        return self._entries

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
            attrib=w.object_attrs(
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


class Table(Element, _IDProperty, _NameProperty, _UcdProperty,
            _DescriptionProperty):
    """
    TABLE_ element: optionally contains data.

    It contains the following publicly-accessible and mutable
    attribute:

        *array*: A Numpy masked array of the data itself, where each
        row is a row of votable data, and columns are named and typed
        based on the <FIELD> elements of the table.  The mask is
        parallel to the data array, except for variable-length fields.
        For those fields, the numpy array's column type is "object"
        (``"O"``), and another masked array is stored there.

    If the Table contains no data, (for example, its enclosing
    :class:`Resource` has :attr:`~Resource.type` == 'meta') *array*
    will have zero-length.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, votable, ID=None, name=None, ref=None, ucd=None,
                 utype=None, nrows=None, id=None, config=None, pos=None,
                 **extra):
        if config is None:
            config = {}
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

        self._fields = HomogeneousList(Field)
        self._params = HomogeneousList(Param)
        self._groups = HomogeneousList(Group)
        self._links  = HomogeneousList(Link)
        self._infos  = HomogeneousList(Info)

        self.array = ma.array([])

        warn_unknown_attrs('TABLE', six.iterkeys(extra), config, pos)

    def __repr__(self):
        return repr(self.to_table())

    def __bytes__(self):
        return bytes(self.to_table())
    if six.PY2:
        __str__ = __bytes__

    def __unicode__(self):
        return six.text_type(self.to_table())
    if six.PY3:
        __str__ = __unicode__

    @property
    def ref(self):
        return self._ref

    @ref.setter
    def ref(self, ref):
        """
        Refer to another TABLE, previously defined, by the *ref* ID_
        for all metadata (FIELD_, PARAM_ etc.) information.
        """
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

    @ref.deleter
    def ref(self):
        self._ref = None

    @property
    def format(self):
        """
        [*required*] The serialization format of the table.  Must be
        one of:

          'tabledata' (TABLEDATA_), 'binary' (BINARY_), 'binary2' (BINARY2_)
          'fits' (FITS_).

        Note that the 'fits' format, since it requires an external
        file, can not be written out.  Any file read in with 'fits'
        format will be read out, by default, in 'tabledata' format.

        See :ref:`votable-serialization`.
        """
        return self._format

    @format.setter
    def format(self, format):
        format = format.lower()
        if format == 'fits':
            vo_raise("fits format can not be written out, only read.",
                     self._config, self._pos, NotImplementedError)
        if format == 'binary2':
            if not self._config['version_1_3_or_later']:
                vo_raise(
                    "binary2 only supported in votable 1.3 or later",
                    self._config, self._pos)
        elif format not in ('tabledata', 'binary'):
            vo_raise("Invalid format '%s'" % format,
                     self._config, self._pos)
        self._format = format

    @property
    def nrows(self):
        """
        [*immutable*] The number of rows in the table, as specified in
        the XML file.
        """
        return self._nrows

    @property
    def fields(self):
        """
        A list of :class:`Field` objects describing the types of each
        of the data columns.
        """
        return self._fields

    @property
    def params(self):
        """
        A list of parameters (constant-valued columns) for the
        table.  Must contain only :class:`Param` objects.
        """
        return self._params

    @property
    def groups(self):
        """
        A list of :class:`Group` objects describing how the columns
        and parameters are grouped.  Currently this information is
        only kept around for round-tripping and informational
        purposes.
        """
        return self._groups

    @property
    def links(self):
        """
        A list of :class:`Link` objects (pointers to other documents
        or servers through a URI) for the table.
        """
        return self._links

    @property
    def infos(self):
        """
        A list of :class:`Info` objects for the table.  Allows for
        post-operational diagnostics.
        """
        return self._infos

    def is_empty(self):
        """
        Returns True if this table doesn't contain any real data
        because it was skipped over by the parser (through use of the
        ``table_number`` kwarg).
        """
        return self._empty

    def create_arrays(self, nrows=0, config=None):
        """
        Create a new array to hold the data based on the current set
        of fields, and store them in the *array* and member variable.
        Any data in the existing array will be lost.

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

            dtype = []
            for x in fields:
                if six.PY3:
                    if x._unique_name == x.ID:
                        id = x.ID
                    else:
                        id = (x._unique_name, x.ID)
                else:
                    if x._unique_name == x.ID:
                        id = x.ID.encode('utf-8')
                    else:
                        id = (x._unique_name.encode('utf-8'),
                              x.ID.encode('utf-8'))
                dtype.append((id, x.converter.format))

            array = np.recarray((nrows,), dtype=np.dtype(dtype))
            descr_mask = []
            for d in array.dtype.descr:
                new_type = (d[1][1] == 'O' and 'O') or 'bool'
                if len(d) == 2:
                    descr_mask.append((d[0], new_type))
                elif len(d) == 3:
                    descr_mask.append((d[0], new_type, d[2]))
            mask = np.zeros((nrows,), dtype=descr_mask)

        self.array = ma.array(array, mask=mask)

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
        return int(np.ceil(size * RESIZE_AMOUNT))

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

    def _add_info(self, iterator, tag, data, config, pos):
        if not config.get('version_1_2_or_later'):
            warn_or_raise(W26, W26, ('INFO', 'TABLE', '1.2'), config, pos)
        info = Info(config=config, pos=pos, **data)
        self.infos.append(info)
        info.parse(iterator, config)

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

        table_id = config.get('table_id')
        if table_id is not None:
            if table_id != self.ID:
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
                            'DATA', six.iterkeys(data), config, pos)
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
                'INFO'        : self._add_info,
                'DESCRIPTION' : self._ignore_add}

            for start, tag, data, pos in iterator:
                if start:
                    if tag == 'DATA':
                        warn_unknown_attrs(
                            'DATA', six.iterkeys(data), config, pos)
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
            colnumbers = list(range(len(fields)))
        else:
            if isinstance(columns, six.string_types):
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
                            'TABLEDATA', six.iterkeys(data), config, pos)
                        self.array = self._parse_tabledata(
                            iterator, colnumbers, fields, config)
                        break
                    elif tag == 'BINARY':
                        warn_unknown_attrs(
                            'BINARY', six.iterkeys(data), config, pos)
                        self.array = self._parse_binary(
                            1, iterator, colnumbers, fields, config, pos)
                        break
                    elif tag == 'BINARY2':
                        if not config['version_1_3_or_later']:
                            warn_or_raise(
                                W52, W52, config['version'], config, pos)
                        self.array = self._parse_binary(
                            2, iterator, colnumbers, fields, config, pos)
                        break
                    elif tag == 'FITS':
                        warn_unknown_attrs(
                            'FITS', six.iterkeys(data), config, pos, ['extnum'])
                        try:
                            extnum = int(data.get('extnum', 0))
                            if extnum < 0:
                                raise ValueError()
                        except ValueError:
                            vo_raise(E17, (), config, pos)
                        self.array = self._parse_fits(
                            iterator, extnum, config)
                        break
                    else:
                        warn_or_raise(W37, W37, tag, config, pos)
                        break

        for start, tag, data, pos in iterator:
            if not start and tag == 'DATA':
                break

        for start, tag, data, pos in iterator:
            if start and tag == 'INFO':
                if not config.get('version_1_2_or_later'):
                    warn_or_raise(
                        W26, W26, ('INFO', 'TABLE', '1.2'), config, pos)
                info = Info(config=config, pos=pos, **data)
                self.infos.append(info)
                info.parse(iterator, config)
            elif not start and tag == 'TABLE':
                break

        return self

    def _parse_tabledata(self, iterator, colnumbers, fields, config):
        # Since we don't know the number of rows up front, we'll
        # reallocate the record array to make room as we go.  This
        # prevents the need to scan through the XML twice.  The
        # allocation is by factors of 1.5.
        invalid = config.get('invalid', 'exception')

        # Need to have only one reference so that we can resize the
        # array
        array = self.array
        del self.array

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
                            tag, six.iterkeys(data), config, pos, ['encoding'])
                    else:
                        if tag == 'TD':
                            if i >= len(fields):
                                vo_raise(E20, len(fields), config, pos)

                            if colnumbers_bits[i]:
                                try:
                                    if binary:
                                        rawdata = base64.b64decode(
                                            data.encode('ascii'))
                                        buf = io.BytesIO(rawdata)
                                        buf.seek(0)
                                        try:
                                            value, mask_value = binparsers[i](
                                                buf.read)
                                        except Exception as e:
                                            vo_reraise(
                                                e, config, pos,
                                                "(in row %d, col '%s')" %
                                                (len(array_chunk),
                                                 fields[i].ID))
                                    else:
                                        try:
                                            value, mask_value = parsers[i](
                                                data, config, pos)
                                        except Exception as e:
                                            vo_reraise(
                                                e, config, pos,
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
                        array = _resize(array, alloc_rows)
                    array[numrows:numrows + chunk_size] = array_chunk
                    array.mask[numrows:numrows + chunk_size] = mask_chunk
                    numrows += chunk_size
                    array_chunk = []
                    mask_chunk = []

            elif not start and tag == 'TABLEDATA':
                break

        # Now, resize the array to the exact number of rows we need and
        # put the last chunk values in there.
        alloc_rows = numrows + len(array_chunk)

        array = _resize(array, alloc_rows)
        array[numrows:] = array_chunk
        if alloc_rows != 0:
            array.mask[numrows:] = mask_chunk
        numrows += len(array_chunk)

        if (self.nrows is not None and
            self.nrows >= 0 and
            self.nrows != numrows):
            warn_or_raise(W18, W18, (self.nrows, numrows), config, pos)
        self._nrows = numrows

        return array

    def _get_binary_data_stream(self, iterator, config):
        have_local_stream = False
        for start, tag, data, pos in iterator:
            if tag == 'STREAM':
                if start:
                    warn_unknown_attrs(
                        'STREAM', six.iterkeys(data), config, pos,
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
            fd = urllib.request.urlopen(href)
            if encoding is not None:
                if encoding == 'gzip':
                    from ...utils.compat import gzip
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

        return careful_read

    def _parse_binary(self, mode, iterator, colnumbers, fields, config, pos):
        fields = self.fields

        careful_read = self._get_binary_data_stream(iterator, config)

        # Need to have only one reference so that we can resize the
        # array
        array = self.array
        del self.array

        binparsers = [field.converter.binparse for field in fields]

        numrows = 0
        alloc_rows = len(array)
        while True:
            # Resize result arrays if necessary
            if numrows >= alloc_rows:
                alloc_rows = self._resize_strategy(alloc_rows)
                array = _resize(array, alloc_rows)

            row_data = []
            row_mask_data = []

            try:
                if mode == 2:
                    mask_bits = careful_read(int((len(fields) + 7) / 8))
                    row_mask_data = list(converters.bitarray_to_bool(
                        mask_bits, len(fields)))
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
                    if mode == 1:
                        row_mask_data.append(value_mask)
                    else:
                        row_mask_data[i] = row_mask_data[i] or value_mask
            except EOFError:
                break

            row = [x.converter.default for x in fields]
            row_mask = [False] * len(fields)
            for i in colnumbers:
                row[i] = row_data[i]
                row_mask[i] = row_mask_data[i]

            array[numrows] = tuple(row)
            array.mask[numrows] = tuple(row_mask)
            numrows += 1

        array = _resize(array, numrows)

        return array

    def _parse_fits(self, iterator, extnum, config):
        for start, tag, data, pos in iterator:
            if tag == 'STREAM':
                if start:
                    warn_unknown_attrs(
                        'STREAM', six.iterkeys(data), config, pos,
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

        fd = urllib.request.urlopen(href)
        if encoding is not None:
            if encoding == 'gzip':
                from ...utils.compat import gzip
                fd = gzip.GzipFile(href, 'r', fileobj=fd)
            elif encoding == 'base64':
                fd = codecs.EncodedFile(fd, 'base64')
            else:
                vo_raise(
                    "Unknown encoding type '%s'" % encoding,
                    self._config, self._pos, NotImplementedError)

        hdulist = fits.open(fd)

        array = hdulist[int(extnum)].data
        if array.dtype != self.array.dtype:
            warn_or_raise(W19, W19, (), self._config, self._pos)

        return array

    def to_xml(self, w, **kwargs):
        specified_format = kwargs.get('tabledata_format')
        if specified_format is not None:
            format = specified_format
        else:
            format = self.format
        if format == 'fits':
            format = 'tabledata'

        with w.tag(
            'TABLE',
            attrib=w.object_attrs(
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
            elif kwargs['version_1_2_or_later']:
                index = list(self._votable.iter_tables()).index(self)
                group = Group(self, ID="_g{0}".format(index))
                group.to_xml(w, **kwargs)

            if len(self.array):
                with w.tag('DATA'):
                    if format == 'tabledata':
                        self._write_tabledata(w, **kwargs)
                    elif format == 'binary':
                        self._write_binary(1, w, **kwargs)
                    elif format == 'binary2':
                        self._write_binary(2, w, **kwargs)

            if kwargs['version_1_2_or_later']:
                for element in self._infos:
                    element.to_xml(w, **kwargs)

    def _write_tabledata(self, w, **kwargs):
        fields = self.fields
        array = self.array

        with w.tag('TABLEDATA'):
            w._flush()
            if (_has_c_tabledata_writer and
                not kwargs.get('_debug_python_based_parser')):
                supports_empty_values = [
                    field.converter.supports_empty_values(kwargs)
                    for field in fields]
                fields = [field.converter.output for field in fields]
                indent = len(w._tags) - 1
                tablewriter.write_tabledata(
                    w.write, array.data, array.mask, fields,
                    supports_empty_values, indent, 1 << 8)
            else:
                write = w.write
                indent_spaces = w.get_indentation_spaces()
                tr_start = indent_spaces + "<TR>\n"
                tr_end = indent_spaces + "</TR>\n"
                td = indent_spaces + " <TD>%s</TD>\n"
                td_empty = indent_spaces + " <TD/>\n"
                fields = [(i, field.converter.output,
                           field.converter.supports_empty_values(kwargs))
                          for i, field in enumerate(fields)]
                for row in xrange(len(array)):
                    write(tr_start)
                    array_row = array.data[row]
                    mask_row = array.mask[row]
                    for i, output, supports_empty_values in fields:
                        data = array_row[i]
                        masked = mask_row[i]
                        if supports_empty_values and np.all(masked):
                            write(td_empty)
                        else:
                            try:
                                val = output(data, masked)
                            except Exception as e:
                                vo_reraise(e,
                                           additional="(in row %d, col '%s')" %
                                           (row, self.fields[i].ID))
                            if len(val):
                                write(td % val)
                            else:
                                write(td_empty)
                    write(tr_end)

    def _write_binary(self, mode, w, **kwargs):
        fields = self.fields
        array = self.array
        if mode == 1:
            tag_name = 'BINARY'
        else:
            tag_name = 'BINARY2'

        with w.tag(tag_name):
            with w.tag('STREAM', encoding='base64'):
                fields_basic = [(i, field.converter.binoutput)
                                for (i, field) in enumerate(fields)]

                data = io.BytesIO()
                for row in xrange(len(array)):
                    array_row = array.data[row]
                    array_mask = array.mask[row]

                    if mode == 2:
                        flattened = np.array([np.all(x) for x in array_mask])
                        data.write(converters.bool_to_bitarray(flattened))

                    for i, converter in fields_basic:
                        try:
                            chunk = converter(array_row[i], array_mask[i])
                            assert type(chunk) == type(b'')
                        except Exception as e:
                            vo_reraise(e,
                                       additional="(in row %d, col '%s')" %
                                       (row, fields[i].ID))
                        data.write(chunk)

                w._flush()
                w.write(base64.b64encode(data.getvalue()).decode('ascii'))

    def to_table(self, use_names_over_ids=False):
        """
        Convert this VO Table to an `astropy.table.Table` instance.

        Parameters
        ----------
        use_names_over_ids : bool, optional
           When `True` use the ``name`` attributes of columns as the
           names of columns in the `astropy.table.Table` instance.
           Since names are not guaranteed to be unique, this may cause
           some columns to be renamed by appending numbers to the end.
           Otherwise (default), use the ID attributes as the column
           names.

        .. warning::
           Variable-length array fields may not be restored
           identically when round-tripping through the
           `astropy.table.Table` instance.
        """
        from ...table import Table

        meta = {}
        for key in ['ID', 'name', 'ref', 'ucd', 'utype', 'description']:
            val = getattr(self, key, None)
            if val is not None:
                meta[key] = val

        if use_names_over_ids:
            names = [field.name for field in self.fields]
            unique_names = []
            for i, name in enumerate(names):
                new_name = name
                i = 2
                while new_name in unique_names:
                    new_name = '{0}{1}'.format(name, i)
                    i += 1
                if six.PY2:
                    new_name = new_name.encode(
                        sys.getdefaultencoding(), 'replace')
                unique_names.append(new_name)
            array = self.array.copy()
            array.dtype.names = unique_names
            names = unique_names
        else:
            array = self.array
            names = [field.ID for field in self.fields]

        table = Table(self.array, meta=meta)

        for name, field in zip(names, self.fields):
            column = table[name]
            field.to_table_column(column)

        return table

    @classmethod
    def from_table(cls, votable, table):
        """
        Create a `Table` instance from a given `astropy.table.Table`
        instance.
        """
        kwargs = {}
        for key in ['ID', 'name', 'ref', 'ucd', 'utype']:
            val = table.meta.get(key)
            if val is not None:
                kwargs[key] = val
        new_table = cls(votable, **kwargs)
        if 'description' in table.meta:
            new_table.description = table.meta['description']

        for colname in table.colnames:
            column = table[colname]
            new_table.fields.append(Field.from_table_column(votable, column))

        if table.mask is None:
            new_table.array = ma.array(np.asarray(table))
        else:
            new_table.array = ma.array(np.asarray(table),
                                       mask=np.asarray(table.mask))

        return new_table

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

    get_field_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_fields_and_params', 'FIELD or PARAM',
        """
        Looks up a FIELD or PARAM element by the given ID.
        """)

    get_field_by_id_or_name = _lookup_by_id_or_name_factory(
        'iter_fields_and_params', 'FIELD or PARAM',
        """
        Looks up a FIELD or PARAM element by the given ID or name.
        """)

    get_fields_by_utype = _lookup_by_attr_factory(
        'utype', False, 'iter_fields_and_params', 'FIELD or PARAM',
        """
        Looks up a FIELD or PARAM element by the given utype and
        returns an iterator emitting all matches.
        """)

    def iter_groups(self):
        """
        Recursively iterate over all GROUP elements in the TABLE.
        """
        for group in self.groups:
            yield group
            for g in group.iter_groups():
                yield g

    get_group_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_groups', 'GROUP',
        """
        Looks up a GROUP element by the given ID.  Used by the group's
        "ref" attribute
        """)

    get_groups_by_utype = _lookup_by_attr_factory(
        'utype', False, 'iter_groups', 'GROUP',
        """
        Looks up a GROUP element by the given utype and returns an
        iterator emitting all matches.
        """)

    def iter_info(self):
        for info in self.infos:
            yield info


class Resource(Element, _IDProperty, _NameProperty, _UtypeProperty,
               _DescriptionProperty):
    """
    RESOURCE_ element: Groups TABLE_ and RESOURCE_ elements.

    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, name=None, ID=None, utype=None, type='results',
                 id=None, config=None, pos=None, **kwargs):
        if config is None:
            config = {}
        self._config           = config
        self._pos              = pos

        Element.__init__(self)
        self.name              = name
        self.ID                = resolve_id(ID, id, config, pos)
        self.utype             = utype
        self.type              = type
        self._extra_attributes = kwargs
        self.description       = None

        self._coordinate_systems = HomogeneousList(CooSys)
        self._params             = HomogeneousList(Param)
        self._infos              = HomogeneousList(Info)
        self._links              = HomogeneousList(Link)
        self._tables             = HomogeneousList(Table)
        self._resources          = HomogeneousList(Resource)

        warn_unknown_attrs('RESOURCE', six.iterkeys(kwargs), config, pos)

    def __repr__(self):
        buff = io.StringIO()
        XMLWriter(buff).element(
            self._element_name,
            attrib=w.object_attrs(self, self._attr_list))
        return buff.getvalue()

    @property
    def type(self):
        """
        [*required*] The type of the resource.  Must be either:

          - 'results': This resource contains actual result values
            (default)

          - 'meta': This resource contains only datatype descriptions
            (FIELD_ elements), but no actual data.
        """
        return self._type

    @type.setter
    def type(self, type):
        if type not in ('results', 'meta'):
            vo_raise(E18, type, self._config, self._pos)
        self._type = type

    @property
    def extra_attributes(self):
        """
        A dictionary of string keys to string values containing any
        extra attributes of the RESOURCE_ element that are not defined
        in the specification.  (The specification explicitly allows
        for extra attributes here, but nowhere else.)
        """
        return self._extra_attributes

    @property
    def coordinate_systems(self):
        """
        A list of coordinate system definitions (COOSYS_ elements) for
        the RESOURCE_.  Must contain only `CooSys` objects.
        """
        return self._coordinate_systems

    @property
    def infos(self):
        """
        A list of informational parameters (key-value pairs) for the
        resource.  Must only contain `Info` objects.
        """
        return self._infos

    @property
    def params(self):
        """
        A list of parameters (constant-valued columns) for the
        resource.  Must contain only `Param` objects.
        """
        return self._params

    @property
    def links(self):
        """
        A list of links (pointers to other documents or servers
        through a URI) for the resource.  Must contain only `Link`
        objects.
        """
        return self._links

    @property
    def tables(self):
        """
        A list of tables in the resource.  Must contain only
        `Table` objects.
        """
        return self._tables

    @property
    def resources(self):
        """
        A list of nested resources inside this resource.  Must contain
        only `Resource` objects.
        """
        return self._resources

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
            elif tag == 'RESOURCE':
                break

        del self._votable

        return self

    def to_xml(self, w, **kwargs):
        attrs = w.object_attrs(self, ('ID', 'type', 'utype'))
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

    def iter_info(self):
        """
        Recursively iterates over all the INFO_ elements in the
        resource and nested resources.
        """
        for info in self.infos:
            yield info
        for table in self.tables:
            for info in table.iter_info():
                yield info
        for resource in self.resources:
            for info in resource.iter_info():
                yield info


class VOTableFile(Element, _IDProperty, _DescriptionProperty):
    """
    VOTABLE_ element: represents an entire file.

    The keyword arguments correspond to setting members of the same
    name, documented below.

    *version* is settable at construction time only, since conformance
    tests for building the rest of the structure depend on it.
    """

    def __init__(self, ID=None, id=None, config=None, pos=None, version="1.2"):
        if config is None:
            config = {}
        self._config             = config
        self._pos                = pos

        Element.__init__(self)
        self.ID                  = resolve_id(ID, id, config, pos)
        self.description         = None

        self._coordinate_systems = HomogeneousList(CooSys)
        self._params             = HomogeneousList(Param)
        self._infos              = HomogeneousList(Info)
        self._resources          = HomogeneousList(Resource)
        self._groups             = HomogeneousList(Group)

        version = str(version)
        assert version in ("1.0", "1.1", "1.2")
        self._version            = version

    def __repr__(self):
        n_tables = len(list(self.iter_tables()))
        return '<VOTABLE>... {0} tables ...</VOTABLE>'.format(n_tables)

    @property
    def version(self):
        """
        The version of the VOTable specification that the file uses.
        """
        return self._version

    @version.setter
    def version(self, version):
        version = str(version)
        if version not in ('1.1', '1.2', '1.3'):
            raise ValueError(
                "astropy.io.votable only supports VOTable versions "
                "1.1, 1.2 and 1.3")
        self._version = version

    @property
    def coordinate_systems(self):
        """
        A list of coordinate system descriptions for the file.  Must
        contain only `CooSys` objects.
        """
        return self._coordinate_systems

    @property
    def params(self):
        """
        A list of parameters (constant-valued columns) that apply to
        the entire file.  Must contain only `Param` objects.
        """
        return self._params

    @property
    def infos(self):
        """
        A list of informational parameters (key-value pairs) for the
        entire file.  Must only contain `Info` objects.
        """
        return self._infos

    @property
    def resources(self):
        """
        A list of resources, in the order they appear in the file.
        Must only contain `Resource` objects.
        """
        return self._resources

    @property
    def groups(self):
        """
        A list of groups, in the order they appear in the file.  Only
        supported as a child of the VOTABLE element in VOTable 1.2 or
        later.
        """
        return self._groups

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
                if tag == 'xml':
                    pass
                elif tag == 'VOTABLE':
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
                        if config['version'] not in ('1.1', '1.2', '1.3'):
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
        config['version_1_3_or_later'] = \
            util.version_compare(config['version'], '1.3') >= 0

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

        if not len(self.resources) and config['version_1_2_or_later']:
            warn_or_raise(W53, W53, (), config, pos)

        return self

    def to_xml(self, fd, write_null_values=False,
               compressed=False, tabledata_format=None,
               _debug_python_based_parser=False,
               _astropy_version=None):
        """
        Write to an XML file.

        Parameters
        ----------
        fd : str path or writable file-like object
            Where to write the file.

        write_null_values : bool, optional
            Deprecated and retained for backward compatibility.  When
            ``write_null_values`` was `False`, invalid VOTable files
            could be generated, so the option has just been removed
            entirely.

        compressed : bool, optional
            When `True`, write to a gzip-compressed file.  (Default:
            `False`)

        tabledata_format : str, optional
            Override the format of the table(s) data to write.  Must
            be one of ``tabledata`` (text representation), ``binary`` or
            ``binary2``.  By default, use the format that was specified
            in each `Table` object as it was created or read in.  See
            :ref:`votable-serialization`.
        """
        if write_null_values != False:
            warnings.warn(
                "write_null_values has been deprecated and has no effect",
                AstropyDeprecationWarning)

        if tabledata_format is not None:
            if tabledata_format.lower() not in (
                    'tabledata', 'binary', 'binary2'):
                raise ValueError("Unknown format type '{0}'".format(format))

        kwargs = {
            'version': self.version,
            'version_1_1_or_later':
                util.version_compare(self.version, '1.1') >= 0,
            'version_1_2_or_later':
                util.version_compare(self.version, '1.2') >= 0,
            'version_1_3_or_later':
                util.version_compare(self.version, '1.3') >= 0,
            'tabledata_format':
                tabledata_format,
            '_debug_python_based_parser': _debug_python_based_parser,
            '_group_number': 1}

        with util.convert_to_writable_filelike(
            fd, compressed=compressed) as fd:
            w = XMLWriter(fd)
            version = self.version
            if _astropy_version is None:
                lib_version = astropy_version
            else:
                lib_version = _astropy_version

            xml_header = """
<?xml version="1.0" encoding="utf-8"?>
<!-- Produced with astropy.io.votable version %(lib_version)s
     http://www.astropy.org/ -->\n"""
            w.write(xml_header.lstrip() % locals())

            with w.tag('VOTABLE',
                       {'version': version,
                        'xmlns:xsi':
                            "http://www.w3.org/2001/XMLSchema-instance",
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

    get_table_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_tables', 'TABLE',
        """
        Looks up a TABLE_ element by the given ID.  Used by the table
        "ref" attribute.
        """)

    get_tables_by_utype = _lookup_by_attr_factory(
        'utype', False, 'iter_tables', 'TABLE',
        """
        Looks up a TABLE_ element by the given utype, and returns an
        iterator emitting all matches.
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
        Recursively iterate over all FIELD_ and PARAM_ elements in the
        VOTABLE_ file.
        """
        for resource in self.resources:
            for field in resource.iter_fields_and_params():
                yield field

    get_field_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_fields_and_params', 'FIELD',
        """
        Looks up a FIELD_ element by the given ID_.  Used by the field's
        "ref" attribute.
        """)

    get_fields_by_utype = _lookup_by_attr_factory(
        'utype', False, 'iter_fields_and_params', 'FIELD',
        """
        Looks up a FIELD_ element by the given utype and returns an
        iterator emitting all matches.
        """)

    get_field_by_id_or_name = _lookup_by_id_or_name_factory(
        'iter_fields_and_params', 'FIELD',
        """
        Looks up a FIELD_ element by the given ID_ or name.
        """)

    def iter_values(self):
        """
        Recursively iterate over all VALUES_ elements in the VOTABLE_
        file.
        """
        for field in self.iter_fields_and_params():
            yield field.values

    get_values_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_values', 'VALUES',
        """
        Looks up a VALUES_ element by the given ID.  Used by the values
        "ref" attribute.
        """)

    def iter_groups(self):
        """
        Recursively iterate over all GROUP_ elements in the VOTABLE_
        file.
        """
        for table in self.iter_tables():
            for group in table.iter_groups():
                yield group

    get_group_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_groups', 'GROUP',
        """
        Looks up a GROUP_ element by the given ID.  Used by the group's
        "ref" attribute
        """)

    get_groups_by_utype = _lookup_by_attr_factory(
        'utype', False, 'iter_groups', 'GROUP',
        """
        Looks up a GROUP_ element by the given utype and returns an
        iterator emitting all matches.
        """)

    def iter_coosys(self):
        """
        Recursively iterate over all COOSYS_ elements in the VOTABLE_
        file.
        """
        for coosys in self.coordinate_systems:
            yield coosys
        for resource in self.resources:
            for coosys in resource.iter_coosys():
                yield coosys

    get_coosys_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_coosys', 'COOSYS',
        """Looks up a COOSYS_ element by the given ID.""")

    def iter_info(self):
        """
        Recursively iterate over all INFO_ elements in the VOTABLE_
        file.
        """
        for info in self.infos:
            yield info
        for resource in self.resources:
            for info in resource.iter_info():
                yield info

    get_info_by_id = _lookup_by_attr_factory(
        'ID', True, 'iter_info', 'INFO',
        """Looks up a INFO element by the given ID.""")

    def set_all_tables_format(self, format):
        """
        Set the output storage format of all tables in the file.
        """
        for table in self.iter_tables():
            table.format = format

    @classmethod
    def from_table(cls, table, table_id=None):
        """
        Create a `VOTableFile` instance from a given
        `astropy.table.Table` instance.

        Parameters
        ----------
        table_id : str, optional
            Set the given ID attribute on the returned Table instance.
        """
        votable_file = cls()
        resource = Resource()
        votable = Table.from_table(votable_file, table)
        if table_id is not None:
            votable.ID = table_id
        resource.tables.append(votable)
        votable_file.resources.append(resource)
        return votable_file
