# Licensed under a 3-clause BSD style license - see PYFITS.rst

import copy
import re
import string
import sys
import warnings

import numpy as np

from .util import _str_to_num, _is_int, maketrans, translate, _words_group
from .verify import _Verify, _ErrList, VerifyError, VerifyWarning

from . import ENABLE_RECORD_VALUED_KEYWORD_CARDS, STRIP_HEADER_WHITESPACE
from ...utils import deprecated


__all__ = ['Card', 'CardList', 'create_card', 'create_card_from_string',
           'upper_key', 'createCard', 'createCardFromString', 'upperKey',
           'Undefined']


FIX_FP_TABLE = maketrans('de', 'DE')
FIX_FP_TABLE2 = maketrans('dD', 'eE')


BLANK_CARD = ' ' * 80


class Undefined:
    """Undefined value."""

    def __init__(self):
        # This __init__ is required to be here for Sphinx documentation
        pass
UNDEFINED = Undefined()


class CardList(list):
    """
    .. deprecated:: 3.1
        `CardList` used to provide the list-like functionality for manipulating
        a header as a list of cards.  This functionality is now subsumed into
        the `Header` class itself, so it is no longer necessary to create or
        use `CardList`\s.
    """

    def __init__(self, cards=[], keylist=None):
        """
        Construct the `CardList` object from a list of `Card` objects.

        `CardList` is now merely a thin wrapper around `Header` to provide
        backwards compatibility for the old API.  This should not be used for
        any new code.

        Parameters
        ----------
        cards
            A list of `Card` objects.
        """

        warnings.warn(
                'The CardList class has been deprecated; all its former '
                'functionality has been subsumed by the Header class, so '
                'CardList objects should not be directly created.  See the '
                'PyFITS 3.1.0 CHANGELOG for more details.',
                DeprecationWarning)

        # This is necessary for now to prevent a circular import
        from .header import Header

        # I'm not sure if they keylist argument here was ever really useful;
        # I'm going to just say don't use it.
        if keylist is not None:
            raise ValueError(
                'The keylist argument to CardList() is no longer supported.')

        if isinstance(cards, Header):
            self._header = cards
        else:
            self._header = Header(cards)

        super(CardList, self).__init__(self._header.cards)

    def __contains__(self, key):
        return key in self._header

    def __iter__(self):
        return iter(self._header.cards)

    def __getitem__(self, key):
        """Get a `Card` by indexing or by the keyword name."""

        if self._header._haswildcard(key):
            return [copy.copy(self._header._cards[idx])
                    for idx in self._header._wildcardmatch(key)]
        elif isinstance(key, slice):
            return CardList(self._header.cards[key])

        idx = self._header._cardindex(key)
        return self._header.cards[idx]

    def __setitem__(self, key, value):
        """Set a `Card` by indexing or by the keyword name."""

        if isinstance(value, tuple) and (1 < len(value) <= 3):
            value = Card(*value)

        if isinstance(value, Card):
            idx = self._header._cardindex(key)
            card = self._header.cards[idx]
            if str(card) != str(value):
                # Replace the existing card at this index by delete/insert
                del self._header[idx]
                self._header.insert(idx, value)
        else:
            raise ValueError('%s is not a Card' % str(value))

    def __delitem__(self, key):
        """Delete a `Card` from the `CardList`."""

        if key not in self._header._keyword_indices:
            raise KeyError("Keyword '%s' not found." % key)
        del self._header[key]

    def __getslice__(self, start, end):
        return CardList(self[slice(start, end)])

    def __repr__(self):
        """Format a list of cards into a string."""

        return str(self._header)

    def __str__(self):
        """Format a list of cards into a printable string."""

        return '\n'.join(str(card) for card in self)

    @deprecated('3.1', alternative=':meth:`Header.copy`', pending=False)
    def copy(self):
        """Make a (deep)copy of the `CardList`."""

        return CardList(self._header.copy())

    @deprecated('3.1', alternative=':meth:`Header.keys`', pending=False)
    def keys(self):
        """
        Return a list of all keywords from the `CardList`.
        """

        return self._header.keys()

    @deprecated('3.1', alternative=':meth:`Header.values`', pending=False)
    def values(self):
        """
        Return a list of the values of all cards in the `CardList`.

        For `RecordValuedKeywordCard` objects, the value returned is
        the floating point value, exclusive of the
        ``field_specifier``.
        """

        return self._header.values()

    @deprecated('3.1', alternative=':meth:`Header.append`', pending=False)
    def append(self, card, useblanks=True, bottom=False):
        """
        Append a `Card` to the `CardList`.

        Parameters
        ----------
        card : `Card` object
            The `Card` to be appended.

        useblanks : bool, optional
            Use any *extra* blank cards?

            If `useblanks` is `True`, and if there are blank cards
            directly before ``END``, it will use this space first,
            instead of appending after these blank cards, so the total
            space will not increase.  When `useblanks` is `False`, the
            card will be appended at the end, even if there are blank
            cards in front of ``END``.

        bottom : bool, optional
           If `False` the card will be appended after the last
           non-commentary card.  If `True` the card will be appended
           after the last non-blank card.
        """

        self._header.append(card, useblanks=useblanks, bottom=bottom)

    @deprecated('3.1', alternative=':meth:`Header.extend`', pending=False)
    def extend(self, cards):
        self._header.extend(cards)

    @deprecated('3.1', alternative=':meth:`Header.insert`', pending=False)
    def insert(self, idx, card, useblanks=True):
        """
        Insert a `Card` to the `CardList`.

        Parameters
        ----------
        pos : int
            The position (index, keyword name will not be allowed) to
            insert. The new card will be inserted before it.

        card : `Card` object
            The card to be inserted.

        useblanks : bool, optional
            If `useblanks` is `True`, and if there are blank cards
            directly before ``END``, it will use this space first,
            instead of appending after these blank cards, so the total
            space will not increase.  When `useblanks` is `False`, the
            card will be appended at the end, even if there are blank
            cards in front of ``END``.
        """

        self._header.insert(idx, card, useblanks=useblanks)

    @deprecated('3.1', alternative=':meth:`Header.remove`')
    def remove(self, card):
        del self._header[self.index(card)]

    @deprecated('3.1', alternative=':meth:`Header.pop`')
    def pop(self, index=-1):
        return self._header.pop(index)

    @deprecated('3.1', alternative=':meth:`Header.index`')
    def index(self, card):
        return self._header._cards.index(card)

    @deprecated('3.1', alternative=':meth:`Header.count`')
    def count(self, card):
        return self._header._cards.count(card)

    @deprecated('3.1', alternative=':meth:`Header.index`', pending=False)
    def index_of(self, key, backward=False):
        """
        Get the index of a keyword in the `CardList`.

        Parameters
        ----------
        key : str or int
            The keyword name (a string) or the index (an integer).

        backward : bool, (optional)
            When `True`, search the index from the ``END``, i.e.,
            backward.

        Returns
        -------
        index : int
            The index of the `Card` with the given keyword.
        """

        # Backward is just ignored now, since the search is not linear anyways

        if _is_int(key) or isinstance(key, basestring):
            return self._header._cardindex(key)
        else:
            raise KeyError('Illegal key data type %s' % type(key))

    @deprecated('3.1', alternative='`header[<wildcard_pattern>]`')
    def filter_list(self, key):
        """
        Construct a `CardList` that contains references to all of the cards in
        this `CardList` that match the input key value including any special
        filter keys (``*``, ``?``, and ``...``).

        Parameters
        ----------
        key : str
            key value to filter the list with

        Returns
        -------
        cardlist :
            A `CardList` object containing references to all the
            requested cards.
        """

        return CardList(self._header[key])

    # For API backwards-compatibility
    @deprecated('3.0', alternative=':meth:`filter_list`', pending=False)
    def filterList(self, key):
        return self.filter_list(key)

    @deprecated('3.1', pending=False)
    def count_blanks(self):
        """
        Returns how many blank cards are *directly* before the ``END``
        card.
        """

        return self._header._countblanks()


class Card(_Verify):
    length = 80

    # String for a FITS standard compliant (FSC) keyword.
    _keywd_FSC_RE = re.compile(r'^[A-Z0-9_-]{0,8}$')
    # This will match any printable ASCII character excluding '='
    _keywd_hierarch_RE = re.compile(r'^(?:HIERARCH )?(?:^[ -<>-~]+ ?)+$', re.I)

    # A number sub-string, either an integer or a float in fixed or
    # scientific notation.  One for FSC and one for non-FSC (NFSC) format:
    # NFSC allows lower case of DE for exponent, allows space between sign,
    # digits, exponent sign, and exponents
    _digits_FSC = r'(\.\d+|\d+(\.\d*)?)([DE][+-]?\d+)?'
    _digits_NFSC = r'(\.\d+|\d+(\.\d*)?) *([deDE] *[+-]? *\d+)?'
    _numr_FSC = r'[+-]?' + _digits_FSC
    _numr_NFSC = r'[+-]? *' + _digits_NFSC

    # This regex helps delete leading zeros from numbers, otherwise
    # Python might evaluate them as octal values (this is not-greedy, however,
    # so it may not strip leading zeros from a float, which is fine)
    _number_FSC_RE = re.compile(r'(?P<sign>[+-])?0*?(?P<digt>%s)'
                                % _digits_FSC)
    _number_NFSC_RE = re.compile(r'(?P<sign>[+-])? *0*?(?P<digt>%s)'
                                 % _digits_NFSC)

    # FSC commentary card string which must contain printable ASCII characters.
    _ascii_text = r'[ -~]*$'
    _comment_FSC_RE = re.compile(_ascii_text)

    # Checks for a valid value/comment string.  It returns a match object
    # for a valid value/comment string.
    # The valu group will return a match if a FITS string, boolean,
    # number, or complex value is found, otherwise it will return
    # None, meaning the keyword is undefined.  The comment field will
    # return a match if the comment separator is found, though the
    # comment maybe an empty string.
    _value_FSC_RE = re.compile(
        r'(?P<valu_field> *'
            r'(?P<valu>'

                #  The <strg> regex is not correct for all cases, but
                #  it comes pretty darn close.  It appears to find the
                #  end of a string rather well, but will accept
                #  strings with an odd number of single quotes,
                #  instead of issuing an error.  The FITS standard
                #  appears vague on this issue and only states that a
                #  string should not end with two single quotes,
                #  whereas it should not end with an even number of
                #  quotes to be precise.
                #
                #  Note that a non-greedy match is done for a string,
                #  since a greedy match will find a single-quote after
                #  the comment separator resulting in an incorrect
                #  match.
                r'\'(?P<strg>([ -~]+?|\'\'|)) *?\'(?=$|/| )|'
                r'(?P<bool>[FT])|'
                r'(?P<numr>' + _numr_FSC + r')|'
                r'(?P<cplx>\( *'
                    r'(?P<real>' + _numr_FSC + r') *, *'
                    r'(?P<imag>' + _numr_FSC + r') *\))'
            r')? *)'
        r'(?P<comm_field>'
            r'(?P<sepr>/ *)'
            r'(?P<comm>[!-~][ -~]*)?'
        r')?$')

    _value_NFSC_RE = re.compile(
        r'(?P<valu_field> *'
            r'(?P<valu>'
                r'\'(?P<strg>([ -~]+?|\'\'|) *?)\'(?=$|/| )|'
                r'(?P<bool>[FT])|'
                r'(?P<numr>' + _numr_NFSC + r')|'
                r'(?P<cplx>\( *'
                    r'(?P<real>' + _numr_NFSC + r') *, *'
                    r'(?P<imag>' + _numr_NFSC + r') *\))'
            r')? *)'
        r'(?P<comm_field>'
            r'(?P<sepr>/ *)'
            r'(?P<comm>(.|\n)*)'
        r')?$')

    _rvkc_identifier = r'[a-zA-Z_]\w*'
    _rvkc_field = _rvkc_identifier + r'(\.\d+)?'
    _rvkc_field_specifier_s = r'%s(\.%s)*' % ((_rvkc_field,) * 2)
    _rvkc_field_specifier_val = (r'(?P<keyword>%s): (?P<val>%s)' %
                                 (_rvkc_field_specifier_s, _numr_FSC))
    _rvkc_keyword_val = r'\'(?P<rawval>%s)\'' % _rvkc_field_specifier_val
    _rvkc_keyword_val_comm = (r' +%s *(/ *(?P<comm>[ -~]*))?$' %
                              _rvkc_keyword_val)

    _rvkc_field_specifier_val_RE = re.compile(_rvkc_field_specifier_val + '$')

    # regular expression to extract the key and the field specifier from a
    # string that is being used to index into a card list that contains
    # record value keyword cards (ex. 'DP1.AXIS.1')
    _rvkc_keyword_name_RE = (
        re.compile(r'(?P<keyword>%s)\.(?P<field_specifier>%s)$' %
                   (_rvkc_identifier, _rvkc_field_specifier_s)))

    # regular expression to extract the field specifier and value and comment
    # from the string value of a record value keyword card
    # (ex "'AXIS.1: 1' / a comment")
    _rvkc_keyword_val_comm_RE = re.compile(_rvkc_keyword_val_comm)

    _commentary_keywords = ['', 'COMMENT', 'HISTORY', 'END']

    def __init__(self, keyword=None, value=None, comment=None, **kwargs):
        # For backwards compatibility, support the 'key' keyword argument:
        if keyword is None and 'key' in kwargs:
            keyword = kwargs['key']

        self._keyword = None
        self._value = None
        self._comment = None
        self._rawvalue = None
        self._image = None

        # This attribute is set to False when creating the card from a card
        # image to ensure that the contents of the image get verified at some
        # point
        self._verified = True

        # A flag to conveniently mark whether or not this was a valid HIERARCH
        # card
        self._hierarch = False

        # If the card could not be parsed according the the FITS standard or
        # any recognized non-standard conventions, this will be True
        self._invalid = False

        self._field_specifier = None
        if not self._check_if_rvkc(keyword, value):
            # If _check_if_rvkc passes, it will handle setting the keyword and
            # value
            if keyword is not None:
                self.keyword = keyword
            if value is not None:
                self.value = value

        if comment is not None:
            self.comment = comment

        self._modified = False
        self._valuestring = None
        self._valuemodified = False

    def __repr__(self):
        return repr((self.keyword, self.value, self.comment))

    def __str__(self):
        return self.image

    def __len__(self):
        return 3

    def __getitem__(self, index):
        return (self.keyword, self.value, self.comment)[index]

    @property
    def keyword(self):
        """Returns the keyword name parsed from the card image."""
        if self._keyword is not None:
            return self._keyword
        elif self._image:
            self._keyword = self._parse_keyword()
            return self._keyword
        else:
            self.keyword = ''
            return ''

    @keyword.setter
    def keyword(self, keyword):
        """Set the key attribute; once set it cannot be modified."""
        if self._keyword is not None:
            raise AttributeError(
                'Once set, the Card keyword may not be modified')
        elif isinstance(keyword, basestring):
            # Be nice and remove trailing whitespace--some FITS code always
            # pads keywords out with spaces; leading whitespace, however,
            # should be strictly disallowed.
            keyword = keyword.rstrip()
            keyword_upper = keyword.upper()
            if len(keyword) <= 8 and self._keywd_FSC_RE.match(keyword_upper):
                # For keywords with length > 8 they will be HIERARCH cards,
                # and can have arbitrary case keywords
                if keyword_upper == 'END':
                    raise ValueError("Keyword 'END' not allowed.")
                keyword = keyword_upper
            elif self._keywd_hierarch_RE.match(keyword):
                # In prior versions of PyFITS HIERARCH cards would only be
                # created if the user-supplied keyword explicitly started with
                # 'HIERARCH '.  Now we will create them automtically for long
                # keywords, but we still want to support the old behavior too;
                # the old behavior makes it possible to create HEIRARCH cards
                # that would otherwise be recognized as RVKCs
                self._hierarch = True

                if keyword_upper[:9] == 'HIERARCH ':
                    # The user explicitly asked for a HIERARCH card, so don't
                    # bug them about it...
                    keyword = keyword[9:]
                else:
                    # We'll gladly create a HIERARCH card, but a warning is
                    # also displayed
                    warnings.warn(
                        'Keyword name %r is greater than 8 characters or '
                        'or contains spaces; a HIERARCH card will be created.' %
                        keyword, VerifyWarning)
            else:
                raise ValueError('Illegal keyword name: %r.' % keyword)
            self._keyword = keyword
            self._modified = True
        else:
            raise ValueError('Keyword name %r is not a string.' % keyword)

    @property
    @deprecated('3.1', alternative='the `.keyword` attribute')
    def key(self):
        return self.keyword

    @property
    def value(self):
        if self.field_specifier:
            return float(self._value)

        if self._value is not None:
            value = self._value
        elif self._valuestring is not None or self._image:
            self._value = self._parse_value()
            value = self._value
        else:
            self._value = value = ''

        if STRIP_HEADER_WHITESPACE() and isinstance(value, basestring):
            value = value.rstrip()

        return value

    @value.setter
    def value(self, value):
        if self._invalid:
            raise ValueError(
                'The value of invalid/unparseable cards cannot set.  Either '
                'delete this card from the header or replace it.')

        if value is None:
            value = ''
        oldvalue = self._value
        if oldvalue is None:
            oldvalue = ''

        if not isinstance(value, (basestring, int, long, float, complex, bool,
                                  Undefined, np.floating, np.integer,
                                  np.complexfloating, np.bool_)):
            raise ValueError('Illegal value: %r.' % value)

        if isinstance(value, unicode):
            try:
                # Any string value must be encodable as ASCII
                value.encode('ascii')
            except UnicodeEncodeError:
                raise ValueError(
                    'FITS header values must contain standard ASCII '
                    'characters; %r contains characters not representable in '
                    'ASCII.' % value)

        if (STRIP_HEADER_WHITESPACE() and
            (isinstance(oldvalue, basestring) and
             isinstance(value, basestring))):
            # Ignore extra whitespace when comparing the new value to the old
            different = oldvalue.rstrip() != value.rstrip()
        else:
            different = oldvalue != value

        if different:
            self._value = value
            self._modified = True
            self._valuestring = None
            self._valuemodified = True
            if self.field_specifier:
                try:
                    self._value = _int_or_float(self._value)
                except ValueError:
                    raise ValueError('value %s is not a float' %
                                     self._value)

    @value.deleter
    def value(self):
        if self._invalid:
            raise ValueError(
                'The value of invalid/unparseable cards cannot deleted.  '
                'Either delete this card from the header or replace it.')

        if not self.field_specifier:
            self.value = ''
        else:
            raise AttributeError('Values cannot be deleted from record-valued '
                                 'keyword cards')

    @property
    def rawvalue(self):
        if self._rawvalue is None and self.value:
            self._rawvalue = self.value
        return self._rawvalue

    @property
    def comment(self):
        """Get the comment attribute from the card image if not already set."""

        if self._comment is not None:
            return self._comment
        elif self._image:
            self._comment = self._parse_comment()
            return self._comment
        else:
            self.comment = ''
            return ''

    @comment.setter
    def comment(self, comment):
        if self._invalid:
            raise ValueError(
                'The comment of invalid/unparseable cards cannot set.  Either '
                'delete this card from the header or replace it.')

        if comment is None:
            comment = ''

        if isinstance(comment, unicode):
            try:
                # Any string value must be encodable as ASCII
                comment.encode('ascii')
            except UnicodeEncodeError:
                raise ValueError(
                    'FITS header comments must contain standard ASCII '
                    'characters; %r contains characters not representable in '
                    'ASCII.' % comment)

        oldcomment = self._comment
        if oldcomment is None:
            oldcomment = ''
        if comment != oldcomment:
            self._comment = comment
            self._modified = True

    @comment.deleter
    def comment(self):
        if self._invalid:
            raise ValueError(
                'The comment of invalid/unparseable cards cannot deleted.  '
                'Either delete this card from the header or replace it.')

        self.comment = ''

    @property
    def field_specifier(self):
        """
        The field-specifier of record-valued keyword cards; always `None` on
        normal cards.
        """

        # Ensure that the keyword exists and has been parsed--the will set the
        # internal _field_specifier attribute if this is a RVKC.
        if self.keyword:
            return self._field_specifier
        else:
            return None

    @field_specifier.setter
    def field_specifier(self, field_specifier):
        if not field_specifier:
            raise ValueError('The field-specifier may not be blank in '
                             'record-valued keyword cards.')
        elif not self.field_specifier:
            raise AttributeError('Cannot coerce cards to be record-valued '
                                 'keyword cards by setting the '
                                 'field_specifier attribute')
        elif field_specifier != self.field_specifier:
            self._field_specifier = field_specifier
            # The keyword need also be updated
            keyword = self._keyword.split('.', 1)[0]
            self._keyword = '.'.join([keyword, field_specifier])
            self._modified = True

    @field_specifier.deleter
    def field_specifier(self):
        raise AttributeError('The field_specifier attribute may not be '
                             'deleted from record-valued keyword cards.')

    @property
    def image(self):
        if self._image and not self._verified:
            self.verify('fix')
        if self._image is None or self._modified:
            self._image = self._format_image()
        return self._image

    @property
    @deprecated('3.1', alternative='the `.image` attribute')
    def cardimage(self):
        return self.image

    @deprecated('3.1', alternative='the `.image` attribute')
    def ascardimage(self, option='silentfix'):
        if not self._verified:
            self.verify(option)
        return self.image

    @classmethod
    def fromstring(cls, image):
        """
        Construct a `Card` object from a (raw) string. It will pad the string
        if it is not the length of a card image (80 columns).  If the card
        image is longer than 80 columns, assume it contains ``CONTINUE``
        card(s).
        """

        card = cls()
        card._image = _pad(image)
        card._verified = False
        return card

    @classmethod
    def normalize_keyword(cls, keyword):
        """
        `classmethod` to convert a keyword value that may contain a
        field-specifier to uppercase.  The effect is to raise the key to
        uppercase and leave the field specifier in its original case.

        Parameters
        ----------
        key : or str
            A keyword value or a ``keyword.field-specifier`` value

        Returns
        -------
            The converted string
        """

        match = cls._rvkc_keyword_name_RE.match(keyword)

        if match:
            return '.'.join((match.group('keyword').strip().upper(),
                             match.group('field_specifier')))
        elif len(keyword) > 9 and keyword[:9].upper() == 'HIERARCH ':
            # Remove 'HIERARCH' from HIERARCH keywords; this could lead to
            # ambiguity if there is actually a keyword card containing
            # "HIERARCH HIERARCH", but shame on you if you do that.
            return keyword[9:].strip()
        else:
            return keyword.strip().upper()

    def _check_if_rvkc(self, *args):
        """
        Determine whether or not the card is a record-valued keyword card.

        If one argument is given, that argument is treated as a full card image
        and parsed as such.  If two arguments are given, the first is treated
        as the card keyword (including the field-specifier if the card is
        intened as a RVKC), and the second as the card value OR the first value
        can be the base keyword, and the second value the 'field-specifier:
        value' string.

        If the check passes the ._keyword, ._value, and .field_specifier
        keywords are set.

        Examples
        --------

        >>> self._check_if_rvkc('DP1', 'AXIS.1: 2')
        >>> self._check_if_rvkc('DP1.AXIS.1', 2)
        >>> self._check_if_rvkc('DP1     = AXIS.1: 2')
        """

        if not ENABLE_RECORD_VALUED_KEYWORD_CARDS():
            return False

        if len(args) == 1:
            image = args[0]
            eq_idx = image.find('=')
            if eq_idx < 0 or eq_idx > 9:
                return False
            keyword, rest = image.split('=', 1)
            if keyword in self._commentary_keywords:
                return False
            match = self._rvkc_keyword_val_comm_RE.match(rest)
            if match:
                field_specifier = match.group('keyword')
                self._keyword = '.'.join((keyword.strip().upper(),
                                          field_specifier))
                self._field_specifier = field_specifier
                self._value = _int_or_float(match.group('val'))
                self._rawvalue = match.group('rawval')
                return True
        elif len(args) == 2:
            keyword, value = args
            if not isinstance(keyword, basestring):
                return False
            if keyword in self._commentary_keywords:
                return False
            match = self._rvkc_keyword_name_RE.match(keyword)
            if match and isinstance(value, (int, float)):
                field_specifier = match.group('field_specifier')
                self._keyword = '.'.join((match.group('keyword').upper(),
                                          field_specifier))
                self._field_specifier = field_specifier
                self._value = value
                return True
            if isinstance(value, basestring):
                match = self._rvkc_field_specifier_val_RE.match(value)
                if match and self._keywd_FSC_RE.match(keyword):
                    field_specifier = match.group('keyword')
                    self._keyword = '.'.join((keyword.upper(),
                                              field_specifier))
                    self._field_specifier = field_specifier
                    self._value = _int_or_float(match.group('val'))
                    self._rawvalue = value
                    return True

    def _parse_keyword(self):
        if self._check_if_rvkc(self._image):
            return self._keyword

        keyword = self._image[:8].strip()
        keyword_upper = keyword.upper()
        value_indicator = self._image.find('= ')

        special = self._commentary_keywords + ['CONTINUE']

        if keyword_upper in special or 0 <= value_indicator <= 8:
            # The value indicator should appear in byte 8, but we are flexible
            # and allow this to be fixed
            if value_indicator >= 0:
                keyword = keyword[:value_indicator]
                keyword_upper = keyword_upper[:value_indicator]

            if keyword_upper != keyword:
                self._modified = True
            return keyword_upper
        elif (keyword_upper == 'HIERARCH' and self._image[8] == ' ' and
              '=' in self._image):
            # This is valid HIERARCH card as described by the HIERARCH keyword
            # convention:
            # http://fits.gsfc.nasa.gov/registry/hierarch_keyword.html
            self._hierarch = True
            return self._image.split('=', 1)[0][9:].rstrip()
        else:
            warnings.warn('The following header keyword is invalid or follows '
                          'an unrecognized non-standard convention:\n%s' %
                          self._image)
            self._invalid = True
            return keyword

    def _parse_value(self):
        """Extract the keyword value from the card image."""

        # for commentary cards, no need to parse further
        # Likewise for invalid cards
        if self.keyword.upper() in self._commentary_keywords or self._invalid:
            return self._image[8:].rstrip()

        if self._check_if_rvkc(self._image):
            return self._value

        if len(self._image) > self.length:
            values = []
            for card in self._itersubcards():
                value = card.value.rstrip().replace("''", "'")
                if value and value[-1] == '&':
                    value = value[:-1]
                values.append(value)

            value = ''.join(values)

            self._valuestring = value
            return value

        m = self._value_NFSC_RE.match(self._split()[1])

        if m is None:
            raise VerifyError("Unparsable card (%s), fix it first with "
                              ".verify('fix')." % self.key)

        if m.group('bool') is not None:
            value = m.group('bool') == 'T'
        elif m.group('strg') is not None:
            value = re.sub("''", "'", m.group('strg'))
        elif m.group('numr') is not None:
            #  Check for numbers with leading 0s.
            numr = self._number_NFSC_RE.match(m.group('numr'))
            digt = translate(numr.group('digt'), FIX_FP_TABLE2, ' ')
            if numr.group('sign') is None:
                sign = ''
            else:
                sign = numr.group('sign')
            value = _str_to_num(sign + digt)

        elif m.group('cplx') is not None:
            #  Check for numbers with leading 0s.
            real = self._number_NFSC_RE.match(m.group('real'))
            rdigt = translate(real.group('digt'), FIX_FP_TABLE2, ' ')
            if real.group('sign') is None:
                rsign = ''
            else:
                rsign = real.group('sign')
            value = _str_to_num(rsign + rdigt)
            imag = self._number_NFSC_RE.match(m.group('imag'))
            idigt = translate(imag.group('digt'), FIX_FP_TABLE2, ' ')
            if imag.group('sign') is None:
                isign = ''
            else:
                isign = imag.group('sign')
            value += _str_to_num(isign + idigt) * 1j
        else:
            value = UNDEFINED

        if not self._valuestring:
            self._valuestring = m.group('valu')
        return value

    def _parse_comment(self):
        """Extract the keyword value from the card image."""

        # for commentary cards, no need to parse further
        # likewise for invalid/unparseable cards
        if self.keyword in Card._commentary_keywords or self._invalid:
            return ''

        if len(self._image) > self.length:
            comments = []
            for card in self._itersubcards():
                if card.comment:
                    comments.append(card.comment)
            comment = '/ ' + ' '.join(comments).rstrip()
            m = self._value_NFSC_RE.match(comment)
        else:
            m = self._value_NFSC_RE.match(self._split()[1])

        if m is not None:
            comment = m.group('comm')
            if comment:
                return comment.rstrip()
        return ''

    def _split(self):
        """
        Split the card image between the keyword and the rest of the card.
        """

        if self._image is not None:
            # If we already have a card image, don't try to rebuild a new card
            # image, which self.image would do
            image = self._image
        else:
            image = self.image

        if self.keyword in self._commentary_keywords + ['CONTINUE']:
            keyword, valuecomment = image.split(' ', 1)
        else:
            try:
                delim_index = image.index('= ')
            except ValueError:
                delim_index = None

            # The equal sign may not be any higher than column 10; anything
            # past that must be considered part of the card value
            if delim_index is None:
                keyword = image[:8]
                valuecomment = image[8:]
            elif delim_index > 10 and image[:9] != 'HIERARCH ':
                keyword = image[:8]
                valuecomment = image[8:]
            else:
                keyword, valuecomment = image.split('= ', 1)
        return keyword.strip(), valuecomment.strip()

    def _fix_keyword(self):
        if self.field_specifier:
            keyword, field_specifier = self._keyword.split('.', 1)
            self._keyword = '.'.join([keyword.upper(), field_specifier])
        else:
            self._keyword = self._keyword.upper()
        self._modified = True

    def _fix_value(self):
        """Fix the card image for fixable non-standard compliance."""

        value = None
        keyword, valuecomment = self._split()
        m = self._value_NFSC_RE.match(valuecomment)

        # for the unparsable case
        if m is None:
            try:
                value, comment = valuecomment.split('/', 1)
                self.value = value.strip()
                self.comment = comment.strip()
            except (ValueError, IndexError):
                self.value = valuecomment
            self._valuestring = self._value
            return
        elif m.group('numr') is not None:
            numr = self._number_NFSC_RE.match(m.group('numr'))
            value = translate(numr.group('digt'), FIX_FP_TABLE, ' ')
            if numr.group('sign') is not None:
                value = numr.group('sign') + value

        elif m.group('cplx') is not None:
            real = self._number_NFSC_RE.match(m.group('real'))
            rdigt = translate(real.group('digt'), FIX_FP_TABLE, ' ')
            if real.group('sign') is not None:
                rdigt = real.group('sign') + rdigt

            imag = self._number_NFSC_RE.match(m.group('imag'))
            idigt = translate(imag.group('digt'), FIX_FP_TABLE, ' ')
            if imag.group('sign') is not None:
                idigt = imag.group('sign') + idigt
            value = '(%s, %s)' % (rdigt, idigt)
        self._valuestring = value
        # The value itself has not been modified, but its serialized
        # representation (as stored in self._valuestring) has been changed, so
        # still set this card as having been modified (see ticket #137)
        self._modified = True

    def _format_keyword(self):
        if self.keyword:
            if self.field_specifier:
                return '%-8s' % self.keyword.split('.', 1)[0]
            elif self._hierarch:
                return 'HIERARCH %s ' % self.keyword
            else:
                return '%-8s' % self.keyword
        else:
            return ' ' * 8

    def _format_value(self):
        # value string
        float_types = (float, np.floating, complex, np.complexfloating)

        # Force the value to be parsed out first
        value = self.value
        # But work with the underlying raw value instead (to preserve
        # whitespace, for now...)
        value = self._value

        if self.keyword in self._commentary_keywords:
            # The value of a commentary card must be just a raw unprocessed
            # string
            value = str(value)
        elif (self._valuestring and not self._valuemodified and
              isinstance(self.value, float_types)):
            # Keep the existing formatting for float/complex numbers
            value = '%20s' % self._valuestring
        elif self.field_specifier:
            value = _format_value(self._value).strip()
            value = "'%s: %s'" % (self.field_specifier, value)
        else:
            value = _format_value(value)

        # For HIERARCH cards the value should be shortened to conserve space
        if not self.field_specifier and len(self.keyword) > 8:
            value = value.strip()

        return value

    def _format_comment(self):
        if not self.comment:
            return ''
        else:
            return ' / %s' % self._comment

    def _format_image(self):
        keyword = self._format_keyword()

        value = self._format_value()
        is_commentary = keyword.strip() in self._commentary_keywords
        if is_commentary:
            comment = ''
        else:
            comment = self._format_comment()

        # equal sign string
        delimiter = '= '
        if is_commentary:
            delimiter = ''

        # put all parts together
        output = ''.join([keyword, delimiter, value, comment])

        # For HIERARCH cards we can save a bit of space if necessary by
        # removing the space between the keyword and the equals sign; I'm
        # guessing this is part of the HIEARCH card specification
        keywordvalue_length = len(keyword) + len(delimiter) + len(value)
        if (keywordvalue_length > self.length and
                keyword.startswith('HIERARCH')):
            if (keywordvalue_length == self.length + 1 and keyword[-1] == ' '):
                output = ''.join([keyword[:-1], delimiter, value, comment])
            else:
                # I guess the HIERARCH card spec is incompatible with CONTINUE
                # cards
                raise ValueError('The keyword %s with its value is too long' %
                                 self.keyword)

        if len(output) <= self.length:
            output = '%-80s' % output
        else:
            # longstring case (CONTINUE card)
            # try not to use CONTINUE if the string value can fit in one line.
            # Instead, just truncate the comment
            if (isinstance(self.value, str) and
                len(value) > (self.length - 10)):
                output = self._format_long_image()
            else:
                warnings.warn('Card is too long, comment will be truncated.',
                              VerifyWarning)
                output = output[:Card.length]
        return output

    def _format_long_image(self):
        """
        Break up long string value/comment into ``CONTINUE`` cards.
        This is a primitive implementation: it will put the value
        string in one block and the comment string in another.  Also,
        it does not break at the blank space between words.  So it may
        not look pretty.
        """

        if self.keyword in Card._commentary_keywords:
            return self._format_long_commentary_image()

        value_length = 67
        comment_length = 64
        output = []

        # do the value string
        value_format = "'%-s&'"
        value = self._value.replace("'", "''")
        words = _words_group(value, value_length)
        for idx, word in enumerate(words):
            if idx == 0:
                headstr = '%-8s= ' % self.keyword
            else:
                headstr = 'CONTINUE  '
            value = value_format % word
            output.append('%-80s' % (headstr + value))

        # do the comment string
        comment_format = "%-s"
        if self.comment:
            words = _words_group(self.comment, comment_length)
            for word in words:
                comment = "CONTINUE  '&' / " + comment_format % word
                output.append('%-80s' % comment)

        return ''.join(output)

    def _format_long_commentary_image(self):
        """
        If a commentary card's value is too long to fit on a single card, this
        will render the card as multiple consecutive commentary card of the
        same type.
        """

        maxlen = Card.length - len(self.keyword) - 1
        value = self._format_value()
        output = []
        idx = 0
        while idx < len(value):
            output.append(str(Card(self.keyword, value[idx:idx + maxlen])))
            idx += maxlen
        return ''.join(output)

    def _verify(self, option='warn'):
        self._verified = True

        errs = _ErrList([])
        fix_text = 'Fixed %r card to meet the FITS standard.' % self.keyword

        # Don't try to verify cards that already don't meet any recognizable
        # standard
        if self._invalid:
            return errs

        # verify the equal sign position
        if (self.keyword not in self._commentary_keywords and
            (self._image and self._image[:8].upper() != 'HIERARCH' and
             self._image.find('=') != 8)):
            errs.append(self.run_option(
                option,
                err_text='Card %r is not FITS standard (equal sign not '
                         'at column 8).' % self.keyword,
                fix_text=fix_text,
                fix=self._fix_value))

        # verify the key, it is never fixable
        # always fix silently the case where "=" is before column 9,
        # since there is no way to communicate back to the _keys.
        if self._image and self._image[:8].upper() == 'HIERARCH':
            pass
        else:
            keyword = self.keyword
            if self.field_specifier:
                keyword = keyword.split('.', 1)[0]

            if keyword != keyword.upper():
            # Keyword should be uppercase unless it's a HIERARCH card
                errs.append(self.run_option(
                    option,
                    err_text='Card keyword %r is not upper case.' % keyword,
                    fix_text=fix_text,
                    fix=self._fix_keyword))
            elif not self._keywd_FSC_RE.match(keyword):
                errs.append(self.run_option(
                    option,
                    err_text='Illegal keyword name %s' % repr(keyword),
                    fixable=False))

        # verify the value, it may be fixable
        keyword, valuecomment = self._split()
        m = self._value_FSC_RE.match(valuecomment)
        if not (m or self.keyword in self._commentary_keywords):
            errs.append(self.run_option(
                option,
                err_text='Card %r is not FITS standard (invalid value '
                         'string: %s).' % (self.keyword, valuecomment),
                fix_text=fix_text,
                fix=self._fix_value))

        # verify the comment (string), it is never fixable
        m = self._value_NFSC_RE.match(valuecomment)
        if m is not None:
            comment = m.group('comm')
            if comment is not None:
                if not self._comment_FSC_RE.match(comment):
                    errs.append(self.run_option(
                        option,
                        err_text='Unprintable string %r' % comment,
                        fixable=False))

        return errs

    def _itersubcards(self):
        """
        If the card image is greater than 80 characters, it should consist of a
        normal card followed by one or more CONTINUE card.  This method returns
        the subcards that make up this logical card.
        """

        ncards = len(self._image) // Card.length

        for idx in xrange(0, Card.length * ncards, Card.length):
            card = Card.fromstring(self._image[idx:idx + Card.length])
            if idx > 0 and card.keyword.upper() != 'CONTINUE':
                raise VerifyError(
                        'Long card images must have CONTINUE cards after '
                        'the first card.')

            if not isinstance(card.value, str):
                raise VerifyError('CONTINUE cards must have string values.')

            yield card


def create_card(key='', value='', comment=''):
    return Card(key, value, comment)
create_card.__doc__ = Card.__init__.__doc__
# For API backwards-compatibility
createCard = deprecated('3.0', name='createCard',
                        alternative=':meth:`Card.__init__`')(create_card)
create_card = deprecated('3.1', name='create_card',
                         alternative=':meth:`Card.__init__`',
                         pending=True)(create_card)


def create_card_from_string(input):
    return Card.fromstring(input)
create_card_from_string.__doc__ = Card.fromstring.__doc__
# For API backwards-compat
createCardFromString = \
    deprecated('3.0', name='createCardFromString',
               alternative=':meth:`Card.fromstring`')(create_card_from_string)
create_card_from_string = deprecated('3.1', name='create_card_from_string',
                                     alternative=':meth:`Card.fromstring`',
                                     pending=True)(create_card_from_string)


def upper_key(key):
    return Card.normalize_keyword(key)
upper_key.__doc__ = Card.normalize_keyword.__doc__
# For API backwards-compat
upperKey = deprecated('3.0', name='upperKey',
                      alternative=':meth:`Card.normalize_keyword`')(upper_key)
upper_key = deprecated('3.1', name='upper_key',
                       alternative=':meth:`Card.normalize_keyword`',
                       pending=True)(upper_key)


def _int_or_float(s):
    """
    Converts an a string to an int if possible, or to a float.

    If the string is neither a string or a float a value error is raised.
    """

    if isinstance(s, float):
        # Already a float so just pass through
        return s

    try:
        return int(s)
    except (ValueError, TypeError):
        try:
            return float(s)
        except (ValueError, TypeError), e:
            raise ValueError(str(e))


def _format_value(value):
    """
    Converts a card value to its appropriate string representation as
    defined by the FITS format.
    """

    # string value should occupies at least 8 columns, unless it is
    # a null string
    if isinstance(value, basestring):
        if value == '':
            return "''"
        else:
            exp_val_str = value.replace("'", "''")
            val_str = "'%-8s'" % exp_val_str
            return '%-20s' % val_str

    # must be before int checking since bool is also int
    elif isinstance(value, (bool, np.bool_)):
        return '%20s' % repr(value)[0]  # T or F

    elif _is_int(value):
        return '%20d' % value

    elif isinstance(value, (float, np.floating)):
        return '%20s' % _format_float(value)

    elif isinstance(value, (complex, np.complexfloating)):
        val_str = '(%s, %s)' % (_format_float(value.real),
                                _format_float(value.imag))
        return '%20s' % val_str

    elif isinstance(value, Undefined):
        return ''
    else:
        return ''


def _format_float(value):
    """Format a floating number to make sure it gets the decimal point."""

    value_str = '%.16G' % value
    if '.' not in value_str and 'E' not in value_str:
        value_str += '.0'
    elif 'E' in value_str:
        # On some Windows builds of Python (and possibly other platforms?) the
        # exponent is zero-padded out to, it seems, three digits.  Normalize
        # the format to pad only to two digits.
        significand, exponent = value_str.split('E')
        if exponent[0] in ('+', '-'):
            sign = exponent[0]
            exponent = exponent[1:]
        else:
            sign = ''
        value_str = '%sE%s%02d' % (significand, sign, int(exponent))

    # Limit the value string to at most 20 characters.
    str_len = len(value_str)

    if str_len > 20:
        idx = value_str.find('E')

        if idx < 0:
            value_str = value_str[:20]
        else:
            value_str = value_str[:20 - (str_len - idx)] + value_str[idx:]

    return value_str


def _pad(input):
    """Pad blank space to the input string to be multiple of 80."""

    _len = len(input)
    if _len == Card.length:
        return input
    elif _len > Card.length:
        strlen = _len % Card.length
        if strlen == 0:
            return input
        else:
            return input + ' ' * (Card.length - strlen)

    # minimum length is 80
    else:
        strlen = _len % Card.length
        return input + ' ' * (Card.length - strlen)
