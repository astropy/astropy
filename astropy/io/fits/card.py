import re
import string
import sys
import warnings

import numpy as np

from pyfits.util import _str_to_num, _is_int, deprecated, maketrans, translate
from pyfits.verify import _Verify, _ErrList


__all__ = ['Card', 'CardList', 'RecordValuedKeywordCard', 'create_card',
           'create_card_from_string', 'upper_key', 'createCard',
           'createCardFromString', 'upperKey', 'Undefined']


FIX_FP_TABLE = maketrans('de', 'DE')
FIX_FP_TABLE2 = maketrans('dD', 'eE')


class Undefined:
    """Undefined value."""

    def __init__(self):
        # This __init__ is required to be here for Sphinx documentation
        pass
UNDEFINED = Undefined()


class Card(_Verify):

    # string length of a card
    length = 80

    # String for a FITS standard compliant (FSC) keyword.
    _keywd_FSC = r'[A-Z0-9_-]* *$'
    _keywd_FSC_RE = re.compile(_keywd_FSC)

    # A number sub-string, either an integer or a float in fixed or
    # scientific notation.  One for FSC and one for non-FSC (NFSC) format:
    # NFSC allows lower case of DE for exponent, allows space between sign,
    # digits, exponent sign, and exponents
    _digits_FSC = r'(\.\d+|\d+(\.\d*)?)([DE][+-]?\d+)?'
    _digits_NFSC = r'(\.\d+|\d+(\.\d*)?) *([deDE] *[+-]? *\d+)?'
    _numr_FSC = r'[+-]?' + _digits_FSC
    _numr_NFSC = r'[+-]? *' + _digits_NFSC

    # This regex helps delete leading zeros from numbers, otherwise
    # Python might evaluate them as octal values.
    _number_FSC_RE = re.compile(r'(?P<sign>[+-])?0*(?P<digt>%s)'
                                % _digits_FSC)
    _number_NFSC_RE = re.compile(r'(?P<sign>[+-])? *0*(?P<digt>%s)'
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
                r'\'(?P<strg>([ -~]+?|\'\'|)) *?\'(?=$|/| )|'
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

    # keys of commentary cards
    _commentary_keys = ['', 'COMMENT', 'HISTORY']

    def __new__(cls, key='', value='', comment=''):
        """
        Return the appropriate Card subclass depending on they key and/or
        value.
        """

        klass = cls
        val_str = _value_to_string(value)
        if len(val_str) > (Card.length - 10):
            klass = _ContinueCard
        elif isinstance(key, basestring) and key.upper().strip() == 'HIERARCH':
            klass = _HierarchCard
        return super(Card, cls).__new__(klass)

    def __init__(self, key='', value='', comment=''):
        """
        Construct a card from `key`, `value`, and (optionally) `comment`.
        Any specifed arguments, except defaults, must be compliant to FITS
        standard.

        Parameters
        ----------
        key : str, optional
            keyword name

        value : str, optional
            keyword value

        comment : str, optional
            comment
        """

        if key != '' or value != '' or comment != '':
            self.key = key
            self._update_value(value)
            self._update_comment(comment)
            # for commentary cards, value can only be strings and there
            # is no comment
            if self.key in Card._commentary_keys:
                if not isinstance(self.value, basestring):
                    raise ValueError('Value in a commentary card must be a '
                                     'string.')
        else:
            self._cardimage = ' ' * 80
        self._modified = False

    def __str__(self):
        return self.cardimage

    @property
    def key(self):
        """Returns the keyword name parsed from the card image."""

        if not hasattr(self, '_key'):
            head = self._get_key_string()
            self._key = head.strip().upper()
        return self._key

    @key.setter
    def key(self, val):
        """Set the key attribute; once set it cannot be modified."""

        if hasattr(self, '_key'):
            raise AttributeError('Keyword name cannot be modified.')

        if isinstance(val, basestring):
            val = val.strip()
            if len(val) <= 8:
                val = val.upper()
                if val == 'END':
                    raise ValueError("Keyword 'END' not allowed.")
                self._check_key(val)
            else:
                if val[:8].upper() == 'HIERARCH':
                    val = val[8:].strip()
                    self.__class__ = _HierarchCard
                else:
                    raise ValueError('Keyword name %s is too long (> 8), use HIERARCH.'
                                     % val)
        else:
            raise ValueError('Keyword name %s is not a string.' % repr(val))
        self._key = val
        self._modified = True

    @property
    def value(self):
        """Get the value attribute from the card image if not already set."""

        if not hasattr(self, '_value'):
            # self._value assigned directly to avoid the checks in
            # _setvalue() (and to avoid setting _value_modified = True)
            self._value = self._extract_value()
        return self._value

    @value.setter
    def value(self, val):
        self._update_value(val)
        if hasattr(self, '_cardimage'):
            # Make sure these are saved from the old cardimage before deleting
            # it
            self.key
            self.comment
            del self._cardimage

    def _update_value(self, val):
        """Set the value attribute."""

        if isinstance(val, (str, int, long, float, complex, bool, Undefined,
                            np.floating, np.integer, np.complexfloating)):
            if isinstance(val, str):
                self._check_text(val)
            if not hasattr(self, '_value') or self._value != val:
                self._modified = True
            self._value_modified = True
            self._value = val
        else:
            raise ValueError('Illegal value %s.' % repr(val))

    @property
    def comment(self):
        """Get the comment attribute from the card image if not already set."""

        if not hasattr(self, '_comment'):
            self._comment = self._extract_comment()
        return self._comment

    @comment.setter
    def comment(self, val):
        self._update_comment(val)
        if hasattr(self, '_cardimage'):
            # Make sure these are saved from the old cardimage before deleting
            # it
            self.key
            self.value
            del self._cardimage

    def _update_comment(self, val):
        """Set the comment attribute."""

        if isinstance(val, basestring):
            self._check_text(val)
        else:
            if val is not None:
                raise ValueError('Comment %s is not a string.' % repr(val))
        if not hasattr(self, '_comment') or self._comment != val:
            self._modified = True
        self._comment = val

    @property
    def cardimage(self):
        if not hasattr(self, '_cardimage'):
            self._ascardimage()
        return self._cardimage

    def _update_cardimage(self):
        """This used to do what _generate_cardimage() does, plus setting te
        value of self._cardimage.  Now it's kept mostly for symmetry.
        """

        self._cardimage = self._generate_cardimage()

    def _generate_cardimage(self):
        """Generate a (new) card image from the attributes: `key`, `value`,
        and `comment`.  Core code for `ascardimage`.
        """

        key_str = self._format_key()
        val_str = self._format_value()
        is_commentary = key_str.strip() in Card._commentary_keys
        if is_commentary:
            comment_str = ''
        else:
            comment_str = self._format_comment()

        # equal sign string
        eq_str = '= '
        if is_commentary:  # not using self.key
            eq_str = ''
            if hasattr(self, '_value'):
                # Value must be a string in commentary cards
                val_str = str(self.value)

        # put all parts together
        output = ''.join([key_str, eq_str, val_str, comment_str])

        # need this in case card-with-continue's value is shortened
        if not isinstance(self, _HierarchCard) and \
           not isinstance(self, RecordValuedKeywordCard):
            self.__class__ = Card
        else:
            key_val_len = len(key_str) + len(eq_str) + len(val_str)
            if key_val_len > Card.length:
                if isinstance(self, _HierarchCard) and \
                   key_val_len == Card.length + 1 and \
                   key_str[-1] == ' ':
                    output = ''.join([key_str[:-1], eq_str, val_str,
                                      comment_str])
                else:
                    raise ValueError('The keyword %s with its value is too '
                                     'long.' % self.key)

        if len(output) <= Card.length:
            output = '%-80s' % output

        # longstring case (CONTINUE card)
        else:
            # try not to use CONTINUE if the string value can fit in one line.
            # Instead, just truncate the comment
            if isinstance(self.value, str) and \
               len(val_str) > (Card.length - 10):
                self.__class__ = _ContinueCard
                output = self._breakup_strings()
            else:
                warnings.warn('Card is too long, comment is truncated.')
                output = output[:Card.length]

        return output

    def _format_key(self):
        # keyword string (check both _key and _cardimage attributes to avoid an
        # infinite loop
        if hasattr(self, '_key') or hasattr(self, '_cardimage'):
            return '%-8s' % self.key
        else:
            return ' ' * 8

    def _format_value(self):
        # value string
        if not (hasattr(self, '_value') or hasattr(self, '_cardimage')):
            return ''
        elif (not hasattr(self, '_value_modified') or
              not self._value_modified) and \
             isinstance(self.value, (float, np.floating, complex,
                                     np.complexfloating)):
            # Keep the existing formatting for float/complex numbers
            return '%20s' % self._valuestring
        else:
            return _value_to_string(self.value)

    def _format_comment(self):
        # comment string
        if hasattr(self, '_comment') or hasattr(self, '_cardimage'):
            if not self.comment:
                return ''
            else:
                return ' / ' + self._comment
        else:
            return ''

    # TODO: Wouldn't 'verification' be a better name for the 'option' keyword
    # argument?  'option' is pretty vague.
    @deprecated(alternative='the .cardimage attribute')
    def ascardimage(self, option='silentfix'):
        return self._ascardimage(option)

    def _ascardimage(self, option='silentfix'):
        """
        Generate a (new) card image from the attributes: `key`, `value`,
        and `comment`, or from raw string.

        Parameters
        ----------
        option : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  See :ref:`verify` for more info.
        """

        # Only if the card image already exist (to avoid infinite loop),
        # fix it first.
        if hasattr(self, '_cardimage'):
            self._check(option)
        self._update_cardimage()
        return self._cardimage

    @classmethod
    def fromstring(cls, cardimage):
        """
        Construct a `Card` object from a (raw) string. It will pad the
        string if it is not the length of a card image (80 columns).
        If the card image is longer than 80 columns, assume it
        contains ``CONTINUE`` card(s).
        """

        cardimage = _pad(cardimage)

        if cardimage[:8].upper() == 'HIERARCH':
            card = _HierarchCard()
        # for card image longer than 80, assume it contains CONTINUE card(s).
        elif len(cardimage) > Card.length:
            card = _ContinueCard()
        else:
            card = cls()

        card._cardimage = cardimage
        return card

    def _check_text(self, val):
        """Verify `val` to be printable ASCII text."""

        if Card._comment_FSC_RE.match(val) is None:
            self._err_text = 'Unprintable string %s' % repr(val)
            self._fixable = False
            raise ValueError(self._err_text)

    def _check_key(self, key):
        """Verify the keyword `key` to be FITS standard."""

        # use repr (not str) in case of control character
        if not Card._keywd_FSC_RE.match(key):
            self._err_text = 'Illegal keyword name %s' % repr(key)
            self._fixable = False
            raise ValueError(self._err_text)

    def _extract_value(self):
        """Extract the keyword value from the card image."""

        # for commentary cards, no need to parse further
        if self.key in Card._commentary_keys:
            return self.cardimage[8:].rstrip()

        valu = self._check(option='parse')

        if valu is None:
            raise ValueError("Unparsable card (%s), fix it first with "
                             ".verify('fix')." % self.key)
        if valu.group('bool') is not None:
            _val = valu.group('bool') == 'T'
        elif valu.group('strg') is not None:
            _val = re.sub("''", "'", valu.group('strg'))
        elif valu.group('numr') is not None:
            #  Check for numbers with leading 0s.
            numr = Card._number_NFSC_RE.match(valu.group('numr'))
            _digt = translate(numr.group('digt'), FIX_FP_TABLE2, ' ')
            if numr.group('sign') is None:
                _sign = ''
            else:
                _sign = numr.group('sign')
            _val = _str_to_num(_sign + _digt)

        elif valu.group('cplx') is not None:
            #  Check for numbers with leading 0s.
            real = Card._number_NFSC_RE.match(valu.group('real'))
            _rdigt = translate(real.group('digt'), FIX_FP_TABLE2, ' ')
            if real.group('sign') is None:
                _rsign = ''
            else:
                _rsign = real.group('sign')
            _val = _str_to_num(_rsign + _rdigt)
            imag  = Card._number_NFSC_RE.match(valu.group('imag'))
            _idigt = translate(imag.group('digt'), FIX_FP_TABLE2, ' ')
            if imag.group('sign') is None:
                _isign = ''
            else:
                _isign = imag.group('sign')
            _val += _str_to_num(_isign + _idigt)*1j
        else:
            _val = UNDEFINED

        if not hasattr(self, '_valuestring'):
            self._valuestring = valu.group('valu')
        if not hasattr(self, '_value_modified'):
            self._value_modified = False
        return _val

    def _extract_comment(self):
        """Extract the keyword value from the card image."""

        # for commentary cards, no need to parse further
        if self.key in Card._commentary_keys:
            return ''

        valu = self._check(option='parse')
        if valu is not None:
            comm = valu.group('comm')
            if isinstance(comm, str):
                return  comm.rstrip()
        return ''

    def _fix_value(self, input):
        """Fix the card image for fixable non-standard compliance."""

        _val_str = None

        # for the unparsable case
        if input is None:
            _tmp = self._get_value_comment_string()
            try:
                slash_loc = _tmp.index("/")
                self._value = _tmp[:slash_loc].strip()
                self._comment = _tmp[slash_loc + 1:].strip()
            except:
                self._value = _tmp.strip()

        elif input.group('numr') is not None:
            numr = Card._number_NFSC_RE.match(input.group('numr'))
            _val_str = translate(numr.group('digt'), FIX_FP_TABLE, ' ')
            if numr.group('sign') is not None:
                _val_str = numr.group('sign')+_val_str

        elif input.group('cplx') is not None:
            real  = Card._number_NFSC_RE.match(input.group('real'))
            _realStr = translate(real.group('digt'), FIX_FP_TABLE, ' ')
            if real.group('sign') is not None:
                _realStr = real.group('sign')+_realStr

            imag  = Card._number_NFSC_RE.match(input.group('imag'))
            _imagStr = translate(imag.group('digt'), FIX_FP_TABLE, ' ')
            if imag.group('sign') is not None:
                _imagStr = imag.group('sign') + _imagStr
            _val_str = '(%s, %s)' % (_realStr, _imagStr)

        self._valuestring = _val_str
        self._update_cardimage()

    def _locate_eq(self):
        """
        Locate the equal sign in the card image before column 10 and
        return its location.  It returns `None` if equal sign is not
        present, or it is a commentary card.
        """

        # no equal sign for commentary cards (i.e. part of the string value)
        _key = self.cardimage[:8].strip().upper()
        if _key in Card._commentary_keys:
            eq_loc = None
        else:
            if _key == 'HIERARCH':
                _limit = Card.length
            else:
                _limit = 10
            try:
                eq_loc = self.cardimage[:_limit].index("=")
            except:
                eq_loc = None
        return eq_loc

    def _get_key_string(self):
        """
        Locate the equal sign in the card image and return the string
        before the equal sign.  If there is no equal sign, return the
        string before column 9.
        """

        eq_loc = self._locate_eq()
        if eq_loc is None:
            eq_loc = 8
        start = 0
        if self.cardimage[:8].upper() == 'HIERARCH':
            start = 8
            self.__class__ = _HierarchCard
        return self.cardimage[start:eq_loc]

    def _get_value_comment_string(self):
        """
        Locate the equal sign in the card image and return the string
        after the equal sign.  If there is no equal sign, return the
        string after column 8.
        """

        eq_loc = self._locate_eq()
        if eq_loc is None:
            eq_loc = 7
        return self.cardimage[eq_loc+1:]

    def _check(self, option='ignore'):
        """Verify the card image with the specified option."""

        self._err_text = ''
        self._fix_text = ''
        self._fixable = True

        if option == 'ignore':
            return
        elif option == 'parse':

            # check the value only, no need to check key and comment for 'parse'
            result = Card._value_NFSC_RE.match(self._get_value_comment_string())

            # if not parsable (i.e. everything else) result = None
            return result
        else:

            # verify the equal sign position
            if self.key not in Card._commentary_keys and \
               self.cardimage.find('=') != 8:
                if option in ['exception', 'warn']:
                    self._err_text = \
                        'Card image is not FITS standard (equal sign not at ' \
                        'column 8).'
                    raise ValueError(self._err_text + '\n%s' % self.cardimage)
                elif option in ['fix', 'silentfix']:
                    result = self._check('parse')
                    self._fix_value(result)
                    if option == 'fix':
                        self._fix_text = \
                            'Fixed card to be FITS standard.: %s' % self.key

            # verify the key, it is never fixable
            # always fix silently the case where "=" is before column 9,
            # since there is no way to communicate back to the _keys.
            self._check_key(self.key)

            # verify the value, it may be fixable
            value = self._get_value_comment_string()
            result = Card._value_FSC_RE.match(value)
            if result is not None or self.key in Card._commentary_keys:
                return result
            else:
                if option in ['fix', 'silentfix']:
                    result = self._check('parse')
                    self._fix_value(result)
                    if option == 'fix':
                        self._fix_text = \
                            'Fixed card to be FITS standard.: %s' % self.key
                else:
                    self._err_text = \
                        'Card image is not FITS standard (unparsable value ' \
                        'string: %s).' % value
                    raise ValueError(self._err_text + '\n%s' % self.cardimage)

            # verify the comment (string), it is never fixable
            if result is not None:
                comm = result.group('comm')
                if comm is not None:
                    self._check_text(comm)

    def _ncards(self):
        return len(self.cardimage) // Card.length

    def _verify(self, option='warn'):
        """Card class verification method."""

        err = _ErrList([])
        try:
            self._check(option)
        except ValueError:
            # Trapping the ValueError raised by _check method.  Want execution to continue while printing
            # exception message.
            pass
        err.append(self.run_option(option, err_text=self._err_text,
                   fix_text=self._fix_text, fixable=self._fixable))

        return err


class RecordValuedKeywordCard(Card):
    """
    Class to manage record-valued keyword cards as described in the
    FITS WCS Paper IV proposal for representing a more general
    distortion model.

    Record-valued keyword cards are string-valued cards where the
    string is interpreted as a definition giving a record field name,
    and its floating point value.  In a FITS header they have the
    following syntax::

        keyword = 'field-specifier: float'

    where `keyword` is a standard eight-character FITS keyword name,
    `float` is the standard FITS ASCII representation of a floating
    point number, and these are separated by a colon followed by a
    single blank.  The grammar for field-specifier is::

        field-specifier:
            field
            field-specifier.field

        field:
            identifier
            identifier.index

    where `identifier` is a sequence of letters (upper or lower case),
    underscores, and digits of which the first character must not be a
    digit, and `index` is a sequence of digits.  No blank characters
    may occur in the field-specifier.  The `index` is provided
    primarily for defining array elements though it need not be used
    for that purpose.

    Multiple record-valued keywords of the same name but differing
    values may be present in a FITS header.  The field-specifier may
    be viewed as part of the keyword name.

    Some examples follow::

        DP1     = 'NAXIS: 2'
        DP1     = 'AXIS.1: 1'
        DP1     = 'AXIS.2: 2'
        DP1     = 'NAUX: 2'
        DP1     = 'AUX.1.COEFF.0: 0'
        DP1     = 'AUX.1.POWER.0: 1'
        DP1     = 'AUX.1.COEFF.1: 0.00048828125'
        DP1     = 'AUX.1.POWER.1: 1'
    """

    #
    # A group of class level regular expression definitions that allow the
    # extraction of the key, field-specifier, value, and comment from a
    # card string.
    #
    identifier = r'[a-zA-Z_]\w*'
    field = identifier + r'(\.\d+)?'
    field_specifier_s = r'%s(\.%s)*' % (field, field)
    field_specifier_val = r'(?P<keyword>%s): (?P<val>%s\s*)' \
                          % (field_specifier_s, Card._numr_FSC)
    field_specifier_NFSC_val = r'(?P<keyword>%s): (?P<val>%s\s*)' \
                               % (field_specifier_s, Card._numr_NFSC)
    keyword_val = r'\'%s\'' % field_specifier_val
    keyword_NFSC_val = r'\'%s\'' % field_specifier_NFSC_val
    keyword_val_comm = r' +%s *(/ *(?P<comm>[ -~]*))?$' % keyword_val
    keyword_NFSC_val_comm = r' +%s *(/ *(?P<comm>[ -~]*))?$' % keyword_NFSC_val
    #
    # regular expression to extract the field specifier and value from
    # a card image (ex. 'AXIS.1: 2'), the value may not be FITS Standard
    # Compliant
    #
    field_specifier_NFSC_image_RE = re.compile(field_specifier_NFSC_val)
    #
    # regular expression to extract the field specifier and value from
    # a card value; the value may not be FITS Standard Compliant
    # (ex. 'AXIS.1: 2.0e5')
    #
    field_specifier_NFSC_val_RE = re.compile(field_specifier_NFSC_val + r'$')
    #
    # regular expression to extract the key and the field specifier from a
    # string that is being used to index into a card list that contains
    # record value keyword cards (ex. 'DP1.AXIS.1')
    #
    keyword_name_RE = re.compile(r'(?P<key>%s)\.(?P<field_spec>%s)$'
                                 % (identifier, field_specifier_s))
    #
    # regular expression to extract the field specifier and value and comment
    # from the string value of a record value keyword card
    # (ex "'AXIS.1: 1' / a comment")
    #
    keyword_val_comm_RE = re.compile(keyword_val_comm)
    #
    # regular expression to extract the field specifier and value and comment
    # from the string value of a record value keyword card  that is not FITS
    # Standard Complient (ex "'AXIS.1: 1.0d12' / a comment")
    #
    keyword_NFSC_val_comm_RE = re.compile(keyword_NFSC_val_comm)

    def __init__(self, key='', value='', comment=''):
        """
        Parameters
        ----------
        key : str, optional
            The key, either the simple key or one that contains
            a field-specifier

        value : str, optional
            The value, either a simple value or one that contains a
            field-specifier

        comment : str, optional
            The comment
        """

        mo = self.keyword_name_RE.match(key)

        if mo:
            self._field_specifier = mo.group('field_spec')
            key = mo.group('key')
        else:
            if isinstance(value, str):
                if value:
                    mo = self.field_specifier_NFSC_val_RE.match(value)

                    if mo:
                        self._field_specifier = mo.group('keyword')
                        value = mo.group('val')
                        # The value should be a float, though we don't coerce
                        # ints into floats.  Anything else should be a value
                        # error
                        try:
                            value = int(value)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                raise ValueError(
                                    "Record-valued keyword card value must be "
                                    "a floating point or integer value.")
                    else:
                        raise ValueError(
                            "Value %s must be in the form "
                            "field_specifier: value (ex. 'NAXIS: 2')" % value)
            else:
                raise ValueError('value %s is not a string' % value)

        super(RecordValuedKeywordCard, self).__init__(key, value, comment)

    @property
    def field_specifier(self):
       self._extract_value()
       return self._field_specifier

    @property
    def raw(self):
        """
        Return this card as a normal Card object not parsed as a record-valued
        keyword card.  Note that this returns a copy, so that modifications to
        it do not update the original record-valued keyword card.
        """

        key = super(RecordValuedKeywordCard, self).key
        return Card(key, self.strvalue(), self.comment)

    @property
    def key(self):
        key = super(RecordValuedKeywordCard, self).key
        if not hasattr(self, '_field_specifier'):
            return key
        return '%s.%s' % (key, self._field_specifier)

    @key.setter
    def key(self, value):
        Card.key.fset(self, value)

    @property
    def value(self):
        """The RVKC value should always be returned as a float."""

        return float(super(RecordValuedKeywordCard, self).value)

    @value.setter
    def value(self, val):
        if not isinstance(val, float):
            try:
                val = int(val)
            except ValueError:
                try:
                    val = float(val)
                except:
                    raise ValueError('value %s is not a float' % val)
        Card.value.fset(self, val)

    #
    # class method definitins
    #

    @classmethod
    def coerce(cls, card):
        """
        Coerces an input `Card` object to a `RecordValuedKeywordCard`
        object if the value of the card meets the requirements of this
        type of card.

        Parameters
        ----------
        card : `Card` object
            A `Card` object to coerce

        Returns
        -------
        card
            - If the input card is coercible:

                a new `RecordValuedKeywordCard` constructed from the
                `key`, `value`, and `comment` of the input card.

            - If the input card is not coercible:

                the input card
        """

        mo = cls.field_specifier_NFSC_val_RE.match(card.value)
        if mo:
            return cls(card.key, card.value, card.comment)
        else:
            return card

    @classmethod
    def upper_key(cls, key):
        """
        `classmethod` to convert a keyword value that may contain a
        field-specifier to uppercase.  The effect is to raise the
        key to uppercase and leave the field specifier in its original
        case.

        Parameters
        ----------
        key : int or str
            A keyword value that could be an integer, a key, or a
            `key.field-specifier` value

        Returns
        -------
        Integer input
            the original integer key

        String input
            the converted string
        """

        if _is_int(key):
            return key

        mo = cls.keyword_name_RE.match(key)

        if mo:
            return mo.group('key').strip().upper() + '.' + \
                   mo.group('field_spec')
        else:
            return key.strip().upper()
    # For API backwards-compatibility
    upperKey = \
        deprecated(name='upperKey', alternative='upper_key()')(upper_key)

    @classmethod
    def valid_key_value(cls, key, value=0):
        """
        Determine if the input key and value can be used to form a
        valid `RecordValuedKeywordCard` object.  The `key` parameter
        may contain the key only or both the key and field-specifier.
        The `value` may be the value only or the field-specifier and
        the value together.  The `value` parameter is optional, in
        which case the `key` parameter must contain both the key and
        the field specifier.

        Parameters
        ----------
        key : str
            The key to parse

        value : str or float-like, optional
            The value to parse

        Returns
        -------
        valid input : A list containing the key, field-specifier, value

        invalid input : An empty list

        Examples
        --------

        >>> valid_key_value('DP1','AXIS.1: 2')
        >>> valid_key_value('DP1.AXIS.1', 2)
        >>> valid_key_value('DP1.AXIS.1')
        """

        rtnKey = rtnFieldSpec = rtnValue = ''
        myKey = cls.upper_key(key)

        if isinstance(myKey, basestring):
            validKey = cls.keyword_name_RE.match(myKey)

            if validKey:
               try:
                   rtnValue = float(value)
               except ValueError:
                   pass
               else:
                   rtnKey = validKey.group('key')
                   rtnFieldSpec = validKey.group('field_spec')
            else:
                if isinstance(value, str) and \
                Card._keywd_FSC_RE.match(myKey) and len(myKey) < 9:
                    validValue = cls.field_specifier_NFSC_val_RE.match(value)
                    if validValue:
                        rtnFieldSpec = validValue.group('keyword')
                        rtnValue = validValue.group('val')
                        rtnKey = myKey

        if rtnFieldSpec:
            return [rtnKey, rtnFieldSpec, rtnValue]
        else:
            return []
    # For API backwards-compatibility
    validKeyValue = \
        deprecated(name='validKeyValue',
                   alternative='valid_key_value()')(valid_key_value)

    @classmethod
    def create(cls, key='', value='', comment=''):
        """
        Create a card given the input `key`, `value`, and `comment`.
        If the input key and value qualify for a
        `RecordValuedKeywordCard` then that is the object created.
        Otherwise, a standard `Card` object is created.

        Parameters
        ----------
        key : str, optional
            The key

        value : str, optional
            The value

        comment : str, optional
            The comment

        Returns
        -------
        card
            Either a `RecordValuedKeywordCard` or a `Card` object.
        """

        if not cls.valid_key_value(key, value):
            # This should be just a normal card
            cls = Card

        return cls(key, value, comment)
    # For API backwards-compatibility
    createCard = deprecated(name='createCard', alternative='create()')(create)

    @classmethod
    def fromstring(cls, input):
        """
        Create a card given the `input` string.  If the `input` string
        can be parsed into a key and value that qualify for a
        `RecordValuedKeywordCard` then that is the object created.
        Otherwise, a standard `Card` object is created.

        Parameters
        ----------
        input : str
            The string representing the card

        Returns
        -------
        card
            either a `RecordValuedKeywordCard` or a `Card` object
        """

        idx1 = input.find("'") + 1
        idx2 = input.rfind("'")

        if idx2 > idx1 and idx1 >= 0 and \
           cls.valid_key_value('', value=input[idx1:idx2]):
            # This calls Card.fromstring, but with the RecordValuedKeywordClass
            # as the cls argument (causing an RVKC to be created)
            return super(RecordValuedKeywordCard, cls).fromstring(input)
        else:
            # This calls Card.fromstring directly, creating a plain Card
            # object.
            return Card.fromstring(input)

    # For API backwards-compatibility
    createCardFromString = deprecated(name='createCardFromString',
                                      alternative='fromstring()')(fromstring)

    def _update_cardimage(self):
        """
        Generate a (new) card image from the attributes: `key`, `value`,
        `field_specifier`, and `comment`.  Core code for `ascardimage`.
        """

        super(RecordValuedKeywordCard, self)._update_cardimage()
        eqloc = self._cardimage.index('=')
        slashloc = self._cardimage.find('/')

        if hasattr(self, '_value_modified') and self._value_modified:
            # Bypass the automatic coertion to float here, so that values like
            # '2' will still be rendered as '2' instead of '2.0'
            value = super(RecordValuedKeywordCard, self).value
            val_str = _value_to_string(value).strip()
        else:
            val_str = self._valuestring

        val_str = "'%s: %s'" % (self._field_specifier, val_str)
        val_str = '%-20s' % val_str

        output = self._cardimage[:eqloc+2] + val_str

        if slashloc > 0:
            output = output + self._cardimage[slashloc-1:]

        if len(output) <= Card.length:
            output = '%-80s' % output

        self._cardimage = output


    def _extract_value(self):
        """Extract the keyword value from the card image."""

        valu = self._check(option='parse')

        if valu is None:
            raise ValueError(
                "Unparsable card, fix it first with .verify('fix').")

        self._field_specifier = valu.group('keyword')

        if not hasattr(self, '_valuestring'):
            self._valuestring = valu.group('val')
        if not hasattr(self, '_value_modified'):
            self._value_modified = False

        return _str_to_num(translate(valu.group('val'), FIX_FP_TABLE2, ' '))

    def strvalue(self):
        """
        Method to extract the field specifier and value from the card
        image.  This is what is reported to the user when requesting
        the value of the `Card` using either an integer index or the
        card key without any field specifier.
        """

        mo = self.field_specifier_NFSC_image_RE.search(self.cardimage)
        return self.cardimage[mo.start():mo.end()]

    def _fix_value(self, input):
        """Fix the card image for fixable non-standard compliance."""

        _val_str = None

        if input is None:
            tmp = self._get_value_comment_string()

            try:
                slash_loc = tmp.index("/")
            except:
                slash_loc = len(tmp)

            self._err_text = 'Illegal value %s' % tmp[:slash_loc]
            self._fixable = False
            raise ValueError(self._err_text)
        else:
            self._valuestring = translate(input.group('val'), FIX_FP_TABLE,
                                          ' ')
            self._update_cardimage()

    def _format_key(self):
        if hasattr(self, '_key') or hasattr(self, '_cardimage'):
            return '%-8s' % super(RecordValuedKeywordCard, self).key
        else:
            return ' ' * 8

    def _check_key(self, key):
        """
        Verify the keyword to be FITS standard and that it matches the
        standard for record-valued keyword cards.
        """

        if '.' in key:
            keyword, field_specifier = key.split('.', 1)
        else:
            keyword, field_specifier = key, None

        super(RecordValuedKeywordCard, self)._check_key(keyword)

        if field_specifier:
            if not re.match(self.field_specifier_s, key):
                self._err_text = 'Illegal keyword name %s' % repr(key)
                # TODO: Maybe fix by treating as normal card and not RVKC?
                self._fixable = False
                raise ValueError(self._err_text)


    def _check(self, option='ignore'):
        """Verify the card image with the specified `option`."""

        self._err_text = ''
        self._fix_text = ''
        self._fixable = True

        if option == 'ignore':
            return
        elif option == 'parse':
            return self.keyword_NFSC_val_comm_RE.match(
                    self._get_value_comment_string())
        else:
            # verify the equal sign position
            if self.cardimage.find('=') != 8:
                if option in ['exception', 'warn']:
                    self._err_text = \
                        'Card image is not FITS standard (equal sign not at ' \
                        'column 8).'
                    raise ValueError(self._err_text + '\n%s' % self.cardimage)
                elif option in ['fix', 'silentfix']:
                    result = self._check('parse')
                    self._fix_value(result)

                    if option == 'fix':
                        self._fix_text = \
                           'Fixed card to be FITS standard. : %s' % self.key

            # verify the key
            self._check_key(self.key)

            # verify the value
            result = \
              self.keyword_val_comm_RE.match(self._get_value_comment_string())

            if result is not None:
                return result
            else:
                if option in ['fix', 'silentfix']:
                    result = self._check('parse')
                    self._fix_value(result)

                    if option == 'fix':
                        self._fix_text = \
                              'Fixed card to be FITS standard.: %s' % self.key
                else:
                    self._err_text = \
                        'Card image is not FITS standard (unparsable value ' \
                        'string).'
                    raise ValueError(self._err_text + '\n%s' % self.cardimage)

            # verify the comment (string), it is never fixable
            if result is not None:
                _str = result.group('comm')
                if _str is not None:
                    self._check_text(_str)


class CardList(list):
    """FITS header card list class."""

    def __init__(self, cards=[], keylist=None):
        """
        Construct the `CardList` object from a list of `Card` objects.

        Parameters
        ----------
        cards
            A list of `Card` objects.
        """

        super(CardList, self).__init__(cards)

        # if the key list is not supplied (as in reading in the FITS file),
        # it will be constructed from the card list.
        if keylist is None:
            self._keys = []
            for c in self:
                if isinstance(c, RecordValuedKeywordCard):
                    self._keys.append(c.raw.key)
                else:
                    self._keys.append(c.key)
        else:
            self._keys = keylist

        # find out how many blank cards are *directly* before the END card
        self._blanks = 0
        self.count_blanks()
        self._mod = False

    def __contains__(self, key):
        return upper_key(key) in self._keys

    def __getitem__(self, key):
        """Get a `Card` by indexing or by the keyword name."""

        if isinstance(key, slice):
            return CardList(super(CardList, self).__getitem__(key),
                            self._keys[key])
        elif isinstance(key, basestring) and self._has_filter_char(key):
            return self.filter_list(key)
        else:
            idx = self.index_of(key)
            return super(CardList, self).__getitem__(idx)

    def __setitem__(self, key, value):
        """Set a `Card` by indexing or by the keyword name."""

        if isinstance(value, Card):
            idx = self.index_of(key)

            # only set if the value is different from the old one
            if str(self[idx]) != str(value):
                super(CardList, self).__setitem__(idx, value)
                if isinstance(value, RecordValuedKeywordCard):
                    self._keys[idx] = value.raw.key.upper()
                else:
                    self._keys[idx] = value.key.upper()
                self.count_blanks()
                self._mod = True
        else:
            raise ValueError('%s is not a Card' % str(value))

    def __delitem__(self, key):
        """Delete a `Card` from the `CardList`."""

        if self._has_filter_char(key):
            cardlist = self.filter_list(key)

            if len(cardlist) == 0:
                raise KeyError("Keyword '%s' not found." % key)

            for card in cardlist:
                key = card.key
                del self[key]
        else:
            idx = self.index_of(key)
            super(CardList, self).__delitem__(idx)
            del self._keys[idx]  # update the keylist
            self.count_blanks()
            self._mod = True

    def __getslice__(self, start, end):
        return self[slice(start, end)]

    def __repr__(self):
        """Format a list of cards into a string."""

        return ''.join(map(str, self))

    def __str__(self):
        """Format a list of cards into a printable string."""
        return '\n'.join(map(str, self))

    @property
    def _mod(self):
        """Has this card list been modified since last write"""

        mod = self.__dict__.get('_mod', False)
        if not mod:
            # See if any of the cards were directly modified
            for card in self:
                if card._modified:
                    self.__dict__['_mod'] = True
                    return True
        return mod

    @_mod.setter
    def _mod(self, value):
        self.__dict__['_mod'] = value
        # If the card list is explicitly set as 'not modified' then make sure
        # the same applies to its underlying cards
        if not value:
            for card in self:
                card._modified = False

    def copy(self):
        """Make a (deep)copy of the `CardList`."""

        return CardList([create_card_from_string(str(c)) for c in self])

    def keys(self):
        """
        Return a list of all keywords from the `CardList`.

        Keywords include ``field_specifier`` for
        `RecordValuedKeywordCard` objects.
        """

        return [card.key for card in self]

    def values(self):
        """
        Return a list of the values of all cards in the `CardList`.

        For `RecordValuedKeywordCard` objects, the value returned is
        the floating point value, exclusive of the
        ``field_specifier``.
        """

        return [c.value for c in self]

    def insert(self, pos, card, useblanks=True):
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

        if isinstance(card, Card):
            super(CardList, self).insert(pos, card)
            if isinstance(card, RecordValuedKeywordCard):
                self._keys.insert(pos, card.raw.key.upper())
            else:
                self._keys.insert(pos, card.key.upper())  # update the keylist
            self.count_blanks()
            if useblanks:
                self._use_blanks(card._ncards())

            self.count_blanks()
            self._mod = True
        else:
            raise ValueError('%s is not a Card' % str(card))

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

        if isinstance(card, Card):
            nc = len(self) - self._blanks
            idx = nc - 1
            if not bottom:
                # locate last non-commentary card
                for idx in range(nc - 1, -1, -1):
                    if self[idx].key not in Card._commentary_keys:
                        break

            super(CardList, self).insert(idx + 1, card)
            if isinstance(card, RecordValuedKeywordCard):
                self._keys.insert(idx + 1, card.raw.key.upper())
            else:
                self._keys.insert(idx + 1, card.key.upper())
            if useblanks:
                self._use_blanks(card._ncards())
            self.count_blanks()
            self._mod = True
        else:
            raise ValueError("%s is not a Card" % str(card))

    def index_of(self, key, backward=False):
        """
        Get the index of a keyword in the `CardList`.

        Parameters
        ----------
        key : str or int
            The keyword name (a string) or the index (an integer).

        backward : bool, optional
            When `True`, search the index from the ``END``, i.e.,
            backward.

        Returns
        -------
        index : int
            The index of the `Card` with the given keyword.
        """

        if _is_int(key):
            return key
        elif isinstance(key, basestring):
            _key = key.strip().upper()
            if _key[:8] == 'HIERARCH':
                _key = _key[8:].strip()
            keys = self._keys
            if backward:
                # We have to search backwards through they key list
                keys = reversed(keys)
            try:
                idx = keys.index(_key)
            except ValueError:
                reqkey = RecordValuedKeywordCard.valid_key_value(key)
                idx = 0

                while reqkey:
                    try:
                        jdx = keys[idx:].index(reqkey[0].upper())
                        idx += jdx
                        card = self[idx]
                        if isinstance(card, RecordValuedKeywordCard) and \
                           reqkey[1] == card.field_specifier:
                            break
                    except ValueError:
                        raise KeyError('Keyword %s not found.' % repr(key))

                    idx += 1
                else:
                    raise KeyError('Keyword %s not found.' % repr(key))

            if backward:
                idx = len(keys) - idx - 1
            return idx
        else:
            raise KeyError('Illegal key data type %s' % type(key))

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

        out_cl = CardList()

        mykey = upper_key(key)
        re_str = mykey.replace('*', '\w*') + '$'
        re_str = re_str.replace('?', '\w')
        re_str = re_str.replace('...', '\S*')
        match_RE = re.compile(re_str)

        for card in self:
            match_str = card.key

            if match_RE.match(match_str):
                out_cl.append(card)

        return out_cl
    filterList = filter_list # For API backwards-compatibility

    def count_blanks(self):
        """
        Returns how many blank cards are *directly* before the ``END``
        card.
        """

        blank = ' ' * Card.length
        for idx in range(1, len(self)):
            if str(self[-idx]) != blank:
                self._blanks = idx - 1
                break

    def _has_filter_char(self, key):
        """
        Return `True` if the input key contains one of the special filtering
        characters (``*``, ``?``, or ...).
        """

        if isinstance(key, basestring) and \
           (key.endswith('...') or key.find('*') > 0 or key.find('?') > 0):
            return True
        else:
            return False

    def _pos_insert(self, card, before, after, useblanks=True):
        """
        Insert a `Card` to the location specified by before or after.

        The argument `before` takes precedence over `after` if both
        specified.  They can be either a keyword name or index.
        """

        if before is not None:
            loc = self.index_of(before)
            self.insert(loc, card, useblanks=useblanks)
        elif after is not None:
            loc = self.index_of(after)
            self.insert(loc + 1, card, useblanks=useblanks)

    def _use_blanks(self, how_many):
        if self._blanks > 0:
            for idx in range(min(self._blanks, how_many)):
                del self[-1] # it also delete the keylist item


def create_card(key='', value='', comment=''):
    return RecordValuedKeywordCard.create(key, value, comment)
create_card.__doc__ = RecordValuedKeywordCard.create.__doc__
# For API backwards-compatibility
createCard = deprecated(name='createCard',
                        alternative='create_card()')(create_card)


def create_card_from_string(input):
    return RecordValuedKeywordCard.fromstring(input)
create_card_from_string.__doc__ = RecordValuedKeywordCard.fromstring.__doc__
# For API backwards-compat
createCardFromString = \
        deprecated(name='createCardFromString',
                   alternative='fromstring()')(create_card_from_string)

def upper_key(key):
    return RecordValuedKeywordCard.upper_key(key)
upper_key.__doc__ = RecordValuedKeywordCard.upper_key.__doc__
# For API backwards-compat
upperKey = deprecated(name='upperKey', alternative='upper_key()')(upper_key)


class _HierarchCard(Card):
    """
    Cards begins with ``HIERARCH`` which allows keyword name longer than 8
    characters.
    """

    @property
    def key(self):
        """Returns the keyword name parsed from the card image."""

        if not hasattr(self, '_key'):
            head = self._get_key_string()
            self._key = head.strip()
        return self._key

    @key.setter
    def key(self, value):
        Card.key.fset(self, value)

    def _format_key(self):
        if hasattr(self, '_key') or hasattr(self, '_cardimage'):
            return 'HIERARCH %s ' % self.key
        else:
            return ' ' * 8

    def _format_value(self):
        val_str = super(_HierarchCard, self)._format_value()
        # conserve space for HIERARCH cards
        return val_str.strip()


    def _verify(self, option='warn'):
        """No verification (for now)."""

        return _ErrList([])


class _ContinueCard(Card):
    """
    Cards having more than one 80-char "physical" cards, the cards after
    the first one must start with ``CONTINUE`` and the whole card must have
    string value.
    """

    def __repr__(self):
        """Format a list of cards into a printable string."""

        output = []
        for idx in range(len(self.cardimage) // 80):
            output.append(self.cardimage[idx*80:(idx+1)*80])
        return '\n'.join(output)

    def _iter_cards(self):
        ncards = self._ncards()
        for idx in range(ncards):
            # take each 80-char card as a regular card and use its methods.
            card = Card.fromstring(self.cardimage[idx*80:(idx+1)*80])
            if idx > 0 and card.key != 'CONTINUE':
                raise ValueError('Long card image must have CONTINUE cards '
                                 'after the first card.')
            if not isinstance(card.value, str):
                raise ValueError(
                    'Cards with CONTINUE must have string value.')
            yield card

    def _extract_value(self):
        """Extract the keyword value from the card image."""

        output = []
        for card in self._iter_cards():
            val = card.value.rstrip().replace("''", "'")
            # drop the ending "&"
            if val and val[-1] == '&':
                val = val[:-1]
            output.append(val)
        return ''.join(output).rstrip()

    def _extract_comment(self):
        """Extract the comment from the card image."""

        output = []
        for card in self._iter_cards():
            comm = card.comment
            if isinstance(comm, str) and comm != '':
                output.append(comm.rstrip() + ' ')
        return ''.join(output).rstrip()

    def _breakup_strings(self):
        """
        Break up long string value/comment into ``CONTINUE`` cards.
        This is a primitive implementation: it will put the value
        string in one block and the comment string in another.  Also,
        it does not break at the blank space between words.  So it may
        not look pretty.
        """

        val_len = 67
        comm_len = 64
        output = []

        # do the value string
        valfmt = "'%-s&'"
        val = self.value.replace("'", "''")
        val_list = self._words_group(val, val_len)
        for idx in range(len(val_list)):
            if idx == 0:
                headstr = "%-8s= " % self.key
            else:
                headstr = "CONTINUE  "
            valstr = valfmt % val_list[idx]
            output.append('%-80s' % (headstr + valstr))

        # do the comment string
        if self.comment is None:
            comm = ''
        else:
            comm = self.comment
        commfmt = "%-s"
        if not comm == '':
            comm_list = self._words_group(comm, comm_len)
            for idx in comm_list:
                commstr = "CONTINUE  '&' / " + commfmt % idx
                output.append('%-80s' % commstr)

        return ''.join(output)

    def _words_group(self, input, strlen):
        """
        Split a long string into parts where each part is no longer
        than `strlen` and no word is cut into two pieces.  But if
        there is one single word which is longer than `strlen`, then
        it will be split in the middle of the word.
        """

        lst = []
        _nblanks = input.count(' ')
        nmax = max(_nblanks, len(input)//strlen+1)
        arr = np.fromstring((input + ' '), dtype=(np.bytes_, 1))

        # locations of the blanks
        blank_loc = np.nonzero(arr == np.bytes_(' '))[0]
        offset = 0
        xoffset = 0
        for idx in range(nmax):
            try:
                loc = np.nonzero(blank_loc >= strlen+offset)[0][0]
                offset = blank_loc[loc-1] + 1
                if loc == 0:
                    offset = -1
            except:
                offset = len(input)

            # check for one word longer than strlen, break in the middle
            if offset <= xoffset:
                offset = xoffset + strlen

            # collect the pieces in a list
            tmp = input[xoffset:offset]
            lst.append(tmp)
            if len(input) == offset:
                break
            xoffset = offset

        return lst


def _value_to_string(value):
    """Converts a card value to its appropriate string representation as
    defined by the FITS format.
    """

    # string value should occupies at least 8 columns, unless it is
    # a null string
    if isinstance(value, str):
        if value == '':
            return "''"
        else:
            exp_val_str = value.replace("'", "''")
            val_str = "'%-8s'" % exp_val_str
            return '%-20s' % val_str

    # must be before int checking since bool is also int
    elif isinstance(value, (bool, np.bool_)):
        return '%20s' % repr(value)[0] # T or F

    elif _is_int(value):
        return '%20d' % value

    # XXX need to consider platform dependence of the format (e.g. E-009 vs. E-09)
    elif isinstance(value, (float, np.floating)):
        return '%20s' % _float_format(value)

    elif isinstance(value, (complex, np.complexfloating)):
        val_str = '(%s, %s)' % (_float_format(value.real),
                                _float_format(value.imag))
        return '%20s' % val_str

    elif isinstance(value, Undefined):
        return ''
    else:
        return ''


def _float_format(value):
    """Format a floating number to make sure it gets the decimal point."""

    value_str = '%.16G' % value
    if '.' not in value_str and 'E' not in value_str:
        value_str += '.0'

    # Limit the value string to at most 20 characters.
    str_len = len(value_str)

    if str_len > 20:
        idx = value_str.find('E')

        if idx < 0:
            value_str = value_str[:20]
        else:
            value_str = value_str[:20-(str_len-idx)] + value_str[idx:]

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
            return input + ' ' * (Card.length-strlen)

    # minimum length is 80
    else:
        strlen = _len % Card.length
        return input + ' ' * (Card.length-strlen)


