import os
import warnings

try:
    from collections import MutableMapping as __HEADERBASE
except ImportError:
    from UserDict import DictMixin
    class __HEADERBASE(DictMixin, object): # object to make a newstyle class
        pass

from pyfits.card import (Card, CardList, RecordValuedKeywordCard,
                         _ContinueCard, _HierarchCard, create_card,
                         create_card_from_string, upper_key, _pad)
from pyfits.util import BLOCK_SIZE, deprecated, _pad_length


class Header(__HEADERBASE):
    """
    FITS header class.

    The purpose of this class is to present the header like a
    dictionary as opposed to a list of cards.

    The attribute `ascard` supplies the header like a list of cards.

    The header class uses the card's keyword as the dictionary key and
    the cards value is the dictionary value.

    The `has_key`, `get`, and `keys` methods are implemented to
    provide the corresponding dictionary functionality.  The header
    may be indexed by keyword value and like a dictionary, the
    associated value will be returned.  When the header contains cards
    with duplicate keywords, only the value of the first card with the
    given keyword will be returned.

    The header may also be indexed by card list index number.  In that
    case, the value of the card at the given index in the card list
    will be returned.

    A delete method has been implemented to allow deletion from the
    header.  When `del` is called, all cards with the given keyword
    are deleted from the header.

    The `Header` class has an associated iterator class `_Header_iter`
    which will allow iteration over the unique keywords in the header
    dictionary.
    """

    # TODO: Allow the header to take a few other types of inputs, for example
    # a list of (key, value) tuples, (key, value, comment) tuples, or a dict
    # of either key: value or key: (value, comment) mappings.  This could all
    # be handled by the underlying CardList I suppose.
    def __init__(self, cards=[], txtfile=None):
        """
        Construct a `Header` from a `CardList` and/or text file.

        Parameters
        ----------
        cards : A list of `Card` objects, optional
            The cards to initialize the header with.

        txtfile : file path, file object or file-like object, optional
            Input ASCII header parameters file.
        """

        # populate the cardlist
        self.ascard = CardList(cards)

        if txtfile:
            # get the cards from the input ASCII file
            self.fromTxtFile(txtfile, not len(self.ascard))
        self._mod = False

    def __len__(self):
        return len(self.ascard)

    def __iter__(self):
        for card in self.ascard:
            yield card.key

    def __contains__(self, key):
        """
        Check for existence of a keyword.

        Parameters
        ----------
        key : str or int
           Keyword name.  If given an index, always returns 0.

        Returns
        -------
        has_key : bool
            Returns `True` if found, otherwise, `False`.
        """

        key = upper_key(key)
        if key[:8] == 'HIERARCH':
            key = key[8:].strip()
        return key in self.ascard
    has_key = deprecated(name='has_key',
                         alternative='`key in header` syntax')(__contains__)

    def __getitem__ (self, key):
        """
        Get a header keyword value.
        """

        card = self.ascard[key]

        if isinstance(card, RecordValuedKeywordCard) and \
           (not isinstance(key, basestring) or '.' not in key):
            return card.strvalue()
        elif isinstance(card, CardList):
            return card
        else:
            return card.value

    def __setitem__ (self, key, value):
        """
        Set a header keyword value.
        """

        self.ascard[key].value = value
        self._mod = 1

    def __delitem__(self, key):
        """
        Delete card(s) with the name `key`.
        """

        # delete ALL cards with the same keyword name
        if isinstance(key, basestring):
            while True:
                try:
                    del self.ascard[key]
                    self._mod = True
                except:
                    return

        # for integer key only delete once
        else:
            del self.ascard[key]
            self._mod = True

    def __repr__(self):
        return repr(self.ascard)

    def __str__(self):
        s = repr(self) + _pad('END')
        return s + _pad_length(len(s)) * ' '

    @classmethod
    def fromstring(cls, data):
        """
        Creates an HDU header from a byte string containing the entire header
        data.

        Parameters
        ----------
        data : str
           String containing the entire header.
        """

        if (len(data) % BLOCK_SIZE) != 0:
            raise ValueError('Header size is not multiple of %d: %d'
                             % (BLOCK_SIZE, len(data)))

        cards = []
        keys = []

        # Split the header into individual cards
        for idx in range(0, len(data), Card.length):
            card = create_card_from_string(data[idx:idx + Card.length])
            key = card.key

            if key == 'END':
                break
            else:
                cards.append(card)
                keys.append(key)

        # Deal with CONTINUE cards
        # if a long string has CONTINUE cards, the "Card" is considered
        # to be more than one 80-char "physical" cards.
        idx = len(cards)
        continueimg = []
        for card in reversed(cards):
            idx -= 1
            if idx != 0 and card.key == 'CONTINUE':
                continueimg.append(card._cardimage)
                del cards[idx]
            elif continueimg:
                continueimg.append(card._cardimage)
                continueimg = ''.join(reversed(continueimg))
                cards[idx] = _ContinueCard.fromstring(continueimg)
                continueimg = []

        return cls(CardList(cards, keylist=keys))

    def keys(self):
        """
        Return a list of keys with duplicates removed.
        """

        retval = []

        for key in self.ascard.keys():
            if not key in retval:
                retval.append(key)

        return retval

    def update(self, key, value, comment=None, before=None, after=None,
               savecomment=False):
        """
        Update one header card.

        If the keyword already exists, it's value and/or comment will
        be updated.  If it does not exist, a new card will be created
        and it will be placed before or after the specified location.
        If no `before` or `after` is specified, it will be appended at
        the end.

        Parameters
        ----------
        key : str
            keyword

        value : str
            value to be used for updating

        comment : str, optional
            to be used for updating, default=None.

        before : str or int, optional
            name of the keyword, or index of the `Card` before which
            the new card will be placed.  The argument `before` takes
            precedence over `after` if both specified.

        after : str or int, optional
            name of the keyword, or index of the `Card` after which
            the new card will be placed.

        savecomment : bool, optional
            When `True`, preserve the current comment for an existing
            keyword.  The argument `savecomment` takes precedence over
            `comment` if both specified.  If `comment` is not
            specified then the current comment will automatically be
            preserved.
        """

        keylist = RecordValuedKeywordCard.valid_key_value(key, value)

        if keylist:
            keyword = keylist[0] + '.' + keylist[1]
        else:
            keyword = key

        if keyword in self:
            j = self.ascard.index_of(keyword)
            if not savecomment and comment is not None:
                _comment = comment
            else:
                _comment = self.ascard[j].comment
            _card = create_card(key, value, _comment)
            if before is not None or after is not None:
                del self.ascard[j]
                self.ascard._pos_insert(_card, before=before, after=after)
            else:
                self.ascard[j] = _card
        elif before is not None or after is not None:
            _card = create_card(key, value, comment)
            self.ascard._pos_insert(_card, before=before, after=after)
        else:
            self.ascard.append(create_card(key, value, comment))

        self._mod = True

    def copy(self, strip=False):
        """
        Make a copy of the `Header`.

        Parameters
        ----------
        strip : bool, optional
           If True, strip any headers that are specific to one of the standard
           HDU types, so that this header can be used in a different HDU.
        """

        tmp = Header(self.ascard.copy())
        if strip:
            tmp._strip()
        return tmp

    @deprecated(alternative='the ascard attribute')
    def ascardlist(self):
        """
        Returns a `CardList` object.
        """

        return self.ascard

    def rename_key(self, oldkey, newkey, force=False):
        """
        Rename a card's keyword in the header.

        Parameters
        ----------
        oldkey : str or int
            old keyword

        newkey : str
            new keyword

        force : bool
            When `True`, if new key name already exists, force to have
            duplicate name.
        """

        oldkey = upper_key(oldkey)
        newkey = upper_key(newkey)

        if newkey == 'CONTINUE':
            raise ValueError('Can not rename to CONTINUE')

        if newkey in Card._commentary_keys or oldkey in Card._commentary_keys:
            if not (newkey in Card._commentary_keys and
                    oldkey in Card._commentary_keys):
                raise ValueError('Regular and commentary keys can not be '
                                 'renamed to each other.')
        elif (force == 0) and newkey in self:
            raise ValueError('Intended keyword %s already exists in header.'
                             % newkey)

        idx = self.ascard.index_of(oldkey)
        comment = self.ascard[idx].comment
        value = self.ascard[idx].value
        self.ascard[idx] = create_card(newkey, value, comment)

    def add_history(self, value, before=None, after=None):
        """
        Add a ``HISTORY`` card.

        Parameters
        ----------
        value : str
            history text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary('history', value, before=before, after=after)

    def add_comment(self, value, before=None, after=None):
        """
        Add a ``COMMENT`` card.

        Parameters
        ----------
        value : str
            text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary('comment', value, before=before, after=after)

    def add_blank(self, value='', before=None, after=None):
        """
        Add a blank card.

        Parameters
        ----------
        value : str, optional
            text to be added.

        before : str or int, optional
            same as in `Header.update`

        after : str or int, optional
            same as in `Header.update`
        """

        self._add_commentary(' ', value, before=before, after=after)

    def get_history(self):
        """
        Get all history cards as a list of string texts.
        """

        return [c for c in self.ascard if c.key == 'HISTORY']

    def get_comment(self):
        """
        Get all comment cards as a list of string texts.
        """

        return [c for c in self.ascard if c.key == 'COMMENT']

    def toTxtFile(self, fileobj, clobber=False):
        """
        Output the header parameters to a file in ASCII format.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            Output header parameters file.

        clobber : bool
            When `True`, overwrite the output file if it exists.
        """

        close_file = False

        # check if the output file already exists
        if isinstance(fileobj, basestring):
            if (os.path.exists(fileobj) and os.path.getsize(fileobj) != 0):
                if clobber:
                    warnings.warn("Overwriting existing file '%s'." % fileobj)
                    os.remove(fileobj)
                else:
                    raise IOError("File '%s' already exist." % fileobj)

            fileobj = open(fileobj, 'w')
            close_file = True

        lines = []   # lines to go out to the header parameters file

        # Add the card image for each card in the header to the lines list

        for j in range(len(self.ascard)):
            lines.append(str(self.ascard[j]) + '\n')

        # Write the header parameter lines out to the ASCII header
        # parameter file
        fileobj.writelines(lines)

        if close_file:
            fileobj.close()

    def fromTxtFile(self, fileobj, replace=False):
        """
        Input the header parameters from an ASCII file.

        The input header cards will be used to update the current
        header.  Therefore, when an input card key matches a card key
        that already exists in the header, that card will be updated
        in place.  Any input cards that do not already exist in the
        header will be added.  Cards will not be deleted from the
        header.

        Parameters
        ----------
        fileobj : file path, file object or file-like object
            Input header parameters file.

        replace : bool, optional
            When `True`, indicates that the entire header should be
            replaced with the contents of the ASCII file instead of
            just updating the current header.
        """

        close_file = False

        if isinstance(fileobj, basestring):
            fileobj = open(fileobj, 'r')
            close_file = True

        lines = fileobj.readlines()

        if close_file:
            fileobj.close()

        if len(self.ascard) > 0 and not replace:
            prevKey = 0
        else:
            if replace:
                self.ascard = CardList([])
            prevKey = 0

        for line in lines:
            card = Card.fromstring(line[:min(80, len(line)-1)])
            card.verify('silentfix')

            if card.key == 'SIMPLE':
                if self.get('EXTENSION'):
                    del self.ascard['EXTENSION']

                self.update(card.key, card.value, card.comment, before=0)
                prevKey = 0
            elif card.key == 'EXTENSION':
                if self.get('SIMPLE'):
                    del self.ascard['SIMPLE']

                self.update(card.key, card.value, card.comment, before=0)
                prevKey = 0
            elif card.key == 'HISTORY':
                if not replace:
                    items = self.items()
                    idx = 0

                    for item in items:
                        if item[0] == card.key and item[1] == card.value:
                            break
                        idx += 1

                    if idx == len(self.ascard):
                        self.add_history(card.value, after=prevKey)
                        prevKey += 1
                else:
                    self.add_history(card.value, after=prevKey)
                    prevKey += 1
            elif card.key == 'COMMENT':
                if not replace:
                    items = self.items()
                    idx = 0

                    for item in items:
                        if item[0] == card.key and item[1] == card.value:
                            break
                        idx += 1

                    if idx == len(self.ascard):
                        self.add_comment(card.value, after=prevKey)
                        prevKey += 1
                else:
                    self.add_comment(card.value, after=prevKey)
                    prevKey += 1
            elif card.key == '        ':
                if not replace:
                    items = self.items()
                    idx = 0

                    for item in items:
                        if item[0] == card.key and item[1] == card.value:
                            break
                        idx += 1

                    if idx == len(self.ascard):
                        self.add_blank(card.value, after=prevKey)
                        prevKey += 1
                else:
                    self.add_blank(card.value, after=prevKey)
                    prevKey += 1
            else:
                if isinstance(card, _HierarchCard):
                    prefix = 'hierarch '
                else:
                    prefix = ''

                self.update(prefix + card.key,
                                     card.value,
                                     card.comment,
                                     after=prevKey)
                prevKey += 1

    def _add_commentary(self, key, value, before=None, after=None):
        """
        Add a commentary card.

        If `before` and `after` are `None`, add to the last occurrence
        of cards of the same name (except blank card).  If there is no
        card (or blank card), append at the end.
        """

        new_card = Card(key, value)
        if before is not None or after is not None:
            self.ascard._pos_insert(new_card, before=before, after=after)
        else:
            if key[0] == ' ':
                useblanks = new_card.cardimage != ' '*80
                self.ascard.append(new_card, useblanks=useblanks, bottom=1)
            else:
                try:
                    _last = self.ascard.index_of(key, backward=1)
                    self.ascard.insert(_last+1, new_card)
                except:
                    self.ascard.append(new_card, bottom=1)

        self._mod = True

    def _strip(self):
        """
        Strip cards specific to a certain kind of header.

        Strip cards like ``SIMPLE``, ``BITPIX``, etc. so the rest of
        the header can be used to reconstruct another kind of header.
        """

        # TODO: Previously this only deleted some cards specific to an HDU if
        # _hdutype matched that type.  But it seemed simple enough to just
        # delete all desired cards anyways, and just ignore the KeyErrors if
        # they don't exist.
        # However, it might be desirable to make this extendable somehow--have
        # a way for HDU classes to specify some headers that are specific only
        # to that type, and should be removed otherwise.

        try:
            if 'NAXIS' in self:
                naxis = self['NAXIS']
            else:
                naxis = 0

            if 'TFIELDS' in self:
                tfields = self['TFIELDS']
            else:
                tfields = 0

            for idx in range(naxis):
                del self['NAXIS' + str(idx + 1)]

            for name in ('TFORM', 'TSCAL', 'TZERO', 'TNULL', 'TTYPE',
                         'TUNIT', 'TDISP', 'TDIM', 'THEAP', 'TBCOL'):
                for idx in range(tfields):
                    del self[name + str(idx + 1)]

            for name in ('SIMPLE', 'XTENSION', 'BITPIX', 'NAXIS', 'EXTEND',
                         'PCOUNT', 'GCOUNT', 'GROUPS', 'BSCALE', 'TFIELDS'):
                del self[name]
        except KeyError:
            pass
