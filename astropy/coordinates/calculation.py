# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


# Standard library
import re
import textwrap
from datetime import datetime
from xml.dom.minidom import parse
from urllib.request import urlopen

# Third-party
from .. import time as atime
from ..utils.console import color_print, _color_text
from . import get_sun

__all__ = []


class HumanError(ValueError): pass


class CelestialError(ValueError): pass


def get_sign(dt):
    """
    """
    if ((int(dt.month) == 12 and int(dt.day) >= 22)or(int(dt.month) == 1 and int(dt.day) <= 19)):
        zodiac_sign = "capricorn"
    elif ((int(dt.month) == 1 and int(dt.day) >= 20)or(int(dt.month) == 2 and int(dt.day) <= 17)):
        zodiac_sign = "aquarius"
    elif ((int(dt.month) == 2 and int(dt.day) >= 18)or(int(dt.month) == 3 and int(dt.day) <= 19)):
        zodiac_sign = "pisces"
    elif ((int(dt.month) == 3 and int(dt.day) >= 20)or(int(dt.month) == 4 and int(dt.day) <= 19)):
        zodiac_sign = "aries"
    elif ((int(dt.month) == 4 and int(dt.day) >= 20)or(int(dt.month) == 5 and int(dt.day) <= 20)):
        zodiac_sign = "taurus"
    elif ((int(dt.month) == 5 and int(dt.day) >= 21)or(int(dt.month) == 6 and int(dt.day) <= 20)):
        zodiac_sign = "gemini"
    elif ((int(dt.month) == 6 and int(dt.day) >= 21)or(int(dt.month) == 7 and int(dt.day) <= 22)):
        zodiac_sign = "cancer"
    elif ((int(dt.month) == 7 and int(dt.day) >= 23)or(int(dt.month) == 8 and int(dt.day) <= 22)):
        zodiac_sign = "leo"
    elif ((int(dt.month) == 8 and int(dt.day) >= 23)or(int(dt.month) == 9 and int(dt.day) <= 22)):
        zodiac_sign = "virgo"
    elif ((int(dt.month) == 9 and int(dt.day) >= 23)or(int(dt.month) == 10 and int(dt.day) <= 22)):
        zodiac_sign = "libra"
    elif ((int(dt.month) == 10 and int(dt.day) >= 23)or(int(dt.month) == 11 and int(dt.day) <= 21)):
        zodiac_sign = "scorpio"
    elif ((int(dt.month) == 11 and int(dt.day) >= 22)or(int(dt.month) == 12 and int(dt.day) <= 21)):
        zodiac_sign = "sagittarius"

    return zodiac_sign


_VALID_SIGNS = ["capricorn", "aquarius", "pisces", "aries", "taurus", "gemini",
                "cancer", "leo", "virgo", "libra", "scorpio", "sagittarius"]
# Some of the constellation names map to different astrological "sign names".
# Astrologers really needs to talk to the IAU...
_CONST_TO_SIGNS = {'capricornus': 'capricorn', 'scorpius': 'scorpio'}


def horoscope(birthday, corrected=True):
    """
    Enter your birthday as an `astropy.time.Time` object and
    receive a mystical horoscope about things to come.

    Parameter
    ---------
    birthday : `astropy.time.Time`
        Your birthday as a `datetime.datetime` or `astropy.time.Time` object.
    corrected : bool
        Whether to account for the precession of the Earth instead of using the
        ancient Greek dates for the signs.  After all, you do want your *real*
        horoscope, not a cheap inaccurate approximation, right?

    Returns
    -------
    Infinite wisdom, condensed into astrologically precise prose.

    Notes
    -----
    This function was implemented on April 1.  Take note of that date.
    """

    special_words = {
        '([sS]tar[s^ ]*)': 'yellow',
        '([yY]ou[^ ]*)': 'magenta',
        '([pP]lay[^ ]*)': 'blue',
        '([hH]eart)': 'red',
        '([fF]ate)': 'lightgreen',
    }

    birthday = atime.Time(birthday)
    today = datetime.now()
    if corrected:
        zodiac_sign = get_sun(birthday).get_constellation().lower()
        zodiac_sign = _CONST_TO_SIGNS.get(zodiac_sign, zodiac_sign)
        if zodiac_sign not in _VALID_SIGNS:
            raise HumanError('On your birthday the sun was in {}, which is not '
                             'a sign of the zodiac.  You must not exist.  Or '
                             'maybe you can settle for '
                             'corrected=False.'.format(zodiac_sign.title()))
    else:
        zodiac_sign = get_sign(birthday.to_datetime())
    url = "http://www.findyourfate.com/rss/dailyhoroscope-feed.php?sign={sign}&id=45"

    f = urlopen(url.format(sign=zodiac_sign.capitalize()))
    try:  # urlopen in py2 is not a decorator
        doc = parse(f)
        item = doc.getElementsByTagName('item')[0]
        desc = item.getElementsByTagName('description')[0].childNodes[0].nodeValue
    except Exception:
        raise CelestialError("Invalid response from celestial gods (failed to load horoscope).")
    finally:
        f.close()

    print("*"*79)
    color_print("Horoscope for {} on {}:".format(zodiac_sign.capitalize(), today.strftime("%Y-%m-%d")),
                'green')
    print("*"*79)
    for block in textwrap.wrap(desc, 79):
        split_block = block.split()
        for i, word in enumerate(split_block):
            for re_word in special_words.keys():
                match = re.search(re_word, word)
                if match is None:
                    continue
                split_block[i] = _color_text(match.groups()[0], special_words[re_word])
        print(" ".join(split_block))


def inject_horoscope():
    import astropy
    astropy._yourfuture = horoscope


inject_horoscope()
