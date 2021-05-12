# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


# Standard library
import re
import textwrap
import warnings
from datetime import datetime
from urllib.request import urlopen, Request

# Third-party
from astropy import time as atime
from astropy.utils.console import color_print, _color_text
from . import get_sun

__all__ = []


class HumanError(ValueError):
    pass


class CelestialError(ValueError):
    pass


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

_ZODIAC = ((1900, "rat"), (1901, "ox"), (1902, "tiger"),
           (1903, "rabbit"), (1904, "dragon"), (1905, "snake"),
           (1906, "horse"), (1907, "goat"), (1908, "monkey"),
           (1909, "rooster"), (1910, "dog"), (1911, "pig"))


# https://stackoverflow.com/questions/12791871/chinese-zodiac-python-program
def _get_zodiac(yr):
    return _ZODIAC[(yr - _ZODIAC[0][0]) % 12][1]


def horoscope(birthday, corrected=True, chinese=False):
    """
    Enter your birthday as an `astropy.time.Time` object and
    receive a mystical horoscope about things to come.

    Parameters
    ----------
    birthday : `astropy.time.Time` or str
        Your birthday as a `datetime.datetime` or `astropy.time.Time` object
        or "YYYY-MM-DD"string.
    corrected : bool
        Whether to account for the precession of the Earth instead of using the
        ancient Greek dates for the signs.  After all, you do want your *real*
        horoscope, not a cheap inaccurate approximation, right?

    chinese : bool
        Chinese annual zodiac wisdom instead of Western one.

    Returns
    -------
    Infinite wisdom, condensed into astrologically precise prose.

    Notes
    -----
    This function was implemented on April 1.  Take note of that date.
    """
    from bs4 import BeautifulSoup

    today = datetime.now()
    err_msg = "Invalid response from celestial gods (failed to load horoscope)."
    headers = {'User-Agent': 'foo/bar'}

    special_words = {
        '([sS]tar[s^ ]*)': 'yellow',
        '([yY]ou[^ ]*)': 'magenta',
        '([pP]lay[^ ]*)': 'blue',
        '([hH]eart)': 'red',
        '([fF]ate)': 'lightgreen',
    }

    if isinstance(birthday, str):
        birthday = datetime.strptime(birthday, '%Y-%m-%d')

    if chinese:
        # TODO: Make this more accurate by using the actual date, not just year
        # Might need third-party tool like https://pypi.org/project/lunardate
        zodiac_sign = _get_zodiac(birthday.year)
        url = ('https://www.horoscope.com/us/horoscopes/yearly/'
               '{}-chinese-horoscope-{}.aspx'.format(today.year, zodiac_sign))
        summ_title_sfx = f'in {today.year}'

        try:
            res = Request(url, headers=headers)
            with urlopen(res) as f:
                try:
                    doc = BeautifulSoup(f, 'html.parser')
                    # TODO: Also include Love, Family & Friends, Work, Money, More?
                    item = doc.find(id='overview')
                    desc = item.getText()
                except Exception:
                    raise CelestialError(err_msg)
        except Exception:
            raise CelestialError(err_msg)

    else:
        birthday = atime.Time(birthday)

        if corrected:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')  # Ignore ErfaWarning
                zodiac_sign = get_sun(birthday).get_constellation().lower()
            zodiac_sign = _CONST_TO_SIGNS.get(zodiac_sign, zodiac_sign)
            if zodiac_sign not in _VALID_SIGNS:
                raise HumanError('On your birthday the sun was in {}, which is not '
                                 'a sign of the zodiac.  You must not exist.  Or '
                                 'maybe you can settle for '
                                 'corrected=False.'.format(zodiac_sign.title()))
        else:
            zodiac_sign = get_sign(birthday.to_datetime())
        url = f"http://www.astrology.com/us/horoscope/daily-overview.aspx?sign={zodiac_sign}"
        summ_title_sfx = f"on {today.strftime('%Y-%m-%d')}"

        res = Request(url, headers=headers)
        with urlopen(res) as f:
            try:
                doc = BeautifulSoup(f, 'html.parser')
                item = doc.find('div', {'id': 'content'})
                desc = item.getText()
            except Exception:
                raise CelestialError(err_msg)

    print("*"*79)
    color_print(f"Horoscope for {zodiac_sign.capitalize()} {summ_title_sfx}:",
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
