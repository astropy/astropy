# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
from datetime import datetime
from xml.dom.minidom import parse
import re
import textwrap

# Third-party
from astropy import time as atime
from astropy.utils.console import color_print, _color_text
from six.moves.urllib.request import urlopen

__all__ = ['horoscope']

def get_sign(dt):
    """
    """
    if ((int(dt.month)==12 and int(dt.day) >= 22)or(int(dt.month)==1 and int(dt.day)<= 19)):
        zodiac_sign = "capricorn"
    elif ((int(dt.month)==1 and int(dt.day) >= 20)or(int(dt.month)==2 and int(dt.day)<= 17)):
        zodiac_sign = "aquarius"
    elif ((int(dt.month)==2 and int(dt.day) >= 18)or(int(dt.month)==3 and int(dt.day)<= 19)):
        zodiac_sign = "pisces"
    elif ((int(dt.month)==3 and int(dt.day) >= 20)or(int(dt.month)==4 and int(dt.day)<= 19)):
        zodiac_sign = "aries"
    elif ((int(dt.month)==4 and int(dt.day) >= 20)or(int(dt.month)==5 and int(dt.day)<= 20)):
        zodiac_sign = "taurus"
    elif ((int(dt.month)==5 and int(dt.day) >= 21)or(int(dt.month)==6 and int(dt.day)<= 20)):
        zodiac_sign = "gemini"
    elif ((int(dt.month)==6 and int(dt.day) >= 21)or(int(dt.month)==7 and int(dt.day)<= 22)):
        zodiac_sign = "cancer"
    elif ((int(dt.month)==7 and int(dt.day) >= 23)or(int(dt.month)==8 and int(dt.day)<= 22)):
        zodiac_sign = "leo"
    elif ((int(dt.month)==8 and int(dt.day) >= 23)or(int(dt.month)==9 and int(dt.day)<= 22)):
        zodiac_sign = "virgo"
    elif ((int(dt.month)==9 and int(dt.day) >= 23)or(int(dt.month)==10 and int(dt.day)<= 22)):
        zodiac_sign = "libra"
    elif ((int(dt.month)==10 and int(dt.day) >= 23)or(int(dt.month)==11 and int(dt.day)<= 21)):
        zodiac_sign = "scorpio"
    elif ((int(dt.month)==11 and int(dt.day) >= 22)or(int(dt.month)==12 and int(dt.day)<= 21)):
        zodiac_sign = "sagittarius"

    return zodiac_sign

def horoscope(birthday):
    """
    Enter your birthday as an `astropy.time.Time` object and
    receive a mystical horoscope about things to come.

    Parameter
    ---------
    birthday : `astropy.time.Time`
        Your birthday as a `datetime.datetime` or `astropy.time.Time` object.

    Returns
    -------
    Infinite wisdom, condensed into astrologically precise prose.

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
    zodiac_sign = get_sign(birthday.to_datetime())
    url = "http://www.findyourfate.com/rss/dailyhoroscope-feed.asp?sign={sign}"

    with urlopen(url.format(sign=zodiac_sign.capitalize())) as f:
        doc = parse(f)

    print(url.format(sign=zodiac_sign.capitalize()))

    try:
        item = doc.getElementsByTagName('item')[0]
    except IndexError:
        raise ValueError("Invalid response from celestial gods (failed to load horoscope).")

    try:
        desc = item.getElementsByTagName('description')[0].childNodes[0].nodeValue
    except (IndexError, AttributeError):
        raise ValueError("Invalid response from celestial gods (failed to load horoscope).")

    print("*"*79)
    color_print("Horoscope for {} on {}:".format(zodiac_sign.capitalize(), today.strftime("%Y-%m-%d")),
                'green')
    print("*"*79)
    for block in textwrap.wrap(desc, 79):
        split_block = block.split()
        for i,word in enumerate(split_block):
            for re_word in special_words.keys():
                match = re.search(re_word, word)
                if match is None:
                    continue
                split_block[i] = _color_text(match.groups()[0], special_words[re_word])
        print(" ".join(split_block))
