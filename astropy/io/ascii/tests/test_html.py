# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``HTML``
reader/writer and aims to document its functionality.

Requires `BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_
to be installed.
"""

from .. import html
from .. import core

from ....tests.helper import pytest
from .common import (raises, assert_equal, assert_almost_equal,
                     assert_true, setup_function, teardown_function)

from itertools import izip

# Check to see if the BeautifulSoup dependency is present.

try:
    from bs4 import BeautifulSoup
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False

@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_soupstring():
    """
    Test to make sure the class SoupString behaves properly.
    """
    
    soup = BeautifulSoup('<html><body><p>foo</p></body></html>')
    soup_str = html.SoupString(soup)
    assert isinstance(soup_str, str)
    assert isinstance(soup_str, html.SoupString)
    assert soup_str == '<html><body><p>foo</p></body></html>'
    assert soup_str.soup is soup

def test_listwriter():
    """
    Test to make sure the class ListWriter behaves properly.
    """
    
    lst = []
    writer = html.ListWriter(lst)

    for i in range(5):
        writer.write(i)
    for ch in 'abcde':
        writer.write(ch)

    assert lst == [0, 1, 2, 3, 4, 'a', 'b', 'c', 'd', 'e']

@pytest.mark.skipif('HAS_BEAUTIFUL_SOUP')
def test_htmlinputter_no_bs4():
    """
    This should return an InconsistentTableError if BeautifulSoup
    is not installed.
    """

    inputter = html.HTMLInputter()
    with pytest.raises(core.InconsistentTableError):
        inputter.process_lines([])
    
@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_htmlinputter():
    """
    Test to ensure that HTMLInputter correctly converts input
    into a list of valid SoupString objects.
    """

    inputter = html.HTMLInputter()
    
    lines = ['<html>', '<body>', '<a href="http://www.astropy.org/">Link',
             '</a><br></br>', '</body>', '<!--Comment-->', '</html>']
    soup_list = inputter.process_lines(lines)
    expected_soups = ['<body><a href="http://www.astropy.org/">Link' \
        '</a><br/></body>', '<a href="http://www.astropy.org/">Link</a>',
        'Link', '<br/>']
    for soup, expected_soup in izip(soup_list, expected_soups):
        assert str(soup).replace('\n', '') == expected_soup
