# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible HTML table reader and writer.

html.py:
  Classes to read and write HTML tables

This requires `BeautifulSoup
        <http://www.crummy.com/software/BeautifulSoup/>`_ to be installed.
"""

##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import absolute_import, division, print_function

from . import core

from bs4 import BeautifulSoup
from bs4.element import Comment

class SoupString(str):
    """
    Allows for strings to hold BeautifulSoup data.
    """

    def __new__(cls, *args, **kwargs):
        return str.__new__(cls, *args, **kwargs)

    def __init__(self, val):
        self.soup = val

class HTMLInputter(core.BaseInputter):
    """
    Input lines of HTML in a valid form.
    """

    def process_lines(self, lines):
        soup = BeautifulSoup('\n'.join(lines))
        soup_list = []
        for x in soup.contents[0].descendants:
            if str(x).strip():
                soup_obj = SoupString(x)
                soup_list.append(soup_obj)
        return soup_list
        
class HTMLSplitter(core.BaseSplitter):
    """
    Split HTML table data.
    """

    def __call__(self, lines):
        """
        Return HTML data from lines as a generator.
        """
        data_found = False
        for line in lines:
            if not hasattr(line, 'soup'):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name in ('table', 'tr'):
                data_found = True
                header_elements = soup.find_all('th')
                if header_elements:
                    yield [el.string.strip() for el in header_elements]
                data_elements = soup.find_all('td')
                if data_elements:
                    yield [el.string.strip() for el in data_elements]
        if not data_found:
            raise core.InconsistentTableError('HTML tables must contain data '
                                              'in a <table> tag')

class HTMLHeader(core.BaseHeader):
    def start_line(self, lines):
        """
        Return the line number at which header data begins.
        """
        for i, line in enumerate(lines):
            if not hasattr(line, 'soup'):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'tr':
                return i
        raise core.InconsistentTableError('HTML tables must contain at least '
                                              'one <tr> tag')

class HTMLData(core.BaseData):
    def start_line(self, lines):
        """
        Return the line number at which table data begins.
        """
        for i, line in enumerate(lines):
            if not hasattr(line, 'soup'):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'tr' and soup.td is not None:
                if soup.th is not None:
                    raise core.InconsistentTableError('HTML tables cannot '
                                'have headings and data in the same row')
                return i
        raise core.InconsistentTableError('HTML tables must contain '
                                            'at least one data row')
    def end_line(self, lines):
        """
        Return the line number at which table data ends.
        """
        last_index = -1
        for i, line in enumerate(lines):
            if not hasattr(line, 'soup'):
                raise TypeError('HTML lines should be of type SoupString')
            soup = line.soup
            if soup.name == 'tr':
                last_index = i
        return last_index + 1

class HTML(core.BaseReader):
    _format_name = 'html'
    _io_registry_format_aliases = ['html']
    _io_registry_suffix = '.html'
    _description = 'HTML table'

    def __init__(self):
        core.BaseReader.__init__(self)
        self.inputter = HTMLInputter()
        self.header = HTMLHeader()
        self.data = HTMLData()
        self.header.splitter = HTMLSplitter()
        self.header.inputter = HTMLInputter()
        self.data.splitter = HTMLSplitter()
        self.data.inputter = HTMLInputter()
        self.data.header = self.header
        self.header.data = self.data
