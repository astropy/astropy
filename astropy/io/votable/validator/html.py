# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import with_statement, absolute_import

# STDLIB
import contextlib
import io
from math import ceil
import os
import re

# ASTROPY
from ....logger import log
from ....utils.xml.writer import XMLWriter, xml_escape
from .... import online_docs_root

# VO
from .. import exceptions

html_header = u"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html
        PUBLIC "-//W3C//DTD XHTML Basic 1.0//EN"
        "http://www.w3.org/TR/xhtml-basic/xhtml-basic10.dtd">
"""

default_style = u"""
body {
font-family: sans-serif
}
a {
text-decoration: none
}
.highlight {
color: red;
font-weight: bold;
text-decoration: underline;
}
.green { background-color: #ddffdd }
.red   { background-color: #ffdddd }
.yellow { background-color: #ffffdd }
tr:hover { background-color: #dddddd }
table {
        border-width: 1px;
        border-spacing: 0px;
        border-style: solid;
        border-color: gray;
        border-collapse: collapse;
        background-color: white;
        padding: 5px;
}
table th {
        border-width: 1px;
        padding: 5px;
        border-style: solid;
        border-color: gray;
}
table td {
        border-width: 1px;
        padding: 5px;
        border-style: solid;
        border-color: gray;
}
"""


@contextlib.contextmanager
def make_html_header(w):
    w.write(html_header)
    with w.tag(u'html', xmlns=u"http://www.w3.org/1999/xhtml", lang=u"en-US"):
        with w.tag(u'head'):
            w.element(u'title', u'VO Validation results')
            w.element(u'style', default_style)

            with w.tag(u'body'):
                yield


def write_source_line(w, line, nchar=0):
    part1 = xml_escape(line[:nchar].decode('utf-8'))
    char = xml_escape(line[nchar:nchar+1].decode('utf-8'))
    part2 = xml_escape(line[nchar+1:].decode('utf-8'))

    w.write(u'  ')
    w.write(part1)
    w.write(u'<span class="highlight">%s</span>' % char)
    w.write(part2)
    w.write(u'\n\n')


def write_warning(w, line, xml_lines):
    warning = exceptions.parse_vowarning(line)
    if not warning['is_something']:
        w.data(line)
    else:
        w.write(u'Line %d: ' % warning['nline'])
        if warning['warning']:
            w.write(u'<a href="%s/%s">%s</a>: ' % (
                online_docs_root, warning['doc_url'], warning['warning']))
        msg = warning['message']
        if not isinstance(warning['message'], unicode):
            msg = msg.decode('utf-8')
        w.write(xml_escape(msg))
        w.write(u'\n')
        try:
            write_source_line(
                w, xml_lines[warning['nline'] - 1], warning['nchar'])
        except IndexError as e:
            log.warn(
                'write_source_line failed for {} with {}'.format(line, str(e)))


def write_votlint_warning(w, line, xml_lines):
    match = re.search("(WARNING|ERROR|INFO) \(l.(?P<line>[0-9]+), c.(?P<column>[0-9]+)\): (?P<rest>.*)", line)
    if match:
        w.write(u'Line %d: %s\n' %
                (int(match.group('line')), xml_escape(match.group('rest'))))
        write_source_line(
            w, xml_lines[int(match.group('line')) - 1],
            int(match.group('column')) - 1)
    else:
        w.data(line)
        w.data('\n')


def write_result(result):
    if 'network_error' in result and result['network_error'] is not None:
        return

    xml = result.get_xml_content()
    xml_lines = xml.splitlines()

    path = os.path.join(result.get_dirpath(), 'index.html')

    with io.open(path, 'w', encoding='utf-8') as fd:
        w = XMLWriter(fd)
        with make_html_header(w):
            with w.tag(u'p'):
                with w.tag(u'a', href=u'vo.xml'):
                    w.data(result.url.decode('ascii'))
            w.element(u'hr')

            with w.tag(u'pre'):
                w._flush()
                for line in result['warnings']:
                    write_warning(w, line, xml_lines)

            if result['xmllint'] is False:
                w.element(u'hr')
                w.element(u'p', 'xmllint results:')
                content = result['xmllint_content']
                if not isinstance(content, unicode):
                    content = content.decode('ascii')
                content = content.replace(result.get_dirpath() + '/', '')
                with w.tag(u'pre'):
                    w.data(content)

            if 'votlint' in result:
                if result['votlint'] is False:
                    w.element(u'hr')
                    w.element(u'p', 'votlint results:')
                    content = result['votlint_content']
                    if not isinstance(content, unicode):
                        content = content.decode('ascii')
                    with w.tag(u'pre'):
                        w._flush()
                        for line in content.splitlines():
                            write_votlint_warning(w, line, xml_lines)


def write_result_row(w, result):
    with w.tag(u'tr'):
        with w.tag(u'td'):
            if ('network_error' in result and
                    result['network_error'] is not None):
                w.data(result.url.decode('ascii'))
            else:
                w.element(u'a', result.url.decode('ascii'),
                          href=u'%s/index.html' % result.get_htmlpath())

        if 'network_error' in result and result['network_error'] is not None:
            w.element(u'td', unicode(result['network_error']),
                      attrib={u'class': u'red'})
            w.element(u'td', '-')
            w.element(u'td', '-')
            w.element(u'td', '-')
            w.element(u'td', '-')
        else:
            w.element(u'td', u'-', attrib={u'class': u'green'})

            if result['nexceptions']:
                cls = u'red'
                msg = u'Fatal'
            elif result['nwarnings']:
                cls = u'yellow'
                msg = unicode(result['nwarnings'])
            else:
                cls = u'green'
                msg = u'-'
            w.element(u'td', msg, attrib={u'class': cls})

            msg = result['version']
            if result['xmllint'] is None:
                cls = u''
            elif result['xmllint'] is False:
                cls = u'red'
            else:
                cls = u'green'
            w.element(u'td', msg, attrib={u'class': cls})

            if result['expected'] == 'good':
                cls = u'green'
                msg = u'-'
            elif result['expected'] == 'broken':
                cls = u'red'
                msg = u'net'
            elif result['expected'] == 'incorrect':
                cls = u'yellow'
                msg = u'invalid'
            w.element(u'td', msg, attrib={u'class': cls})

            if 'votlint' in result:
                if result['votlint']:
                    cls = u'green'
                    msg = u'Passed'
                else:
                    cls = u'red'
                    msg = u'Failed'
            else:
                cls = u''
                msg = '?'
            w.element(u'td', msg, attrib={u'class': cls})


def write_table(basename, name, results, root="results", chunk_size=500):
    def write_page_links(j):
        if npages <= 1:
            return
        with w.tag(u'center'):
            if j > 0:
                w.element(u'a', u'<< ', href=u'%s_%02d.html' % (basename, j-1))
            for i in xrange(npages):
                if i == j:
                    w.data(unicode(i+1))
                else:
                    w.element(
                        u'a', unicode(i+1),
                        href=u'%s_%02d.html' % (basename, i))
                w.data(' ')
            if j < npages - 1:
                w.element(u'a', u'>>', href=u'%s_%02d.html' % (basename, j+1))

    npages = int(ceil(float(len(results)) / chunk_size))

    for i, j in enumerate(xrange(0, max(len(results), 1), chunk_size)):
        subresults = results[j:j+chunk_size]
        path = os.path.join(root, '%s_%02d.html' % (basename, i))
        with io.open(path, 'w', encoding='utf-8') as fd:
            w = XMLWriter(fd)
            with make_html_header(w):
                write_page_links(i)

                w.element('h2', name)

                with w.tag(u'table'):
                    with w.tag(u'tr'):
                        w.element(u'th', u'URL')
                        w.element(u'th', u'Network')
                        w.element(u'th', u'Warnings')
                        w.element(u'th', u'Schema')
                        w.element(u'th', u'Expected')
                        w.element(u'th', u'votlint')

                    for result in subresults:
                        write_result_row(w, result)

                write_page_links(i)


def add_subset(w, basename, name, subresults, inside=['p'], total=None):
    with w.tag('tr'):
        subresults = list(subresults)
        if total is None:
            total = len(subresults)
        percentage = (float(len(subresults)) / total)
        with w.tag('td'):
            for element in inside:
                w.start(element)
            w.element(u'a', name, href=u'%s_00.html' % basename)
            for element in reversed(inside):
                w.end(element)
        numbers = '%d (%.2f%%)' % (len(subresults), percentage * 100.0)
        with w.tag('td'):
            w.data(numbers)


def write_index(subsets, results, root='results'):
    path = os.path.join(root, 'index.html')
    with io.open(path, 'w', encoding='utf-8') as fd:
        w = XMLWriter(fd)
        with make_html_header(w):
            w.element(u'h1', u'VO Validation results')

            with w.tag('table'):
                for subset in subsets:
                    add_subset(w, *subset, total=len(results))


def write_index_table(root, basename, name, subresults, inside=None,
                      total=None, chunk_size=500):
    if total is None:
        total = len(subresults)
    percentage = (float(len(subresults)) / total)
    numbers = '%d (%.2f%%)' % (len(subresults), percentage * 100.0)
    write_table(basename, name + ' ' + numbers, subresults, root, chunk_size)
