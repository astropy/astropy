# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
import contextlib
import os
import re
from math import ceil

from astropy import online_docs_root
from astropy.io.votable import exceptions
from astropy.utils.xml.writer import XMLWriter, xml_escape

html_header = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html
        PUBLIC "-//W3C//DTD XHTML Basic 1.0//EN"
        "http://www.w3.org/TR/xhtml-basic/xhtml-basic10.dtd">
"""

default_style = """
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
    with w.tag("html", xmlns="http://www.w3.org/1999/xhtml", lang="en-US"):
        with w.tag("head"):
            w.element("title", "VO Validation results")
            w.element("style", default_style)

            with w.tag("body"):
                yield


def write_source_line(w, line, nchar=0):
    part1 = xml_escape(line[:nchar].decode("utf-8"))
    char = xml_escape(line[nchar : nchar + 1].decode("utf-8"))
    part2 = xml_escape(line[nchar + 1 :].decode("utf-8"))

    w.write("  ")
    w.write(part1)
    w.write(f'<span class="highlight">{char}</span>')
    w.write(part2)
    w.write("\n\n")


def write_warning(w, line, xml_lines):
    warning = exceptions.parse_vowarning(line)
    if not warning["is_something"]:
        w.data(line)
    else:
        w.write(f"Line {warning['nline']:d}: ")
        if warning["warning"]:
            w.write(
                '<a href="{}/{}">{}</a>: '.format(
                    online_docs_root, warning["doc_url"], warning["warning"]
                )
            )
        msg = warning["message"]
        if not isinstance(warning["message"], str):
            msg = msg.decode("utf-8")
        w.write(xml_escape(msg))
        w.write("\n")
        if 1 <= warning["nline"] < len(xml_lines):
            write_source_line(w, xml_lines[warning["nline"] - 1], warning["nchar"])


def write_votlint_warning(w, line, xml_lines):
    match = re.search(
        r"(WARNING|ERROR|INFO) \(l.(?P<line>[0-9]+), c.(?P<column>[0-9]+)\):"
        r" (?P<rest>.*)",
        line,
    )
    if match:
        w.write(f"Line {int(match.group('line'))}: {xml_escape(match.group('rest'))}\n")
        write_source_line(
            w, xml_lines[int(match.group("line")) - 1], int(match.group("column")) - 1
        )
    else:
        w.data(line)
        w.data("\n")


def write_result(result):
    if "network_error" in result and result["network_error"] is not None:
        return

    xml = result.get_xml_content()
    xml_lines = xml.splitlines()

    path = os.path.join(result.get_dirpath(), "index.html")

    with open(path, "w", encoding="utf-8") as fd:
        w = XMLWriter(fd)
        with make_html_header(w):
            with w.tag("p"):
                with w.tag("a", href="vo.xml"):
                    w.data(result.url.decode("ascii"))
            w.element("hr")

            with w.tag("pre"):
                w._flush()
                for line in result["warnings"]:
                    write_warning(w, line, xml_lines)

            if result["xmllint"] is False:
                w.element("hr")
                w.element("p", "xmllint results:")
                content = result["xmllint_content"]
                if not isinstance(content, str):
                    content = content.decode("ascii")
                content = content.replace(result.get_dirpath() + "/", "")
                with w.tag("pre"):
                    w.data(content)

            if "votlint" in result:
                if result["votlint"] is False:
                    w.element("hr")
                    w.element("p", "votlint results:")
                    content = result["votlint_content"]
                    if not isinstance(content, str):
                        content = content.decode("ascii")
                    with w.tag("pre"):
                        w._flush()
                        for line in content.splitlines():
                            write_votlint_warning(w, line, xml_lines)


def write_result_row(w, result):
    with w.tag("tr"):
        with w.tag("td"):
            if "network_error" in result and result["network_error"] is not None:
                w.data(result.url.decode("ascii"))
            else:
                w.element(
                    "a",
                    result.url.decode("ascii"),
                    href=f"{result.get_htmlpath()}/index.html",
                )

        if "network_error" in result and result["network_error"] is not None:
            w.element("td", str(result["network_error"]), attrib={"class": "red"})
            w.element("td", "-")
            w.element("td", "-")
            w.element("td", "-")
            w.element("td", "-")
        else:
            w.element("td", "-", attrib={"class": "green"})

            if result["nexceptions"]:
                cls = "red"
                msg = "Fatal"
            elif result["nwarnings"]:
                cls = "yellow"
                msg = str(result["nwarnings"])
            else:
                cls = "green"
                msg = "-"
            w.element("td", msg, attrib={"class": cls})

            msg = result["version"]
            if result["xmllint"] is None:
                cls = ""
            elif result["xmllint"] is False:
                cls = "red"
            else:
                cls = "green"
            w.element("td", msg, attrib={"class": cls})

            if result["expected"] == "good":
                cls = "green"
                msg = "-"
            elif result["expected"] == "broken":
                cls = "red"
                msg = "net"
            elif result["expected"] == "incorrect":
                cls = "yellow"
                msg = "invalid"
            w.element("td", msg, attrib={"class": cls})

            if "votlint" in result:
                if result["votlint"]:
                    cls = "green"
                    msg = "Passed"
                else:
                    cls = "red"
                    msg = "Failed"
            else:
                cls = ""
                msg = "?"
            w.element("td", msg, attrib={"class": cls})


def write_table(basename, name, results, root="results", chunk_size=500):
    def write_page_links(j):
        if npages <= 1:
            return
        with w.tag("center"):
            if j > 0:
                w.element("a", "<< ", href=f"{basename}_{j - 1:02d}.html")
            for i in range(npages):
                if i == j:
                    w.data(str(i + 1))
                else:
                    w.element("a", str(i + 1), href=f"{basename}_{i:02d}.html")
                w.data(" ")
            if j < npages - 1:
                w.element("a", ">>", href=f"{basename}_{j + 1:02d}.html")

    npages = ceil(float(len(results)) / chunk_size)

    for i, j in enumerate(range(0, max(len(results), 1), chunk_size)):
        subresults = results[j : j + chunk_size]
        path = os.path.join(root, f"{basename}_{i:02d}.html")
        with open(path, "w", encoding="utf-8") as fd:
            w = XMLWriter(fd)
            with make_html_header(w):
                write_page_links(i)

                w.element("h2", name)

                with w.tag("table"):
                    with w.tag("tr"):
                        w.element("th", "URL")
                        w.element("th", "Network")
                        w.element("th", "Warnings")
                        w.element("th", "Schema")
                        w.element("th", "Expected")
                        w.element("th", "votlint")

                    for result in subresults:
                        write_result_row(w, result)

                write_page_links(i)


def add_subset(w, basename, name, subresults, inside=["p"], total=None):
    with w.tag("tr"):
        subresults = list(subresults)
        if total is None:
            total = len(subresults)
        if total == 0:  # pragma: no cover
            percentage = 0.0
        else:
            percentage = float(len(subresults)) / total
        with w.tag("td"):
            for element in inside:
                w.start(element)
            w.element("a", name, href=f"{basename}_00.html")
            for element in reversed(inside):
                w.end(element)
        numbers = f"{len(subresults):d} ({percentage:.2%})"
        with w.tag("td"):
            w.data(numbers)


def write_index(subsets, results, root="results"):
    path = os.path.join(root, "index.html")
    with open(path, "w", encoding="utf-8") as fd:
        w = XMLWriter(fd)
        with make_html_header(w):
            w.element("h1", "VO Validation results")

            with w.tag("table"):
                for subset in subsets:
                    add_subset(w, *subset, total=len(results))


def write_index_table(
    root, basename, name, subresults, inside=None, total=None, chunk_size=500
):
    if total is None:
        total = len(subresults)
    percentage = float(len(subresults)) / total
    numbers = f"{len(subresults):d} ({percentage:.2%})"
    write_table(basename, name + " " + numbers, subresults, root, chunk_size)
