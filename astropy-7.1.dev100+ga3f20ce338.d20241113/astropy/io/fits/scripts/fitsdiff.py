# Licensed under a 3-clause BSD style license - see LICENSE.rst
import argparse
import glob
import logging
import os
import sys

from astropy import __version__
from astropy.io import fits
from astropy.io.fits.util import fill

log = logging.getLogger("fitsdiff")


DESCRIPTION = """
Compare two FITS image files and report the differences in header keywords and
data.

    fitsdiff [options] filename1 filename2

where filename1 filename2 are the two files to be compared.  They may also be
wild cards, in such cases, they must be enclosed by double or single quotes, or
they may be directory names.  If both are directory names, all files in each of
the directories will be included; if only one is a directory name, then the
directory name will be prefixed to the file name(s) specified by the other
argument.  for example::

    fitsdiff "*.fits" "/machine/data1"

will compare all FITS files in the current directory to the corresponding files
in the directory /machine/data1.

This script is part of the Astropy package. See
https://docs.astropy.org/en/latest/io/fits/usage/scripts.html#fitsdiff
for further documentation.
""".strip()


EPILOG = fill(
    """
If the two files are identical within the specified conditions, it will report
"No difference is found." If the value(s) of -c and -k takes the form
'@filename', list is in the text file 'filename', and each line in that text
file contains one keyword.

Example
-------

    fitsdiff -k filename,filtnam1 -n 5 -r 1.e-6 test1.fits test2

This command will compare files test1.fits and test2.fits, report maximum of 5
different pixels values per extension, only report data values larger than
1.e-6 relative to each other, and will neglect the different values of keywords
FILENAME and FILTNAM1 (or their very existence).

fitsdiff command-line arguments can also be set using the environment variable
FITSDIFF_SETTINGS.  If the FITSDIFF_SETTINGS environment variable is present,
each argument present will override the corresponding argument on the
command-line unless the --exact option is specified.  The FITSDIFF_SETTINGS
environment variable exists to make it easier to change the
behavior of fitsdiff on a global level, such as in a set of regression tests.
""".strip(),
    width=80,
)


class StoreListAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super().__init__(option_strings, dest, nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, [])
        # Accept either a comma-separated list or a filename (starting with @)
        # containing a value on each line
        if values and values[0] == "@":
            value = values[1:]
            if not os.path.exists(value):
                log.warning(f"{self.dest} argument {value} does not exist")
                return
            try:
                values = [v.strip() for v in open(value)]
                setattr(namespace, self.dest, values)
            except OSError as exc:
                log.warning(
                    f"reading {value} for {self.dest} failed: {exc}; "
                    "ignoring this argument"
                )
                del exc
        else:
            setattr(namespace, self.dest, [v.strip() for v in values.split(",")])


def handle_options(argv=None):
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "fits_files", metavar="file", nargs="+", help=".fits files to process."
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Produce no output and just return a status code.",
    )

    parser.add_argument(
        "-n",
        "--num-diffs",
        type=int,
        default=10,
        dest="numdiffs",
        metavar="INTEGER",
        help=(
            "Max number of data differences (image pixel or table element) "
            "to report per extension (default %(default)s)."
        ),
    )

    parser.add_argument(
        "-r",
        "--rtol",
        "--relative-tolerance",
        type=float,
        default=None,
        dest="rtol",
        metavar="NUMBER",
        help=(
            "The relative tolerance for comparison of two numbers, "
            "specifically two floating point numbers.  This applies to data "
            "in both images and tables, and to floating point keyword values "
            "in headers (default %(default)s)."
        ),
    )

    parser.add_argument(
        "-a",
        "--atol",
        "--absolute-tolerance",
        type=float,
        default=None,
        dest="atol",
        metavar="NUMBER",
        help=(
            "The absolute tolerance for comparison of two numbers, "
            "specifically two floating point numbers.  This applies to data "
            "in both images and tables, and to floating point keyword values "
            "in headers (default %(default)s)."
        ),
    )

    parser.add_argument(
        "-b",
        "--no-ignore-blanks",
        action="store_false",
        dest="ignore_blanks",
        default=True,
        help=(
            "Don't ignore trailing blanks (whitespace) in string values.  "
            "Otherwise trailing blanks both in header keywords/values and in "
            "table column values) are not treated as significant i.e., "
            "without this option 'ABCDEF   ' and 'ABCDEF' are considered "
            "equivalent. "
        ),
    )

    parser.add_argument(
        "--no-ignore-blank-cards",
        action="store_false",
        dest="ignore_blank_cards",
        default=True,
        help=(
            "Don't ignore entirely blank cards in headers.  Normally fitsdiff "
            "does not consider blank cards when comparing headers, but this "
            "will ensure that even blank cards match up. "
        ),
    )

    parser.add_argument(
        "--exact",
        action="store_true",
        dest="exact_comparisons",
        default=False,
        help=(
            "Report ALL differences, "
            "overriding command-line options and FITSDIFF_SETTINGS. "
        ),
    )

    parser.add_argument(
        "-o",
        "--output-file",
        metavar="FILE",
        help="Output results to this file; otherwise results are printed to stdout.",
    )

    parser.add_argument(
        "-u",
        "--ignore-hdus",
        action=StoreListAction,
        default=[],
        dest="ignore_hdus",
        metavar="HDU_NAMES",
        help=(
            "Comma-separated list of HDU names not to be compared.  HDU "
            "names may contain wildcard patterns."
        ),
    )

    group = parser.add_argument_group("Header Comparison Options")

    group.add_argument(
        "-k",
        "--ignore-keywords",
        action=StoreListAction,
        default=[],
        dest="ignore_keywords",
        metavar="KEYWORDS",
        help=(
            "Comma-separated list of keywords not to be compared.  Keywords "
            "may contain wildcard patterns.  To exclude all keywords, use "
            '"*"; make sure to have double or single quotes around the '
            "asterisk on the command-line."
        ),
    )

    group.add_argument(
        "-c",
        "--ignore-comments",
        action=StoreListAction,
        default=[],
        dest="ignore_comments",
        metavar="COMMENTS",
        help=(
            "Comma-separated list of keywords whose comments will not be "
            "compared.  Wildcards may be used as with --ignore-keywords."
        ),
    )

    group = parser.add_argument_group("Table Comparison Options")

    group.add_argument(
        "-f",
        "--ignore-fields",
        action=StoreListAction,
        default=[],
        dest="ignore_fields",
        metavar="COLUMNS",
        help=(
            "Comma-separated list of fields (i.e. columns) not to be "
            'compared.  All columns may be excluded using "*" as with '
            "--ignore-keywords."
        ),
    )

    options = parser.parse_args(argv)

    # Determine which filenames to compare
    if len(options.fits_files) != 2:
        parser.error(
            "\nfitsdiff requires two arguments; see `fitsdiff --help` for more details."
        )

    return options


def setup_logging(outfile=None):
    log.setLevel(logging.INFO)
    error_handler = logging.StreamHandler(sys.stderr)
    error_handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    error_handler.setLevel(logging.WARNING)
    log.addHandler(error_handler)

    if outfile is not None:
        output_handler = logging.FileHandler(outfile)
    else:
        output_handler = logging.StreamHandler()

        class LevelFilter(logging.Filter):
            """Log only messages matching the specified level."""

            def __init__(self, name="", level=logging.NOTSET):
                logging.Filter.__init__(self, name)
                self.level = level

            def filter(self, rec):
                return rec.levelno == self.level

        # File output logs all messages, but stdout logs only INFO messages
        # (since errors are already logged to stderr)
        output_handler.addFilter(LevelFilter(level=logging.INFO))

    output_handler.setFormatter(logging.Formatter("%(message)s"))
    log.addHandler(output_handler)


def match_files(paths):
    if os.path.isfile(paths[0]) and os.path.isfile(paths[1]):
        # shortcut if both paths are files
        return [paths]

    dirnames = [None, None]
    filelists = [None, None]

    for i, path in enumerate(paths):
        if glob.has_magic(path):
            files = [os.path.split(f) for f in glob.glob(path)]
            if not files:
                log.error("Wildcard pattern %r did not match any files.", path)
                sys.exit(2)

            dirs, files = list(zip(*files))
            if len(set(dirs)) > 1:
                log.error("Wildcard pattern %r should match only one directory.", path)
                sys.exit(2)

            dirnames[i] = set(dirs).pop()
            filelists[i] = sorted(files)
        elif os.path.isdir(path):
            dirnames[i] = path
            filelists[i] = [
                f
                for f in sorted(os.listdir(path))
                if os.path.isfile(os.path.join(path, f))
            ]
        elif os.path.isfile(path):
            dirnames[i] = os.path.dirname(path)
            filelists[i] = [os.path.basename(path)]
        else:
            log.error(
                "%r is not an existing file, directory, or wildcard "
                "pattern; see `fitsdiff --help` for more usage help.",
                path,
            )
            sys.exit(2)

        dirnames[i] = os.path.abspath(dirnames[i])

    filematch = set(filelists[0]) & set(filelists[1])

    for a, b in [(0, 1), (1, 0)]:
        if len(filelists[a]) > len(filematch) and not os.path.isdir(paths[a]):
            for extra in sorted(set(filelists[a]) - filematch):
                log.warning("%r has no match in %r", extra, dirnames[b])

    return [
        (os.path.join(dirnames[0], f), os.path.join(dirnames[1], f)) for f in filematch
    ]


def main(args=None):
    args = args or sys.argv[1:]

    if "FITSDIFF_SETTINGS" in os.environ:
        args = os.environ["FITSDIFF_SETTINGS"].split() + args

    opts = handle_options(args)

    if opts.rtol is None:
        opts.rtol = 0.0
    if opts.atol is None:
        opts.atol = 0.0

    if opts.exact_comparisons:
        # override the options so that each is the most restrictive
        opts.ignore_keywords = []
        opts.ignore_comments = []
        opts.ignore_fields = []
        opts.rtol = 0.0
        opts.atol = 0.0
        opts.ignore_blanks = False
        opts.ignore_blank_cards = False

    if not opts.quiet:
        setup_logging(opts.output_file)
    files = match_files(opts.fits_files)

    close_file = False
    if opts.quiet:
        out_file = None
    elif opts.output_file:
        out_file = open(opts.output_file, "w")
        close_file = True
    else:
        out_file = sys.stdout

    identical = []
    try:
        for a, b in files:
            # TODO: pass in any additional arguments here too
            diff = fits.diff.FITSDiff(
                a,
                b,
                ignore_hdus=opts.ignore_hdus,
                ignore_keywords=opts.ignore_keywords,
                ignore_comments=opts.ignore_comments,
                ignore_fields=opts.ignore_fields,
                numdiffs=opts.numdiffs,
                rtol=opts.rtol,
                atol=opts.atol,
                ignore_blanks=opts.ignore_blanks,
                ignore_blank_cards=opts.ignore_blank_cards,
            )

            diff.report(fileobj=out_file)
            identical.append(diff.identical)

        return int(not all(identical))
    finally:
        if close_file:
            out_file.close()
        # Close the file if used for the logging output, and remove handlers to
        # avoid having them multiple times for unit tests.
        for handler in log.handlers:
            if isinstance(handler, logging.FileHandler):
                handler.close()
            log.removeHandler(handler)
