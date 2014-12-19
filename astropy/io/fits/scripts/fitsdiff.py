# Licensed under a 3-clause BSD style license - see LICENSE.rst
import glob
import logging
import optparse
import os
import sys
import textwrap

from ... import fits
from ..util import fill


log = logging.getLogger('fitsdiff')


USAGE = """
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

This script is part of the Astropy package. 
More documentation can be found here:
http://astropy.readthedocs.org/en/latest/io/fits/usage/scripts.html
""".strip()


EPILOG = """
If the two files are identical within the specified conditions, it will report
"No difference is found." If the value(s) of -c and -k takes the form
'@filename', list is in the text file 'filename', and each line in that text
file contains one keyword.

Example
-------

    fitsdiff -k filename,filtnam1 -n 5 -d 1.e-6 test1.fits test2

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
""".strip()


class HelpFormatter(optparse.TitledHelpFormatter):
    def format_epilog(self, epilog):
        return '\n%s\n' % fill(epilog, self.width)


def handle_options(argv=None):
    # This is a callback--less trouble than actually adding a new action type
    def store_list(option, opt, value, parser):
        setattr(parser.values, option.dest, [])
        # Accept either a comma-separated list or a filename (starting with @)
        # containing a value on each line
        if value and value[0] == '@':
            value = value[1:]
            if not os.path.exists(value):
                log.warning('%s argument %s does not exist' % (opt, value))
                return
            try:
                values = [v.strip() for v in open(value, 'r').readlines()]
                setattr(parser.values, option.dest, values)
            except IOError as exc:
                log.warning('reading %s for %s failed: %s; ignoring this '
                            'argument' % (value, opt, exc))
                del exc
        else:
            setattr(parser.values, option.dest,
                    [v.strip() for v in value.split(',')])

    parser = optparse.OptionParser(usage=USAGE, epilog=EPILOG,
                                   formatter=HelpFormatter())

    parser.add_option(
        '-q', '--quiet', action='store_true',
        help='Produce no output and just return a status code.')

    parser.add_option(
        '-n', '--num-diffs', type='int', default=10, dest='numdiffs',
        metavar='INTEGER',
        help='Max number of data differences (image pixel or table element) '
             'to report per extension (default %default).')

    parser.add_option(
        '-d', '--difference-tolerance', type='float', default=0.,
        dest='tolerance', metavar='NUMBER',
        help='The relative tolerance for comparison of two numbers, '
             'specifically two floating point numbers.  This applies to data '
             'in both images and tables, and to floating point keyword values '
             'in headers (default %default).')

    parser.add_option(
        '-b', '--no-ignore-blanks', action='store_false',
        dest='ignore_blanks', default=True,
        help="Don't ignore trailing blanks (whitespace) in string values.  "
             "Otherwise trailing blanks both in header keywords/values and in "
             "table column values) are not treated as significant i.e., "
             "without this option 'ABCDEF   ' and 'ABCDEF' are considered "
             "equivalent. ")

    parser.add_option(
        '--no-ignore-blank-cards', action='store_false',
        dest='ignore_blank_cards', default=True,
        help="Don't ignore entirely blank cards in headers.  Normally fitsdiff "
             "does not consider blank cards when comparing headers, but this "
             "will ensure that even blank cards match up. ")

    parser.add_option(
        '--exact', action='store_true',
        dest='exact_comparisons', default=False,
        help="Report ALL differences, "
             "overriding command-line options and FITSDIFF_SETTINGS. ")

    parser.add_option(
        '-o', '--output-file', metavar='FILE',
        help='Output results to this file; otherwise results are printed to '
             'stdout.')

    group = optparse.OptionGroup(parser, 'Header Comparison Options')

    group.add_option(
        '-k', '--ignore-keywords', action='callback', callback=store_list,
        nargs=1, type='str', default=[], dest='ignore_keywords',
        metavar='KEYWORDS',
        help='Comma-separated list of keywords not to be compared.  Keywords '
             'may contain wildcard patterns.  To exclude all keywords, use '
             '"*"; make sure to have double or single quotes around the '
             'asterisk on the command-line.')

    group.add_option(
        '-c', '--ignore-comments', action='callback', callback=store_list,
        nargs=1, type='str', default=[], dest='ignore_comments',
        metavar='KEYWORDS',
        help='Comma-separated list of keywords whose comments will not be '
             'compared.  Wildcards may be used as with --ignore-keywords.')

    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, 'Table Comparison Options')

    group.add_option(
        '-f', '--ignore-fields', action='callback', callback=store_list,
        nargs=1, type='str', default=[], dest='ignore_fields',
        metavar='COLUMNS',
        help='Comma-separated list of fields (i.e. columns) not to be '
             'compared.  All columns may be excluded using "*" as with '
             '--ignore-keywords.')

    parser.add_option_group(group)
    options, args = parser.parse_args(argv)

    # Determine which filenames to compare
    if len(args) != 2:
        parser.error('\n' + textwrap.fill(
            'fitsdiff requires two arguments; see `fitsdiff --help` for more '
            'details.', parser.formatter.width))

    return options, args


def setup_logging(outfile=None):
    log.setLevel(logging.INFO)
    error_handler = logging.StreamHandler(sys.stderr)
    error_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    error_handler.setLevel(logging.WARNING)
    log.addHandler(error_handler)

    if outfile is not None:
        output_handler = logging.FileHandler(outfile)
    else:
        output_handler = logging.StreamHandler()

        class LevelFilter(logging.Filter):
            """Log only messages matching the specified level."""

            def __init__(self, name='', level=logging.NOTSET):
                logging.Filter.__init__(self, name)
                self.level = level

            def filter(self, rec):
                return rec.levelno == self.level

        # File output logs all messages, but stdout logs only INFO messages
        # (since errors are already logged to stderr)
        output_handler.addFilter(LevelFilter(level=logging.INFO))

    output_handler.setFormatter(logging.Formatter('%(message)s'))
    log.addHandler(output_handler)


def match_files(paths):
    filelists = []

    for path in paths:
        if glob.has_magic(path):
            files = [os.path.abspath(f) for f in glob.glob(path)]
            if not files:
                log.error(
                    'Wildcard pattern %r did not match any files.' % path)
                sys.exit(2)
            filelists.append(files)
        elif os.path.isdir(path):
            filelists.append([os.path.abspath(f) for f in os.listdir(path)])
        elif os.path.isfile(path):
            filelists.append([path])
        else:
            log.error(
                '%r is not an existing file, directory, or wildcard pattern; '
                'see `fitsdiff --help` for more usage help.' % path)
            sys.exit(2)

    filelists[0].sort()
    filelists[1].sort()

    for a, b in [(0, 1), (1, 0)]:
        if len(filelists[a]) > len(filelists[b]):
            for extra in filelists[a][len(filelists[b]):]:
                log.warning('%r has no match in %r' % (extra, paths[b]))
            filelists[a] = filelists[a][:len(filelists[b])]
            break

    return zip(*filelists)


def main():
    if 'FITSDIFF_SETTINGS' in os.environ:
        argv = os.environ['FITSDIFF_SETTINGS'].split() + sys.argv[1:]
    else:
        argv = sys.argv[1:]

    opts, args = handle_options(argv)

    if opts.exact_comparisons:
        # override the options so that each is the most restrictive
        opts.ignore_keywords = []
        opts.ignore_comments = []
        opts.ignore_fields = []
        opts.tolerance = 0.0
        opts.ignore_blanks = False
        opts.ignore_blank_cards = False

    if not opts.quiet:
        setup_logging(opts.output_file)
    files = match_files(args)

    close_file = False
    if opts.quiet:
        out_file = None
    elif opts.output_file:
        out_file = open(opts.output_file, 'wb')
        close_file = True
    else:
        out_file = sys.stdout

    identical = []
    try:
        for a, b in files:
            # TODO: pass in any additonal arguments here too
            diff = fits.diff.FITSDiff(
                a, b,
                ignore_keywords=opts.ignore_keywords,
                ignore_comments=opts.ignore_comments,
                ignore_fields=opts.ignore_fields,
                numdiffs=opts.numdiffs,
                tolerance=opts.tolerance,
                ignore_blanks=opts.ignore_blanks,
                ignore_blank_cards=opts.ignore_blank_cards)
            diff.report(fileobj=out_file)
            identical.append(diff.identical)

        return int(not all(identical))
    finally:
        if close_file:
            out_file.close()
