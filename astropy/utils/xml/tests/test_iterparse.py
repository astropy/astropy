# Licensed under a 3-clause BSD style license - see LICENSE.rst

# LOCAL
from ....utils.xml.iterparser import _fast_iterparse

# SYSTEM
import io
import zlib

# The C-based XML parser for VOTables previously used fixed-sized
# buffers (allocated at __init__() time).  This test will
# only pass with the patch that allows a dynamic realloc() of
# the queue.  This addresses the bugs:
#
# - "RuntimeError: XML queue overflow"
#   https://github.com/astropy/astropy/issues/5824
#   (Kudos to Stefan Becker---ARI/ZAH Heidelberg)
#
# - "iterparse.c: add queue_realloc() + move 'buffersize / 2' logic there"
#   https://github.com/astropy/astropy/issues/5869
#
# This test code can emulate a combination of network buffering and
# gzip decompression---with different request sizes, it can be used to
# demonstrate both under-reading and over-reading.
#
# Using the 512-tag VOTABLE XML sample input, and various combinations
# of minimum/maximum fetch sizes, the following situations can be
# generated:
#
# maximum_fetch =  1 (ValueError, no element found) still within gzip headers
# maximum_fetch = 80 (ValueError, unclosed token) short read
# maximum_fetch =217 passes, because decompressed_length > requested
#                            && <512 tags in a single parse
# maximum_fetch =218 (RuntimeError, XML queue overflow)
#
# The test provided here covers the over-reading identified in #5824
# (equivalent to the 217).

# Firstly, assemble a minimal VOTABLE header, table contents and footer.
# This is done in textual form, as the aim is to only test the parser, not
# the outputter!
HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE>
 <RESOURCE type="results">
  <TABLE>
   <FIELD ID="foo" name="foo" datatype="int" arraysize="1"/>
    <DATA>
     <TABLEDATA>
"""

ROW = """<TR><TD>0</TD></TR>
"""

FOOTER = """
    </TABLEDATA>
   </DATA>
  </TABLE>
 </RESOURCE>
</VOTABLE>
"""

# minimum passable buffer size => 1024
# 1024 / 2 => 512 tags for overflow
# 512 - 7 tags in header, - 5 tags in footer = 500 tags required for overflow
# 500 / 4 tags (<tr><td></td></tr>) per row == 125 rows required for overflow
VOTABLE_XML = HEADER + 125*ROW + FOOTER

# UngzipFileWrapper() wraps an existing file-like Object,
# decompressing the content and returning the plaintext.
# This therefore emulates the behavior of the Python 'requests'
# library when transparently decompressing Gzip HTTP responses.
#
# The critical behavior is that---because of the
# decompression---read() can return considerably more
# bytes than were requested!  (But, read() can also return less).
#
# inspiration:
# http://stackoverflow.com/questions/4013843/how-to-wrap-file-object-read-and-write-operation-which-are-readonly


class UngzipFileWrapper:
    def __init__(self, fd, **kwargs):
        self._file = fd
        self._z = zlib.decompressobj(16 + zlib.MAX_WBITS)

    def read(self, requested_length):
        # emulate network buffering dynamics by clamping the read size
        clamped_length = max(1, min(1 << 24, requested_length))
        compressed = self._file.read(clamped_length)
        plaintext = self._z.decompress(compressed)
        # Only for real local files---just for the testcase
        if len(compressed) == 0:
            self.close()
        return plaintext

    def __getattr__(self, attr):
        return getattr(self._file, attr)

# test_iterparser_over_read_simple() is a very cut down test,
# of the original more flexible test-case, but without external
# dependencies.  The plaintext is compressed and then decompressed
# to provide a better emulation of the original situation where
# the bug was observed.
#
# If a dependency upon 'zlib' is not desired, it would be possible to
# simplify this testcase by replacing the compress/decompress with a
# read() method emulation that always returned more from a buffer tha
# was requested.


def test_iterparser_over_read_simple():
    # Take the plaintext of 512 tags, and compression it with a
    # Gzip-style header (+16), to most closely emulate the behavior
    # of most HTTP servers.
    zlib_GZIP_STYLE_HEADER = 16
    compo = zlib.compressobj(zlib.Z_BEST_COMPRESSION,
                             zlib.DEFLATED,
                             zlib.MAX_WBITS + zlib_GZIP_STYLE_HEADER)

    # Bytes vs. String  .encode()/.decode() for compatibility with Python 3.5.
    s = compo.compress(VOTABLE_XML.encode())
    s = s + compo.flush()
    fd = io.BytesIO(s)
    fd.seek(0)

    # Finally setup the test of the C-based '_fast_iterparse()' iterator
    # and a situation in which it can be called a-la the VOTable Parser.
    MINIMUM_REQUESTABLE_BUFFER_SIZE = 1024
    uncompressed_fd = UngzipFileWrapper(fd)
    iterable = _fast_iterparse(uncompressed_fd.read,
                               MINIMUM_REQUESTABLE_BUFFER_SIZE)
    list(iterable)
