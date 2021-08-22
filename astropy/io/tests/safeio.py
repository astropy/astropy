# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io


class CatchZeroByteWriter(io.BufferedWriter):
    """File handle to intercept 0-byte writes"""

    def write(self, buffer):
        nbytes = super().write(buffer)
        if nbytes == 0:
            raise ValueError("This writer does not allow empty writes")
        return nbytes
