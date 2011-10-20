.. _verify:

********************
Verification options
********************

There are 5 options for the `output_verify` argument of the following methods:
:meth:`close`, :meth:`writeto`, and :meth:`flush`. In these cases, they are
passed to a :meth:`verify` call within these methods.

exception
=========

This option will raise an exception if any FITS standard is violated. This is
the default option for output (i.e. when :meth:`writeto`, :meth:`close`, or
:meth:`flush` is called. If a user wants to overwrite this default on output,
the other options listed below can be used.

ignore
======

This option will ignore any FITS standard violation. On output, it will write
the HDU List content to the output FITS file, whether or not it is conforming
to FITS standard.

The `ignore` option is useful in these situations, for example:

  1. An input FITS file with non-standard is read and the user wants to copy or
     write out after some modification to an output file. The non-standard will
     be preserved in such output file.

  2. A user wants to create a non-standard FITS file on purpose, possibly for
     testing purpose.

No warning message will be printed out. This is like a silent warn (see below)
option.

fix
===

This option wil try to fix any FITS standard violations. It is not always
possible to fix such violations. In general, there are two kinds of FITS
standard violation: fixable and not fixable. For example, if a keyword has a
floating number with an exponential notation in lower case ’e’ (e.g. 1.23e11)
instead of the upper case ’E’ as required by the FITS standard, it's a fixable
violation. On the other hand, a keyword name like ``P.I.`` is not fixable,
since it will not know what to use to replace the disallowed periods. If a
violation is fixable, this option will print out a message noting it is fixed.
If it is not fixable, it will throw an exception.

The principle behind the fixing is do no harm. For example, it is plausible to
’fix’ a `Card` with a keyword name like ``P.I.`` by deleting it, but AstroPy
will not take such action to hurt the integrity of the data.

Not all fixes may be the "correct" fix, but at least AstroPy will try to make
the fix in such a way that it will not throw off other FITS readers.

silentfix
=========

Same as fix, but will not print out informative messages. This may be useful in
a large script where the user does not want excessive harmless messages. If the
violation is not fixable, it will still throw an exception.

warn
====

This option is the same as the ignore option but will send warning messages. It
will not try to fix any FITS standard violations whether fixable or not.
