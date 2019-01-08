.. _bitmask_details:

*********************************************************
Utility functions for handling bit masks and mask arrays.
*********************************************************

It is common to use `bit fields <https://en.wikipedia.org/wiki/Bit_field>`_ - \
e.g., integer variables whose individual bits
represent some attributes - to characterize the state of data. For example,
HST uses arrays of bit fields to characterize data quality (DQ) of HST images,
see, e.g., DQ field values for
`WFPC2 image data <http://documents.stsci.edu/hst/wfpc2/documents/handbooks/dhb/wfpc2_ch34.html#1971480>`_
and `WFC3 image data <http://www.stsci.edu/hst/wfc3/documents/handbooks/currentDHB/Chapter2_data_structure3.html#567105>`_.
As one can see, the meaning assigned to various *bit flags* in for the two
instruments is generally different.

Bit fields can be thought of as tightly packed collections of bit flags. Using
`masking <https://en.wikipedia.org/wiki/Mask_(computing)>`_ we can "inspect"
the status of individual bits.

One common operation performed on bit field arrays is their conversion to
boolean masks, for example by simply assigning boolean `True` (in the boolean
mask) to those elements that correspond to non-zero-valued bit fields
(bit fields with at least one bit set to ``1``) or, oftentimes, by assigning
`True` to elements whose corresponding bit fields have only *specific fields*
set (to ``1``). This more sophisticated analysis of bit fields can be
accomplished using *bit masks* and the aforementioned masking operation.

The `~astropy.nddata.bitmask` module provides two functions that facilitate
conversion of bit field arrays (i.e., DQ arrays) to boolean masks:
`~astropy.nddata.bitmask.bitfield_to_boolean_mask` to convert an input bit
fields array to a boolean mask using an input bit mask (or list of individual
bit flags) and `~astropy.nddata.bitmask.interpret_bit_flags` to create bit mask
from input list of individual bit flags.

Creating boolean masks
**********************


Overview
========

`~astropy.nddata.bitmask.bitfield_to_boolean_mask` by default assumes that
all input bit fields that have at least one bit turned "ON" correspond to
"bad" data (i.e., pixels) and converts them to boolean `True` in the output
boolean mask (otherwise output boolean mask values are set to `False`).

Often, for specific algorithms and situations, some bit flags are OK and
can be ignored. `~astropy.nddata.bitmask.bitfield_to_boolean_mask` accepts
lists of bit flags that *by default must be ignored* in the input bit fields
when creating boolean masks.

Fundamentally, *by default*, `~astropy.nddata.bitmask.bitfield_to_boolean_mask`
performs the following operation:

.. _main_eq:

``(1)    boolean_mask = (bitfield & ~bit_mask) != 0``

(here ``&`` is bitwise ``and`` and ``~`` is the bitwise ``not`` operations).
In the previous formula, ``bit_mask`` is a bit mask created from individual
bit flags that need to be ignored in the bit field.

.. _table1:

.. table:: Table 1: Examples of Boolean Mask Computations \
           (default parameters and 8-bit data type)

    +--------------+--------------+--------------+--------------+------------+
    | Bit Field    |  Bit Mask    | ~(Bit Mask)  | Bit Field &  |Boolean Mask|
    |              |              |              | ~(Bit Mask)  |            |
    +==============+==============+==============+==============+============+
    |11011001 (217)|01010000 (80) |10101111 (175)|10001001 (137)|   True     |
    +--------------+--------------+--------------+--------------+------------+
    |11011001 (217)|10101111 (175)|01010000 (80) |01010000 (80) |   True     |
    +--------------+--------------+--------------+--------------+------------+
    |00001001 (9)  |01001001 (73) |10110110 (182)|00000000 (0)  |   False    |
    +--------------+--------------+--------------+--------------+------------+
    |00001001 (9)  |00000000 (0)  |11111111 (255)|00001001 (9)  |   True     |
    +--------------+--------------+--------------+--------------+------------+
    |00001001 (9)  |11111111 (255)|00000000 (0)  |00000000 (0)  |   False    |
    +--------------+--------------+--------------+--------------+------------+


Specifying bit flags
====================

`~astropy.nddata.bitmask.bitfield_to_boolean_mask` accepts either an integer
bit mask or lists of bit flags. Lists of bit flags will be combined into a
bit mask and can be provided either as a Python list of
**integer bit flag values** or as a comma-separated (or ``+``-separated)
list of integer bit flag values. Consider the bit mask from the first example
in `Table 1 <table1_>`_. In this case ``ignore_flags`` can be set either to:

    - an integer value bit mask 80, or
    - a Python list indicating individual non-zero
      *bit flag values:* ``[16, 64]``, or
    - a string of comma-separated *bit flag values*: ``'16,64'``, or
    - a string of ``+``-separated *bit flag values*: ``'16+64'``

For example,

    >>> from astropy.nddata import bitmask
    >>> import numpy as np
    >>> bitmask.bitfield_to_boolean_mask(217, ignore_flags=80)
    array(True...)
    >>> bitmask.bitfield_to_boolean_mask(217, ignore_flags='16,64')
    array(True...)
    >>> bitmask.bitfield_to_boolean_mask(217, ignore_flags=[16, 64])
    array(True...)
    >>> bitmask.bitfield_to_boolean_mask(9, ignore_flags=[1, 8, 64])
    array(False...)
    >>> bitmask.bitfield_to_boolean_mask([9, 10, 73, 217], ignore_flags='1,8,64')
    array([False,  True, False,  True]...)

It is also possible to specify the type of the output mask:

    >>> bitmask.bitfield_to_boolean_mask([9, 10, 73, 217], ignore_flags='1,8,64', dtype=np.uint8)
    array([0, 1, 0, 1], dtype=uint8)


Modifying the Formula for Creating Boolean Masks
================================================

`~astropy.nddata.bitmask.bitfield_to_boolean_mask` provides several parameters
that can be used to modify the formula used to create boolean masks.


Inverting Bit Mask
------------------

Sometimes it is more convenient to be able to specify those bit
flags that *must be considered* when creating the boolean mask and all other
flags should be ignored. In `~astropy.nddata.bitmask.bitfield_to_boolean_mask`
this can be accomplished by setting parameter ``flip_bits`` to `True`.
This effectively modifies `equation (1) <main_eq_>`_ to:

.. _modif_eq2:

``(2)    boolean_mask = (bitfield & bit_mask) != 0``

So, instead of

    >>> bitmask.bitfield_to_boolean_mask([9, 10, 73, 217], ignore_flags=[1, 8, 64])
    array([False,  True, False,  True]...)

one can obtain the same result as

    >>> bitmask.bitfield_to_boolean_mask(
    ...     [9, 10, 73, 217], ignore_flags=[2, 4, 16, 32, 128], flip_bits=True
    ... )
    array([False,  True, False,  True]...)

Note however, when ``ignore_flags`` is a comma-separated list of bit flag
values, ``flip_bits`` cannot be set to neither `True` or `False`. Instead,
to flip bits of the bit mask formed from a string list of comma-separated
bit flag values, one can prepend a single ``~`` to the list:

    >>> bitmask.bitfield_to_boolean_mask([9, 10, 73, 217], ignore_flags='~2+4+16+32+128')
    array([False,  True, False,  True]...)


Inverting Boolean Mask
----------------------

Other times, it may be more convenient to obtain an inverted mask in which
flagged data are converted to `False` instead of `True`:

.. _modif_eq3:

``(3)    boolean_mask = (bitfield & ~bit_mask) == 0``

This can be accomplished by changing ``good_mask_value`` parameter from
its default value (`False`) to `True`. For example,

    >>> bitmask.bitfield_to_boolean_mask([9, 10, 73, 217], ignore_flags=[1, 8, 64],
    ...                                  good_mask_value=True)
    array([ True, False,  True, False]...)
