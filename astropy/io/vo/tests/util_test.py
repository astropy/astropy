# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

"""
A set of tests for the util.py module
"""
from __future__ import absolute_import, print_function

# THIRD-PARTY
from numpy.testing import assert_array_equal, assert_raises

#LOCAL
from astropy.io.vo import util

def test_range_list():
    assert util.coerce_range_list_param((5,)) == ("5.0", 1)

def test_range_list2():
    assert util.coerce_range_list_param((5e-7,8e-7)) == ("5e-07,8e-07", 2)

def test_range_list3():
    assert util.coerce_range_list_param((5e-7,8e-7,"FOO")) == ("5e-07,8e-07;FOO", 3)

def test_range_list4():
    def raises():
        util.coerce_range_list_param((5e-7,(None,8e-7),(4,None),(4,5),"J","FOO"))
    assert_raises(ValueError, raises)

    print(util.coerce_range_list_param((5e-7,(None,8e-7),(4,None),(4,5),"J","FOO"), numeric=False))
    assert util.coerce_range_list_param((5e-7,(None,8e-7),(4,None),(4,5),"J","FOO"), numeric=False) == ("5e-07,/8e-07,4/,4/5,J;FOO", 6)

def test_range_list5():
    def raises():
        util.coerce_range_list_param(('FOO',))
    assert_raises(ValueError, raises)

def test_range_list6():
    def raises():
        print(util.coerce_range_list_param((5,'FOO'), util.stc_reference_frames))
    assert_raises(ValueError, raises)

def test_range_list7():
    assert util.coerce_range_list_param(("J",), numeric=False) == ("J", 1)

def test_range_list8():
    for s in ["5.0",
              "5e-07,8e-07",
              "5e-07,8e-07;FOO",
              "5e-07,/8e-07,4.0/,4.0/5.0;FOO",
              "J"]:
        assert util.coerce_range_list_param(s, numeric=False)[0] == s

def test_range_list9():
    def raises():
        util.coerce_range_list_param("52,-27.8;FOO", util.stc_reference_frames)
    assert_raises(ValueError, raises)

    assert util.coerce_range_list_param("52,-27.8;GALACTIC", util.stc_reference_frames)

