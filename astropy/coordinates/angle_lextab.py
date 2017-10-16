# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

_tabversion   = '3.8'
_lextokens    = set(('UINT', 'SIMPLE_UNIT', 'DEGREE', 'MINUTE', 'HOUR', 'COLON', 'UFLOAT', 'SIGN', 'SECOND'))
_lexreflags   = 0
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.\\d*)|(\\.\\d+))([eE][+-−]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+−-])|(?P<t_SIMPLE_UNIT>(?:karcsec)|(?:uarcsec)|(?:Earcmin)|(?:Zdeg)|(?:crad)|(?:cycle)|(?:hectoradian)|(?:Yarcmin)|(?:kiloarcsecond)|(?:zeptoarcminute)|(?:adeg)|(?:darcmin)|(?:ddeg)|(?:exaradian)|(?:parcsec)|(?:yoctoradian)|(?:arcsecond)|(?:petadegree)|(?:petaarcminute)|(?:microarcsecond)|(?:mas)|(?:parcmin)|(?:hdeg)|(?:narcmin)|(?:attodegree)|(?:kilodegree)|(?:zettaradian)|(?:fdeg)|(?:zeptoradian)|(?:microradian)|(?:Gdeg)|(?:hectodegree)|(?:attoarcsecond)|(?:Marcmin)|(?:exadegree)|(?:femtodegree)|(?:yottaradian)|(?:pdeg)|(?:zarcmin)|(?:kiloarcminute)|(?:urad)|(?:teraarcsecond)|(?:nrad)|(?:carcsec)|(?:Pdeg)|(?:Yrad)|(?:yrad)|(?:picoarcsecond)|(?:aarcsec)|(?:dekaradian)|(?:Zrad)|(?:femtoradian)|(?:yarcsec)|(?:arcmin)|(?:arcsec)|(?:yottadegree)|(?:drad)|(?:dekadegree)|(?:zdeg)|(?:zeptoarcsecond)|(?:farcmin)|(?:Parcmin)|(?:decaarcminute)|(?:nanoarcminute)|(?:nanoarcsecond)|(?:Tdeg)|(?:decaarcsecond)|(?:nanodegree)|(?:farcsec)|(?:femtoarcminute)|(?:microdegree)|(?:deciarcsecond)|(?:deciarcminute)|(?:attoradian)|(?:dadeg)|(?:decidegree)|(?:hectoarcminute)|(?:milliarcsecond)|(?:femtoarcsecond)|(?:megaarcminute)|(?:yoctoarcminute)|(?:zrad)|(?:hectoarcsecond)|(?:frad)|(?:centiarcsecond)|(?:carcmin)|(?:Garcmin)|(?:decadegree)|(?:Grad)|(?:petaarcsecond)|(?:gigaarcsecond)|(?:megaradian)|(?:Tarcsec)|(?:Prad)|(?:zettadegree)|(?:yottaarcminute)|(?:mrad)|(?:yottaarcsecond)|(?:exaarcminute)|(?:harcmin)|(?:dekaarcsecond)|(?:cy)|(?:ndeg)|(?:teraradian)|(?:teradegree)|(?:Zarcsec)|(?:gigadegree)|(?:Mdeg)|(?:Mrad)|(?:centiarcminute)|(?:uarcmin)|(?:picoradian)|(?:radian)|(?:ydeg)|(?:milliarcminute)|(?:deciradian)|(?:narcsec)|(?:Trad)|(?:picodegree)|(?:yoctodegree)|(?:zettaarcminute)|(?:daarcmin)|(?:arcminute)|(?:yarcmin)|(?:kdeg)|(?:Earcsec)|(?:Edeg)|(?:harcsec)|(?:rad)|(?:centidegree)|(?:Garcsec)|(?:marcsec)|(?:megaarcsecond)|(?:attoarcminute)|(?:cdeg)|(?:Erad)|(?:kiloradian)|(?:daarcsec)|(?:Parcsec)|(?:megadegree)|(?:millidegree)|(?:centiradian)|(?:uas)|(?:teraarcminute)|(?:prad)|(?:yoctoarcsecond)|(?:hrad)|(?:picoarcminute)|(?:petaradian)|(?:Marcsec)|(?:marcmin)|(?:Tarcmin)|(?:zeptodegree)|(?:Yarcsec)|(?:gigaarcminute)|(?:Zarcmin)|(?:arad)|(?:karcmin)|(?:darcsec)|(?:exaarcsecond)|(?:nanoradian)|(?:udeg)|(?:zarcsec)|(?:Ydeg)|(?:decaradian)|(?:milliradian)|(?:aarcmin)|(?:zettaarcsecond)|(?:darad)|(?:microarcminute)|(?:mdeg)|(?:dekaarcminute)|(?:krad)|(?:gigaradian))|(?P<t_MINUTE>m(in(ute(s)?)?)?|′|\\\'|ᵐ)|(?P<t_SECOND>s(ec(ond(s)?)?)?|″|\\"|ˢ)|(?P<t_DEGREE>d(eg(ree(s)?)?)?|°)|(?P<t_HOUR>hour(s)?|h(r)?|ʰ)|(?P<t_COLON>:)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_SIMPLE_UNIT', 'SIMPLE_UNIT'), (None, 'MINUTE'), None, None, None, (None, 'SECOND'), None, None, None, (None, 'DEGREE'), None, None, None, (None, 'HOUR'), None, None, (None, 'COLON')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
