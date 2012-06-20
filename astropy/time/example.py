import astrotime

times = ['1999-01-01 00:00:00.123456789', '2010-01-01 00:00:00']
t = astrotime.Time(times, format='iso', system='utc')
t
t.jd1
t.jd2

print 'Set system to tai'
t.set_system('tai')

t
t.jd1
t.jd2

t.tt  # system property returns a new Time object
t.cxcsec  # format property returns an array or scalar of that representation

# Use properties to convert systems and formats

t = astrotime.Time('2010-01-01 00:00:00', format='iso', system='utc')
t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the transformation
t.jd        # JD representation of time in current system (UTC)
t.iso       # ISO representation of time in current system (UTC)
t.tt.iso    # ISO representation of time in TT system
t.tai.iso   # ISO representation of time in TAI system
t.utc.jd    # JD representation of time in UTC system
t.ut1.jd    # JD representation of time in UT1 system
t.tcg.isot  # ISO time with a "T" in the middle
t.unix      # seconds since 1970.0 (utc) excluding leapseconds
t.cxcsec    # SI seconds since 1998.0 (tt)

t.precision = 9
t.iso

# Use a time scale that requires latitude and longitude
lat = 19.48125
lon = -155.933222
hm = 0.0

t = astrotime.Time('2006-01-15 21:24:37.5', format='iso', system='utc',
                   precision=6, lat=lat, lon=lon)
t.set_delta_ut1_utc(0.3341)  # Explicitly set one part of the transformation
t.utc.iso
t.ut1.iso
t.tai.iso
t.tt.iso
t.tcg.iso
t.tdb.iso
t.tcb.iso
