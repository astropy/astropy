#Test:
#	leap year
#	last day of year
#	integer
#	float
#	list
#	single number
#	numpy array
#	string (should fail)

from astropy.time import Time
def test_basic_invocation():
	t1 = Time(2012.0, format = 'decimalyear'), 'cannot initialize decimalyear object'
	
def test_basic_conversion():
	t1 = Time('2014-01-23', format = 'iso')
	t1.decimalyear

def test_first_day():
	'''
	Test that this works for the first day of the year 2012.0
	'''
	t1 = Time(2012.0, format = 'decimalyear')
	yday_split = t1.yday.split(':')
	assert float(yday_split[0]) == 2012.0, 'Year is not 2012'
	assert float(yday_split[1]) == 1.0, 'Day is not 001'
	assert float(yday_split[2]) == 0.0, 'hour is not 0'
	assert float(yday_split[3]) == 0.0, 'minute is not 0'
	assert float(yday_split[4]) == 0.0, 'second is not 0'
	
def test_leap_year():
	'''
	Test that leap year is handled correctly
	'''
	time_decimalyear = Time('2008-02-29')
	