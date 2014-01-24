Test:
	leap year
	first day of year
	last day of year
	integer
	float
	list
	single number
	numpy array
	string (should fail)

def test_first_day():
	'''
	Test that this works for the first day of the year 2012.0
	'''
	assert t1 = Time(2012.0, format = 'decimalyear', scale = 'utc'), 'Time object not initiating with first day of year'
	yday_split = t1.yday.split()
	assert float(yday_split[0]) == 2012.0, 'Year is not 2012'
	assert float(yday_split[1]) == 1.0, 'Day is not 001'
	assert float(yday_split[2]) == 0.0, 'hour is not 0'
	assert float(yday_split[3]) == 0.0, 'minute is not 0'
	assert float(yday_split[4]) == 0.0, 'second is not 0'
	