import numpy as np
import pytest

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
	Test conversion to and from leap year
	'''
	time_decimalyear = Time('2008-02-29', format = 'iso').decimalyear
	assert round(time_decimalyear%1.0, 8) == round((31.0 + 28.0)/366.0, 8) #the very beginning of day 29
	time_datetime = Time(2008.1613, format = 'decimalyear').datetime
	assert time_datetime.day == 29
	
def test_last_day_of_year():
	'''
	Test last day of year is handled correctly
	'''
	t1 = Time('2012-12-31 23:59:59.0', format = 'iso').decimalyear
	assert t1 < 2013.0, 'Rounding to next year'
	
def test_integer_input():
	'''
	Test that everything behaves correctly when decimalyear is an integer
	'''
	t1 = Time(2010, format = 'decimalyear')
	t2 = Time(2010.0, format = 'decimalyear')
	assert t1 == t2, 'integer test failed'
	
def test_float32_input():
	'''
	Test that everything behaves correctly when a float32 is input
	'''
	t1 = Time(np.float32(2011.0), format = 'decimalyear').decimalyear
	t2 = Time(2011.0, format = 'decimalyear').decimalyear
	assert t1 == t2, 'float 32 test failed'
	
def test_list_input():
	'''
	Test that conversion works for list object. Just make sure it runs
	'''
	
	t_list = [2010.0, 2011.2, 2014.1]
	decimalyear_list = Time(t_list, format = 'decimalyear')
	iso_list = decimalyear_list.iso
	
def test_array_input():
	'''
	Test that conversion works for np.array object. Just make sure it runs
	'''
	t_list = np.array([2010.0, 2011.2, 2014.1])
	decimalyear_list = Time(t_list, format = 'decimalyear')
	iso_list = decimalyear_list.iso
	
def test_string_input():
	'''
	Test that if a string is input, a Time object is not formed
	'''
	with pytest.raises(ValueError):
		Time('2010.3', format = 'decimalyear')
