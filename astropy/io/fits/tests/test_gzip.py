from astropy.io import fits
import os
import shutil
import astropy
from astropy.table import Table

def test_gzip():
	print astropy.__file__
	data_dir = os.path.join(os.path.dirname(__file__), 'data')
	shutil.copy(os.path.join(data_dir, 'table.fits'), os.path.join(data_dir, 'table2.fits'))
	os.system('gzip '+os.path.join(data_dir, 'table2.fits'))
	try:
		Table.read(os.path.join(data_dir, 'table2.fits.gz'))
		passed = True
	except ValueError as e:
		passed = False
	os.remove(os.path.join(data_dir, 'table2.fits.gz'))
	if not passed:
		raise(e)

	
