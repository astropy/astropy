from ..wcs import transform_coord_meta_from_wcs


# def test_coord_type_from_ctype():
#     assert transform_coord_meta_from_wcs(' LON') == ('longitude', None, None)
#     assert transform_coord_meta_from_wcs(' LAT') == ('latitude', None, None)
#     assert transform_coord_meta_from_wcs('HPLN') == ('longitude', u.arcsec, 180.)
#     assert transform_coord_meta_from_wcs('HPLT') == ('latitude', u.arcsec, None)
#     assert transform_coord_meta_from_wcs('RA--') == ('longitude', u.hourangle, None)
#     assert transform_coord_meta_from_wcs('DEC-') == ('latitude', None, None)
#     assert transform_coord_meta_from_wcs('spam') == ('scalar', None, None)
