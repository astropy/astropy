from pathlib import Path


def test_include():
    from astropy.wcs import get_include

    include = get_include()
    assert type(include) is str
    include_path = Path(include)
    assert include_path.is_dir()
    assert include_path.joinpath("astropy_wcs_api.h").is_file()
