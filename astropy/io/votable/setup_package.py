# Licensed under a 3-clause BSD style license - see LICENSE.rst

from pathlib import Path

from setuptools import Extension


def get_extensions(build_type="release"):
    VO_DIR = Path(__file__).parent.joinpath("src").relative_to(Path.cwd())

    return [
        Extension(
            "astropy.io.votable.tablewriter",
            [str(VO_DIR / "tablewriter.c")],
            include_dirs=[str(VO_DIR)],
        )
    ]
