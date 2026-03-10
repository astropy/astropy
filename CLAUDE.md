# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Installation

Astropy contains C and Cython extensions that must be compiled before use:

```bash
# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Build in-place (for C/Cython extension development)
python setup.py build_ext --inplace
```

## Running Tests

```bash
# Run all tests
pytest --pyargs astropy

# Run tests for a specific subpackage
pytest astropy/coordinates/

# Run a single test file
pytest astropy/coordinates/tests/test_sky_coordinate.py

# Run a single test by name
pytest astropy/coordinates/tests/test_sky_coordinate.py::TestSkyCoord::test_init_empty

# Run with remote data (tests marked @pytest.mark.remote_data)
pytest --remote-data astropy/

# Run doctests only
pytest --doctest-rst docs/

# Run with coverage
pytest --pyargs astropy --cov astropy --cov-report xml:coverage.xml

# Run via tox (e.g., recommended deps, Python 3.12)
tox -e py312-test-recdeps
```

Tests treat warnings as errors by default (configured in `pyproject.toml`). Tests that require network access are skipped unless `--remote-data` is passed.

## Linting and Code Style

```bash
# Run ruff linter
ruff check astropy/

# Run ruff formatter check
ruff format --check astropy/

# Auto-fix linting issues
ruff check --fix astropy/

# Run codestyle tox env
tox -e codestyle
```

Ruff is configured in `pyproject.toml` (base rules) and `.ruff.toml` (additional ignores). The project uses numpy-style docstrings.

## Changelog Entries

Every PR requires a changelog fragment in `docs/changes/<subpackage>/`. Files are named `<PR_NUMBER>.<TYPE>.rst` where type is one of: `feature`, `api`, `bugfix`, `perf`, `other`. Type `other` is not allowed in subpackage subdirectories.

Example: `docs/changes/coordinates/12345.bugfix.rst`

## Architecture Overview

Astropy is structured as a collection of semi-independent subpackages under `astropy/`:

- **`units`** — Physical units and quantities. The `Quantity` class wraps arrays with unit metadata. Most numerical subpackages depend on this.
- **`coordinates`** — Sky/spatial coordinate systems and transformations. Uses a frame-based architecture: `SkyCoord` is the user-facing class, backed by frame classes in `builtin_frames/`. The `representation/` subpackage handles spherical, cartesian, and other representations.
- **`time`** — Time representations and conversions across many formats (ISO, JD, MJD, etc.) and scales (UTC, TAI, TDB, etc.).
- **`table`** — Flexible tabular data with heterogeneous columns. Supports mixin columns (e.g., `Quantity`, `SkyCoord`, `Time` as table columns). The `io/` subpackage provides unified I/O.
- **`modeling`** — Framework for defining and fitting 1D/2D models. Models are composable (additive, multiplicative, compound). Fitters are separate from models.
- **`io/fits`** — FITS file reading/writing (the `astropy.io.fits` module, previously standalone as PyFITS).
- **`io/ascii`** — Reading/writing ASCII tables (CSV, ECSV, VOTable, etc.).
- **`io/registry`** — Unified I/O registry that powers `Table.read()`, `Table.write()` etc. via format auto-detection.
- **`cosmology`** — Standard cosmological models (ΛCDM, wCDM, etc.) with I/O support.
- **`wcs`** — World Coordinate System transformations, wrapping the C WCSLIB library.
- **`nddata`** — N-dimensional data containers with uncertainty, mask, and WCS.
- **`stats`** — Statistical functions (sigma clipping, biweight, etc.).
- **`visualization`** — Image normalization and stretch for display; includes `ZScaleInterval`, `astropy_mpl_style`.
- **`utils`** — Shared infrastructure used across subpackages: `data` (file download/caching), `decorators`, `exceptions`, `metadata`, `console`, IERS data tables.
- **`constants`** — Physical and astronomical constants as `Constant` objects (subclass of `Quantity`).

### Key Design Patterns

- **Lazy subpackage imports**: Subpackages are lazily imported at the `astropy` top level via `__init__.py` using `importlib`.
- **I/O Registry**: All read/write operations go through `astropy.io.registry`. New formats register via `io_registry.register_reader/writer/identifier`.
- **Mixin columns**: `astropy/table/mixins/` defines how non-standard column types (Quantity, SkyCoord, Time, etc.) serialize into tables.
- **C/Cython extensions**: Performance-critical code lives in `.c`/`.pyx` files alongside the Python code. `cextern/` contains vendored C libraries (ERFA, CFITSIO expat).
- **Coordinate frames**: New coordinate frames subclass `BaseCoordinateFrame` and define `frame_specific_representation_info` and transformation functions registered via `frame_transform_graph`.

### Testing Conventions

- Tests live in `tests/` subdirectories within each subpackage.
- `conftest.py` files at multiple levels provide fixtures.
- Remote data tests use `@pytest.mark.remote_data`.
- Image comparison tests use `@pytest.mark.mpl_image_compare`.
- `pytest-doctestplus` runs doctests in `.rst` docs files and module docstrings.
- `filterwarnings = ["error"]` in `pyproject.toml` means unhandled warnings fail tests.
