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

---

## In-Progress Work: Issue #10776 — `Quantity.__array_ufunc__` NotImplemented fix

**Branch:** `scouting-report-priority-2-1692904929402517997`

**Venv:** `.venv/` in repo root (created with `python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"`).
Run tests with `.venv/bin/pytest`.

### What was done

**Fix in `astropy/units/quantity.py` (line ~690):**
Removed `None` from the `ignored_ufunc` tuple in the `except` block of `Quantity.__array_ufunc__`.

Before:
```python
ignored_ufunc = (
    None,
    np.ndarray.__array_ufunc__,
    type(self).__array_ufunc__,
)
```
After:
```python
ignored_ufunc = (
    np.ndarray.__array_ufunc__,
    type(self).__array_ufunc__,
)
```

**Why:** When an external object (like `Time`) has no `__array_ufunc__` defined, `getattr(type(obj), "__array_ufunc__", None)` returns `None`. Previously `None` was in `ignored_ufunc`, so the caught exception was re-raised instead of returning `NotImplemented`. Removing `None` from the tuple causes `Quantity.__array_ufunc__` to correctly return `NotImplemented` for all non-ndarray/non-Quantity objects.

**Tests in `astropy/units/tests/test_quantity_ufuncs.py`:**

1. Added `TestUfuncFallbackForUnknownUnitClass` class with two tests:
   - `test_radd_fallback_no_array_ufunc` — object has `unit=u.s` (incompatible with meters for add), no `__array_ufunc__`; verifies `__radd__` is NOT invoked but... **this test currently FAILS** (see below).
   - `test_rmul_fallback_array_ufunc_none` — object has `__array_ufunc__ = None`; verifies `__rmul__` is called. **PASSES** (but see caveat below).

2. Updated the existing `TestUfuncReturnsNotImplemented::TestBinaryUfuncs::test_basic` regex to also match numpy's new "all returned NotImplemented" error message.

### Current test status

```
.venv/bin/pytest astropy/units/tests/test_quantity_ufuncs.py::TestUfuncFallbackForUnknownUnitClass
.venv/bin/pytest astropy/units/tests/test_quantity_ufuncs.py::TestUfuncReturnsNotImplemented
```

- `test_rmul_fallback_array_ufunc_none` — PASSES
- `test_radd_fallback_no_array_ufunc` — FAILS with `TypeError: operand type(s) all returned NotImplemented from __array_ufunc__(...)`
- `TestUfuncReturnsNotImplemented` — all PASS (with updated regex)

### Key insight the next instance needs to understand

**The `test_rmul_fallback_array_ufunc_none` test passes WITHOUT the fix too.** When an object has `__array_ufunc__ = None`, numpy raises `TypeError` *before* calling `Quantity.__array_ufunc__` at all. NumPy's C-level ndarray operator catches this TypeError and returns `NotImplemented` to Python, which then tries `obj.__rmul__`. So our fix doesn't affect this path.

**The real behavioral change** from our fix: for objects *without* `__array_ufunc__` (like `DuckQuantity1`, `DuckQuantity2`, `Time`), Quantity's `__array_ufunc__` used to re-raise the caught exception. Now it returns `NotImplemented`, causing numpy to raise "all returned NotImplemented" TypeError. Either way an exception is raised — only the message changes. The fallback to `__radd__`/`__rmul__` does NOT work via this path because numpy's "all returned NotImplemented" TypeError propagates through the C-level operator (unlike the `__array_ufunc__ = None` path where C catches it).

### What still needs to be done

1. **Fix or replace `test_radd_fallback_no_array_ufunc`**: The current test is wrong — returning `NotImplemented` from Quantity's `__array_ufunc__` does NOT enable `__radd__` fallback for objects without `__array_ufunc__`. Options:
   - **Option A**: Call `Quantity.__array_ufunc__` *directly* and assert the return value is `NotImplemented` (instead of testing end-to-end arithmetic). This directly tests the fix.
   - **Option B**: Test end-to-end using an object with `unit` that causes the try block to fail AND has `__array_ufunc__ = None` (not just "not defined"), so numpy's C-level catches the TypeError and falls back. But as shown, when `__array_ufunc__ = None`, Quantity's `__array_ufunc__` isn't even called.
   - **Recommended**: Use Option A — directly test the return value of `Quantity.__array_ufunc__` when called with an external object. Example:
     ```python
     class ObjNoUfunc:
         unit = u.s
     obj = ObjNoUfunc()
     q = 1.0 * u.m
     result = u.Quantity.__array_ufunc__(q, np.add, '__call__', q, obj)
     assert result is NotImplemented
     ```

2. **Add a changelog entry**: `docs/changes/units/<PR_NUMBER>.bugfix.rst` — needs actual PR number from GitHub.

3. **Verify the full test suite** still passes:
   ```bash
   .venv/bin/pytest astropy/units/tests/test_quantity_ufuncs.py -x
   .venv/bin/pytest astropy/units/tests/test_quantity.py -x
   ```

4. **Consider if the fix is sufficient**: Our fix changes the error message for external objects without `__array_ufunc__`. It also correctly makes `Quantity.__array_ufunc__` return `NotImplemented` when called directly (e.g., by `super().__array_ufunc__()` from a subclass). This is the correct and important use case that the scout report describes.
