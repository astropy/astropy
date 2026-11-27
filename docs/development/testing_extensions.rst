========================================
Testing C and Cython Extensions Directly
========================================

Welcome to the low-level testing guide for Astropy's compiled extensions!
If you're reading this, you are probably about to write tests for C or Cython modules.

While Astropy's public API is heavily tested, our underlying compiled building blocks (the "hot-paths")
need their own isolated testing. This guide explains how to operate one layer below the public API,
right at the boundary between Python and compiled code. Doing this ensures our C-engine safely handles
memory sizing, index mapping, and string parsing on its own—an essential prerequisite for milestones
like the APE split and making our core architecture rock solid.

Why We Test the Raw Boundary
============================

Our primary goal here is to test raw Cython and C-extensions in pure isolation. We want to completely
bypass high-level Python wrappers (like ``astropy.table.Column``). By stripping away the Python wrappers,
we can directly verify that the compiled engine correctly handles memory allocation, index mapping,
and boolean masking without the Python layer hiding potential regressions.

Think of this test layer as the bridge between the high-level user API and the low-level C code.

Testing Architecture & Setup Guide
==================================

1. Bypass the APIs
------------------
Construct the minimal raw numpy arrays, structs, or byte streams required to trigger the compiled logic directly.

2. Isolate C-slots
----------------------------
When testing C-slots (like ``__getitem__``), use minimal shim classes and ``.view()`` casting to map
directly onto raw memory buffers rather than using standard Astropy objects.

3. Target the Boundary
----------------------------
Differentiate between pure C math engines and the Python-toC wrapper. The tests and type stubs should
specifically target the wrapper boundary where Python objects are unpacked into C memory pointers.


Writing the Tests
=================

Data Structures and Parametrization
-----------------------------------

To maintain clean and readable continuous integration (CI) logs, avoid using large tuples in ``@pytest.mark.parametrize``. Instead, use strictly typed dataclasses to hold expected inputs, outputs, and boolean masks.

.. code-block:: python

    from dataclasses import dataclass

    @dataclass(kw_only=True, slots=True, frozen=True)
    class ExpectedResults:
        ...

Ensure every parametrized test uses explicit IDs via ``pytest.param(..., id="name_of_case")``.


Strict Typing and Assertions
============================

1. Type Assertions (No ``isinstance``)
--------------------------------------
Avoid using ``isinstance()`` for return type validation at the C-boundary, as standard array slicing
can strip subclasses while masked array slicing preserves them.
Use strict ``type(x) is y`` assertions to verify the exact return type expected from the C-engine.

2. Strict Type Hinting
----------------------
Do not use generic ``list`` or ``np.ndarray`` type hints. You must use ``numpy.typing.NDArray`` with
exact 1D shapes and dtypes (e.g., ``np.ndarray[tuple[int], np.dtype[np.intp]]``).

3. Semantic Hinting
-------------------
When passing dictionaries that the C-engine will strictly read and not mutate, use ``collections.abc.Mapping``
instead of ``dict``. This signals your exact intent to static type checkers.

Data Validation and Edge Cases
==============================

* **Exhaustive Validation:** Don't just check memory boundaries (e.g., ``assert type(event) is tuple``). Implement deep, strict equality checks to prove the C-engine's internal calculations (like coordinate tuples or state flags) are completely flawless.
* **Strict Iterables:** Whenever using ``zip()`` to compare expected versus actual runtime events, you must include ``strict=True`` to prevent the silent truncation of mismatched lists.
* **Isolate Edge Cases:** Keep standard tests (where clean set-theory math applies) separated from edge cases like Cartesian expansions (which cause explosive O(N^2) memory allocations) or ``np.nan`` identities (which break standard logical assertions).
* **Memory Efficiency:** If you need a boolean mask array of all ``False``, don't type out ``[False, False]``. Use ``np.full(N, False)`` for better memory efficiency.

Catching C-Level Exceptions
===========================

When testing how the C-engine handles truncated byte streams or malformed data, it should safely raise an error across the Python boundary rather than crashing the interpreter.

Use standard Python terminology for exceptions (e.g., the engine "raises" a ``ValueError``). Because C-level errors often prefix unpredictable memory addresses or line numbers, use substring matching rather than exact string comparisons to validate the error:

.. code-block:: python

    import pytest

    with pytest.raises(ValueError) as exc_info:
        # Trigger C-extension error
        ...

    assert "no element found" in str(exc_info.value).lower()

Type Stubs (.pyi) Architecture
==============================
Static type checkers (like Pyright or mypy) cannot easily natively parse ``.pyx`` files.
We bridge this gap using ``.pyi`` stub files.

* **Locate the True Boundary:** Distinguish between pure C math engines and the actual Python wrapper. We only type the wrapper.
* **Argument Enforcement:** If a hand-written C-API parser doesn't use a ``kwlist``, it doesn't inherently support keyword arguments. You must use the ``/`` marker in your ``.pyi`` file to enforce positional-only arguments.
* **Dimensionality Locking:** If the C-engine expects multiple arrays that must match in dimensions (like an image and a kernel), use bounded ``TypeVar`` aliases mapped to ``np.ndarray`` to mathematically enforce dimensionality matching at the static level.
* **Linter Compliance:** To pass Ruff's private export rules (PYI001), any internal ``TypeVar`` or ``TypeAlias`` used strictly for the stub file must be prefixed with an underscore (e.g., ``_D = TypeVar("_D")``).
