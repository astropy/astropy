# Licensed under a 3-clause BSD style license - see LICENSE.rst

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from threading import Barrier
from typing import Generic, TypeVar
from uuid import uuid4

import pytest

from astropy.config.paths import _DirType

N_THREADS = 10

T = TypeVar("T")


@dataclass(slots=True, frozen=True, kw_only=True)
class Result(Generic[T]):
    actual: T
    expected: T

    def match(self) -> bool:
        return self.actual == self.expected


def closure_base(
    dirtype: _DirType,
    *,
    directory: Path,
    barrier: Barrier,
) -> Result[Path]:
    local_path = directory / str(uuid4())
    local_path.mkdir()

    barrier.wait()
    with dirtype.legacy_context_manager(local_path):
        res = dirtype.path_getter("astropy")

    return Result(actual=res, expected=local_path / "astropy")


def assert_valid_results(results: list[Result[T]]) -> None:
    __tracebackhide__ = True
    assert len(results) == N_THREADS
    assert len(set(results)) == N_THREADS
    assert all(r.match() for r in results)


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_set_temp_dir(dirtype, tmp_path):
    closure = partial(
        closure_base,
        dirtype,
        directory=tmp_path,
        barrier=Barrier(N_THREADS),
    )
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = [executor.submit(closure) for _ in range(N_THREADS)]

    assert_valid_results([f.result() for f in futures])


@pytest.mark.parametrize("dt_in", _DirType)
@pytest.mark.parametrize("dt_out", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_nesting(dt_out: _DirType, dt_in: _DirType, tmp_path: Path):
    # check that nesting doesn't deadlock
    barrier = Barrier(N_THREADS)

    def closure() -> Result[tuple[Path, Path, Path, Path]]:
        local_path_1 = tmp_path / str(uuid4())
        local_path_1.mkdir()
        local_path_2 = tmp_path / str(uuid4())
        local_path_2.mkdir()

        barrier.wait()
        with dt_out.legacy_context_manager(local_path_1):
            res0 = dt_out.path_getter()
            with dt_in.legacy_context_manager(local_path_2):
                res1 = dt_in.path_getter()
                res2 = dt_out.path_getter()
            res3 = dt_out.path_getter()
        return Result(
            actual=(res0, res1, res2, res3),
            expected=(
                local_path_1 / "astropy",
                local_path_2 / "astropy",
                (local_path_1 if dt_in != dt_out else local_path_2) / "astropy",
                local_path_1 / "astropy",
            ),
        )

    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = [executor.submit(closure) for _ in range(N_THREADS)]

    assert_valid_results([f.result() for f in futures])


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_mixed_settings(tmp_path):
    barrier = Barrier(N_THREADS)

    cache_setter = partial(
        closure_base, _DirType.CACHE, directory=tmp_path, barrier=barrier
    )
    config_setter = partial(
        closure_base, _DirType.CONFIG, directory=tmp_path, barrier=barrier
    )

    closures = [cache_setter if n % 2 else config_setter for n in range(N_THREADS)]
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = [executor.submit(c) for c in closures]

    assert_valid_results([f.result() for f in futures])
