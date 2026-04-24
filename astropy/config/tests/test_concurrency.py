# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from threading import Barrier
from typing import Generic, TypeAlias, TypeVar
from uuid import uuid4

import pytest

from astropy.config import (
    get_cache_dir_path,
    get_config_dir_path,
    set_temp_cache,
    set_temp_config,
)
from astropy.config.paths import _SetTempPath

N_THREADS = 10

T = TypeVar("T")


@dataclass(slots=True, frozen=True, kw_only=True)
class Result(Generic[T]):
    actual: T
    expected: T

    def match(self) -> bool:
        return self.actual == self.expected


PathGetter: TypeAlias = Callable[[], Path]


def getter_from_manager(cm: type[_SetTempPath]) -> PathGetter:
    # this is functionally equivalent to an immutable dict
    # it could be refactored into a frozendict when Python 3.14 is unsupported
    if cm is set_temp_cache:
        return get_cache_dir_path
    elif cm is set_temp_config:
        return get_config_dir_path
    else:
        raise ValueError


def closure_base(
    ctx_manager: type[_SetTempPath],
    *,
    directory: Path,
    barrier: Barrier,
) -> Result[Path]:
    local_path = directory / str(uuid4())
    local_path.mkdir()
    getter = getter_from_manager(ctx_manager)

    barrier.wait()
    with ctx_manager(local_path):
        res = getter()

    return Result(actual=res, expected=local_path / "astropy")


def assert_valid_results(results: list[Result[T]]) -> None:
    __tracebackhide__ = True
    assert len(results) == N_THREADS
    assert len(set(results)) == N_THREADS
    assert all(r.match() for r in results)


@pytest.mark.parametrize("ctx_manager", [set_temp_cache, set_temp_config])
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_set_temp_dir(ctx_manager, tmp_path):
    closure = partial(
        closure_base,
        ctx_manager,
        directory=tmp_path,
        barrier=Barrier(N_THREADS),
    )
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = [executor.submit(closure) for _ in range(N_THREADS)]

    assert_valid_results([f.result() for f in futures])


@pytest.mark.parametrize(
    "cm_out, cm_in",
    [
        pytest.param(set_temp_cache, set_temp_cache, id="cache-cache"),
        pytest.param(set_temp_cache, set_temp_config, id="cache-config"),
        pytest.param(set_temp_config, set_temp_cache, id="config-cache"),
        pytest.param(set_temp_config, set_temp_config, id="config-config"),
    ],
)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_nesting(cm_out: type[_SetTempPath], cm_in: type[_SetTempPath], tmp_path: Path):
    # check that nesting doesn't deadlock
    barrier = Barrier(N_THREADS)
    getter_out = getter_from_manager(cm_out)
    getter_in = getter_from_manager(cm_in)

    def closure() -> Result[tuple[Path, Path, Path, Path]]:
        local_path_1 = tmp_path / str(uuid4())
        local_path_1.mkdir()
        local_path_2 = tmp_path / str(uuid4())
        local_path_2.mkdir()

        barrier.wait()
        with cm_out(local_path_1):
            res0 = getter_out()
            with cm_in(local_path_2):
                res1 = getter_in()
                res2 = getter_out()
            res3 = getter_out()
        return Result(
            actual=(res0, res1, res2, res3),
            expected=(
                local_path_1 / "astropy",
                local_path_2 / "astropy",
                (local_path_1 if cm_in != cm_out else local_path_2) / "astropy",
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
        closure_base, set_temp_cache, directory=tmp_path, barrier=barrier
    )
    config_setter = partial(
        closure_base, set_temp_config, directory=tmp_path, barrier=barrier
    )

    closures = [cache_setter if n % 2 else config_setter for n in range(N_THREADS)]
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        futures = [executor.submit(c) for c in closures]

    assert_valid_results([f.result() for f in futures])
