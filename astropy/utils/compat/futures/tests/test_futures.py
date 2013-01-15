import logging
import multiprocessing
import re
import sys
import threading
import time
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

if sys.version_info < (3, 0):
    next = lambda x: x.next()

if sys.platform.startswith('win'):
    import ctypes
    import ctypes.wintypes

from ... import futures
from .._base import (
    PENDING, RUNNING, CANCELLED, CANCELLED_AND_NOTIFIED, FINISHED, Future,
    LOGGER, STDERR_HANDLER)


def create_future(state=PENDING, exception=None, result=None):
    f = Future()
    f._state = state
    f._exception = exception
    f._result = result
    return f


PENDING_FUTURE = create_future(state=PENDING)
RUNNING_FUTURE = create_future(state=RUNNING)
CANCELLED_FUTURE = create_future(state=CANCELLED)
CANCELLED_AND_NOTIFIED_FUTURE = create_future(state=CANCELLED_AND_NOTIFIED)
EXCEPTION_FUTURE = create_future(state=FINISHED, exception=IOError())
SUCCESSFUL_FUTURE = create_future(state=FINISHED, result=42)


def mul(x, y):
    return x * y


class Call(object):
    """A call that can be submitted to a future.Executor for testing.

    The call signals when it is called and waits for an event before finishing.
    """
    CALL_LOCKS = {}

    def _create_event(self):
        if sys.platform.startswith('win'):
            class SECURITY_ATTRIBUTES(ctypes.Structure):
                _fields_ = [("nLength", ctypes.wintypes.DWORD),
                            ("lpSecurityDescriptor", ctypes.wintypes.LPVOID),
                            ("bInheritHandle", ctypes.wintypes.BOOL)]

            s = SECURITY_ATTRIBUTES()
            s.nLength = ctypes.sizeof(s)
            s.lpSecurityDescriptor = None
            s.bInheritHandle = True

            handle = ctypes.windll.kernel32.CreateEventA(ctypes.pointer(s),
                                                         True,
                                                         False,
                                                         None)
            assert handle is not None
            return handle
        else:
            event = multiprocessing.Event()
            self.CALL_LOCKS[id(event)] = event
            return id(event)

    def _wait_on_event(self, handle):
        if sys.platform.startswith('win'):
            r = ctypes.windll.kernel32.WaitForSingleObject(handle, 5 * 1000)
            assert r == 0
        else:
            self.CALL_LOCKS[handle].wait()

    def _signal_event(self, handle):
        if sys.platform.startswith('win'):
            r = ctypes.windll.kernel32.SetEvent(handle)
            assert r != 0
        else:
            self.CALL_LOCKS[handle].set()

    def __init__(self, manual_finish=False, result=42):
        self._called_event = self._create_event()
        self._can_finish = self._create_event()

        self._result = result

        if not manual_finish:
            self._signal_event(self._can_finish)

    def wait_on_called(self):
        self._wait_on_event(self._called_event)

    def set_can(self):
        self._signal_event(self._can_finish)

    def __call__(self):
        self._signal_event(self._called_event)
        self._wait_on_event(self._can_finish)

        return self._result

    def close(self):
        self.set_can()
        if sys.platform.startswith('win'):
            ctypes.windll.kernel32.CloseHandle(self._called_event)
            ctypes.windll.kernel32.CloseHandle(self._can_finish)
        else:
            del self.CALL_LOCKS[self._called_event]
            del self.CALL_LOCKS[self._can_finish]


class ExceptionCall(Call):
    def __call__(self):
        self._signal_event(self._called_event)
        self._wait_on_event(self._can_finish)
        raise ZeroDivisionError()


class MapCall(Call):
    def __init__(self, result=42):
        super(MapCall, self).__init__(manual_finish=True, result=result)

    def __call__(self, manual_finish):
        if manual_finish:
            super(MapCall, self).__call__()
        return self._result


class ExecutorShutdownTest(unittest.TestCase):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    def test_run_after_shutdown(self):
        self.executor.shutdown()
        self.assertRaises(RuntimeError,
                          self.executor.submit,
                          pow, 2, 5)

    def _start_some_futures(self):
        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        call3 = Call(manual_finish=True)

        try:
            self.executor.submit(call1)
            self.executor.submit(call2)
            self.executor.submit(call3)

            call1.wait_on_called()
            call2.wait_on_called()
            call3.wait_on_called()

            call1.set_can()
            call2.set_can()
            call3.set_can()
        finally:
            call1.close()
            call2.close()
            call3.close()


class TestThreadPoolShutdown(ExecutorShutdownTest):
    def setUp(self):
        self.executor = futures.ThreadPoolExecutor(max_workers=5)

    def tearDown(self):
        self.executor.shutdown(wait=True)

    def test_threads_terminate(self):
        self._start_some_futures()
        self.assertEqual(len(self.executor._threads), 3)
        self.executor.shutdown()
        for t in self.executor._threads:
            t.join()

    def test_context_manager_shutdown(self):
        with futures.ThreadPoolExecutor(max_workers=5) as e:
            executor = e
            self.assertEqual(list(e.map(abs, range(-5, 5))),
                             [5, 4, 3, 2, 1, 0, 1, 2, 3, 4])

        for t in executor._threads:
            t.join()

    def test_del_shutdown(self):
        executor = futures.ThreadPoolExecutor(max_workers=5)
        executor.map(abs, range(-5, 5))
        threads = executor._threads
        del executor

        for t in threads:
            t.join()


class TestProcessPoolShutdown(ExecutorShutdownTest):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=5)

    def tearDown(self):
        self.executor.shutdown(wait=True)

    def test_processes_terminate(self):
        self._start_some_futures()
        self.assertEqual(len(self.executor._processes), 5)
        processes = self.executor._processes
        self.executor.shutdown()

        for p in processes:
            p.join()

    def test_context_manager_shutdown(self):
        with futures.ProcessPoolExecutor(max_workers=5) as e:
            executor = e
            self.assertEqual(list(e.map(abs, range(-5, 5))),
                             [5, 4, 3, 2, 1, 0, 1, 2, 3, 4])

        for p in self.executor._processes:
            p.join()

    def test_del_shutdown(self):
        executor = futures.ProcessPoolExecutor(max_workers=5)
        list(executor.map(abs, range(-5, 5)))
        queue_management_thread = executor._queue_management_thread
        processes = executor._processes
        del executor

        queue_management_thread.join()
        for p in processes:
            p.join()


class WaitTests(unittest.TestCase):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    def test_first_completed(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()

        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)

            t = threading.Thread(target=wait_test)
            t.start()
            done, not_done = futures.wait(
                [CANCELLED_FUTURE, future1, future2],
                return_when=futures.FIRST_COMPLETED)

            self.assertEquals(set([future1]), done)
            self.assertEquals(set([CANCELLED_FUTURE, future2]), not_done)
        finally:
            call1.close()
            call2.close()

    def test_first_completed_one_already_completed(self):
        call1 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)

            finished, pending = futures.wait(
                [SUCCESSFUL_FUTURE, future1],
                return_when=futures.FIRST_COMPLETED)

            self.assertEquals(set([SUCCESSFUL_FUTURE]), finished)
            self.assertEquals(set([future1]), pending)
        finally:
            call1.close()

    def test_first_exception(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()
            call2.set_can()

        call1 = Call(manual_finish=True)
        call2 = ExceptionCall(manual_finish=True)
        call3 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)
            future3 = self.executor.submit(call3)

            t = threading.Thread(target=wait_test)
            t.start()
            finished, pending = futures.wait(
                [future1, future2, future3],
                return_when=futures.FIRST_EXCEPTION)

            self.assertEquals(set([future1, future2]), finished)
            self.assertEquals(set([future3]), pending)
        finally:
            call1.close()
            call2.close()
            call3.close()

    def test_first_exception_some_already_complete(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()

        call1 = ExceptionCall(manual_finish=True)
        call2 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)

            t = threading.Thread(target=wait_test)
            t.start()
            finished, pending = futures.wait(
                [SUCCESSFUL_FUTURE,
                 CANCELLED_FUTURE,
                 CANCELLED_AND_NOTIFIED_FUTURE,
                 future1, future2],
                return_when=futures.FIRST_EXCEPTION)

            self.assertEquals(set([SUCCESSFUL_FUTURE,
                                   CANCELLED_AND_NOTIFIED_FUTURE,
                                   future1]), finished)
            self.assertEquals(set([CANCELLED_FUTURE, future2]), pending)
        finally:
            call1.close()
            call2.close()

    def test_first_exception_one_already_failed(self):
        call1 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)

            finished, pending = futures.wait(
                [EXCEPTION_FUTURE, future1],
                return_when=futures.FIRST_EXCEPTION)

            self.assertEquals(set([EXCEPTION_FUTURE]), finished)
            self.assertEquals(set([future1]), pending)
        finally:
            call1.close()

    def test_all_completed(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()
            call2.set_can()

        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)

            t = threading.Thread(target=wait_test)
            t.start()
            finished, pending = futures.wait(
                [future1, future2],
                return_when=futures.ALL_COMPLETED)

            self.assertEquals(set([future1, future2]), finished)
            self.assertEquals(set(), pending)
        finally:
            call1.close()
            call2.close()

    def test_all_completed_some_already_completed(self):
        def wait_test():
            while not future1._waiters:
                pass

            future4.cancel()
            call1.set_can()
            call2.set_can()
            call3.set_can()

        self.assertTrue(
            futures.process.EXTRA_QUEUED_CALLS <= 1,
            'this test assumes that future4 will be cancelled before it is '
            'queued to run - which might not be the case if '
            'ProcessPoolExecutor is too aggresive in scheduling futures')
        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        call3 = Call(manual_finish=True)
        call4 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)
            future3 = self.executor.submit(call3)
            future4 = self.executor.submit(call4)

            t = threading.Thread(target=wait_test)
            t.start()
            finished, pending = futures.wait(
                [SUCCESSFUL_FUTURE,
                 CANCELLED_AND_NOTIFIED_FUTURE,
                 future1, future2, future3, future4],
                return_when=futures.ALL_COMPLETED)

            self.assertEquals(set([SUCCESSFUL_FUTURE,
                                   CANCELLED_AND_NOTIFIED_FUTURE,
                                   future1, future2, future3, future4]),
                              finished)
            self.assertEquals(set(), pending)
        finally:
            call1.close()
            call2.close()
            call3.close()
            call4.close()

    def test_timeout(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()

        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)

            t = threading.Thread(target=wait_test)
            t.start()
            finished, pending = futures.wait(
                [CANCELLED_AND_NOTIFIED_FUTURE,
                 EXCEPTION_FUTURE,
                 SUCCESSFUL_FUTURE,
                 future1, future2],
                timeout=1,
                return_when=futures.ALL_COMPLETED)

            self.assertEquals(set([CANCELLED_AND_NOTIFIED_FUTURE,
                                   EXCEPTION_FUTURE,
                                   SUCCESSFUL_FUTURE,
                                   future1]), finished)
            self.assertEquals(set([future2]), pending)
        finally:
            call1.close()
            call2.close()


class TestThreadPoolWait(WaitTests):
    def setUp(self):
        self.executor = futures.ThreadPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class TestProcessPoolWait(WaitTests):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class AsCompletedTests(unittest.TestCase):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    # TODO(brian@sweetapp.com): Should have a test with a non-zero timeout.
    def test_no_timeout(self):
        def wait_test():
            while not future1._waiters:
                pass
            call1.set_can()
            call2.set_can()

        call1 = Call(manual_finish=True)
        call2 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            future2 = self.executor.submit(call2)

            t = threading.Thread(target=wait_test)
            t.start()
            completed = set(futures.as_completed(
                [CANCELLED_AND_NOTIFIED_FUTURE,
                 EXCEPTION_FUTURE,
                 SUCCESSFUL_FUTURE,
                 future1, future2]))
            self.assertEquals(set(
                [CANCELLED_AND_NOTIFIED_FUTURE,
                 EXCEPTION_FUTURE,
                 SUCCESSFUL_FUTURE,
                 future1, future2]), completed)
        finally:
            call1.close()
            call2.close()

    def test_zero_timeout(self):
        call1 = Call(manual_finish=True)
        try:
            future1 = self.executor.submit(call1)
            completed_futures = set()
            try:
                for future in futures.as_completed(
                        [CANCELLED_AND_NOTIFIED_FUTURE,
                         EXCEPTION_FUTURE,
                         SUCCESSFUL_FUTURE,
                         future1],
                        timeout=0):
                    completed_futures.add(future)
            except futures.TimeoutError:
                pass

            self.assertEquals(set([CANCELLED_AND_NOTIFIED_FUTURE,
                                   EXCEPTION_FUTURE,
                                   SUCCESSFUL_FUTURE]),
                              completed_futures)
        finally:
            call1.close()


class TestThreadPoolAsCompleted(AsCompletedTests):
    def setUp(self):
        self.executor = futures.ThreadPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class TestProcessPoolAsCompleted(AsCompletedTests):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class ExecutorTest(unittest.TestCase):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    # Executor.shutdown() and context manager usage is tested by
    # ExecutorShutdownTest.
    def test_submit(self):
        future = self.executor.submit(pow, 2, 8)
        self.assertEquals(256, future.result())

    def test_submit_keyword(self):
        future = self.executor.submit(mul, 2, y=8)
        self.assertEquals(16, future.result())

    def test_map(self):
        self.assertEqual(
            list(self.executor.map(pow, range(10), range(10))),
            list(map(pow, range(10), range(10))))

    def test_map_exception(self):
        i = self.executor.map(divmod, [1, 1, 1, 1], [2, 3, 0, 5])
        self.assertEqual(next(i), (0, 1))
        self.assertEqual(next(i), (0, 1))
        self.assertRaises(ZeroDivisionError, next, i)

    def test_map_timeout(self):
        results = []
        timeout_call = MapCall()
        try:
            try:
                for i in self.executor.map(timeout_call,
                                           [False, False, True],
                                           timeout=1):
                    results.append(i)
            except futures.TimeoutError:
                pass
            else:
                self.fail('expected TimeoutError')
        finally:
            timeout_call.close()

        self.assertEquals([42, 42], results)


class TestThreadPoolExecutor(ExecutorTest):
    def setUp(self):
        self.executor = futures.ThreadPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class TestProcessPoolExecutor(ExecutorTest):
    def setUp(self):
        self.executor = futures.ProcessPoolExecutor(max_workers=1)

    def tearDown(self):
        self.executor.shutdown(wait=True)


class TestFuture(unittest.TestCase):
    def test_done_callback_with_result(self):
        self.callback_result = None

        def fn(callback_future):
            self.callback_result = callback_future.result()

        f = Future()
        f.add_done_callback(fn)
        f.set_result(5)
        self.assertEquals(5, self.callback_result)

    def test_done_callback_with_exception(self):
        self.callback_exception = None

        def fn(callback_future):
            self.callback_exception = callback_future.exception()

        f = Future()
        f.add_done_callback(fn)
        f.set_exception(Exception('test'))
        self.assertEquals(('test',), self.callback_exception.args)

    def test_done_callback_with_cancel(self):
        self.was_cancelled = None

        def fn(callback_future):
            self.was_cancelled = callback_future.cancelled()

        f = Future()
        f.add_done_callback(fn)
        self.assertTrue(f.cancel())
        self.assertTrue(self.was_cancelled)

    def test_done_callback_raises(self):
        LOGGER.removeHandler(STDERR_HANDLER)
        logging_stream = StringIO()
        handler = logging.StreamHandler(logging_stream)
        LOGGER.addHandler(handler)
        try:
            self.raising_was_called = False
            self.fn_was_called = False

            def raising_fn(callback_future):
                self.raising_was_called = True
                raise Exception('doh!')

            def fn(callback_future):
                self.fn_was_called = True

            f = Future()
            f.add_done_callback(raising_fn)
            f.add_done_callback(fn)
            f.set_result(5)
            self.assertTrue(self.raising_was_called)
            self.assertTrue(self.fn_was_called)
            self.assertTrue('Exception: doh!' in logging_stream.getvalue())
        finally:
            LOGGER.removeHandler(handler)
            LOGGER.addHandler(STDERR_HANDLER)

    def test_done_callback_already_successful(self):
        self.callback_result = None

        def fn(callback_future):
            self.callback_result = callback_future.result()

        f = Future()
        f.set_result(5)
        f.add_done_callback(fn)
        self.assertEquals(5, self.callback_result)

    def test_done_callback_already_failed(self):
        self.callback_exception = None

        def fn(callback_future):
            self.callback_exception = callback_future.exception()

        f = Future()
        f.set_exception(Exception('test'))
        f.add_done_callback(fn)
        self.assertEquals(('test',), self.callback_exception.args)

    def test_done_callback_already_cancelled(self):
        self.was_cancelled = None

        def fn(callback_future):
            self.was_cancelled = callback_future.cancelled()

        f = Future()
        self.assertTrue(f.cancel())
        f.add_done_callback(fn)
        self.assertTrue(self.was_cancelled)

    def test_repr(self):
        self.assertTrue(re.match('<Future at 0x[0-9a-f]+L? state=pending>',
                                 repr(PENDING_FUTURE)))
        self.assertTrue(re.match('<Future at 0x[0-9a-f]+L? state=running>',
                                 repr(RUNNING_FUTURE)))
        self.assertTrue(re.match('<Future at 0x[0-9a-f]+L? state=cancelled>',
                                 repr(CANCELLED_FUTURE)))
        self.assertTrue(re.match('<Future at 0x[0-9a-f]+L? state=cancelled>',
                                 repr(CANCELLED_AND_NOTIFIED_FUTURE)))
        self.assertTrue(re.match(
            '<Future at 0x[0-9a-f]+L? state=finished raised IOError>',
            repr(EXCEPTION_FUTURE)))
        self.assertTrue(re.match(
            '<Future at 0x[0-9a-f]+L? state=finished returned int>',
            repr(SUCCESSFUL_FUTURE)))

    def test_cancel(self):
        f1 = create_future(state=PENDING)
        f2 = create_future(state=RUNNING)
        f3 = create_future(state=CANCELLED)
        f4 = create_future(state=CANCELLED_AND_NOTIFIED)
        f5 = create_future(state=FINISHED, exception=IOError())
        f6 = create_future(state=FINISHED, result=5)

        self.assertTrue(f1.cancel())
        self.assertEquals(f1._state, CANCELLED)

        self.assertFalse(f2.cancel())
        self.assertEquals(f2._state, RUNNING)

        self.assertTrue(f3.cancel())
        self.assertEquals(f3._state, CANCELLED)

        self.assertTrue(f4.cancel())
        self.assertEquals(f4._state, CANCELLED_AND_NOTIFIED)

        self.assertFalse(f5.cancel())
        self.assertEquals(f5._state, FINISHED)

        self.assertFalse(f6.cancel())
        self.assertEquals(f6._state, FINISHED)

    def test_cancelled(self):
        self.assertFalse(PENDING_FUTURE.cancelled())
        self.assertFalse(RUNNING_FUTURE.cancelled())
        self.assertTrue(CANCELLED_FUTURE.cancelled())
        self.assertTrue(CANCELLED_AND_NOTIFIED_FUTURE.cancelled())
        self.assertFalse(EXCEPTION_FUTURE.cancelled())
        self.assertFalse(SUCCESSFUL_FUTURE.cancelled())

    def test_done(self):
        self.assertFalse(PENDING_FUTURE.done())
        self.assertFalse(RUNNING_FUTURE.done())
        self.assertTrue(CANCELLED_FUTURE.done())
        self.assertTrue(CANCELLED_AND_NOTIFIED_FUTURE.done())
        self.assertTrue(EXCEPTION_FUTURE.done())
        self.assertTrue(SUCCESSFUL_FUTURE.done())

    def test_running(self):
        self.assertFalse(PENDING_FUTURE.running())
        self.assertTrue(RUNNING_FUTURE.running())
        self.assertFalse(CANCELLED_FUTURE.running())
        self.assertFalse(CANCELLED_AND_NOTIFIED_FUTURE.running())
        self.assertFalse(EXCEPTION_FUTURE.running())
        self.assertFalse(SUCCESSFUL_FUTURE.running())

    def test_result_with_timeout(self):
        self.assertRaises(futures.TimeoutError,
                          PENDING_FUTURE.result, timeout=0)
        self.assertRaises(futures.TimeoutError,
                          RUNNING_FUTURE.result, timeout=0)
        self.assertRaises(futures.CancelledError,
                          CANCELLED_FUTURE.result, timeout=0)
        self.assertRaises(futures.CancelledError,
                          CANCELLED_AND_NOTIFIED_FUTURE.result, timeout=0)
        self.assertRaises(IOError, EXCEPTION_FUTURE.result, timeout=0)
        self.assertEqual(SUCCESSFUL_FUTURE.result(timeout=0), 42)

    def test_result_with_success(self):
        # TODO(brian@sweetapp.com): This test is timing dependant.
        def notification():
            # Wait until the main thread is waiting for the result.
            time.sleep(1)
            f1.set_result(42)

        f1 = create_future(state=PENDING)
        t = threading.Thread(target=notification)
        t.start()

        self.assertEquals(f1.result(timeout=5), 42)

    def test_result_with_cancel(self):
        # TODO(brian@sweetapp.com): This test is timing dependant.
        def notification():
            # Wait until the main thread is waiting for the result.
            time.sleep(1)
            f1.cancel()

        f1 = create_future(state=PENDING)
        t = threading.Thread(target=notification)
        t.start()

        self.assertRaises(futures.CancelledError, f1.result, timeout=5)

    def test_exception_with_timeout(self):
        self.assertRaises(futures.TimeoutError,
                          PENDING_FUTURE.exception, timeout=0)
        self.assertRaises(futures.TimeoutError,
                          RUNNING_FUTURE.exception, timeout=0)
        self.assertRaises(futures.CancelledError,
                          CANCELLED_FUTURE.exception, timeout=0)
        self.assertRaises(futures.CancelledError,
                          CANCELLED_AND_NOTIFIED_FUTURE.exception, timeout=0)
        self.assertTrue(isinstance(EXCEPTION_FUTURE.exception(timeout=0),
                                   IOError))
        self.assertEqual(SUCCESSFUL_FUTURE.exception(timeout=0), None)

    def test_exception_with_success(self):
        def notification():
            # Wait until the main thread is waiting for the exception.
            time.sleep(1)
            with f1._condition:
                f1._state = FINISHED
                f1._exception = IOError()
                f1._condition.notify_all()

        f1 = create_future(state=PENDING)
        t = threading.Thread(target=notification)
        t.start()

        self.assertTrue(isinstance(f1.exception(timeout=5), IOError))
