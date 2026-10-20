# -*- coding: utf-8 -*-
import concurrent.futures as cf
import os
import signal
import subprocess
import sys
import time

import psutil
import pytest

from automlsa2.helper_functions import kill_child_processes, worker_init


def test_kill_child_processes():
    # sh -c keeps a grandchild (sleep) running under the child
    proc = subprocess.Popen(['sh', '-c', 'sleep 60 & wait'])
    time.sleep(0.2)
    descendants = psutil.Process().children(recursive=True)
    assert len(descendants) >= 2
    kill_child_processes()
    proc.wait(timeout=10)
    assert not any(p.is_running() and p.status() != 'zombie' for p in descendants)


def _interrupt_self():
    os.kill(os.getpid(), signal.SIGINT)
    time.sleep(0.1)
    return os.getpid()


@pytest.mark.skipif(sys.platform == 'win32', reason='POSIX signals')
def test_worker_ignores_sigint():
    with cf.ProcessPoolExecutor(1, initializer=worker_init) as pool:
        pid = pool.submit(_interrupt_self).result(timeout=30)
    assert pid != os.getpid()


def test_worker_init_restores_default_for_exec():
    old = signal.getsignal(signal.SIGINT)
    try:
        worker_init()
        # A Python-level handler (not SIG_IGN) is reset to the default on exec,
        # so external programs started by workers still respond to Ctrl-C.
        assert callable(signal.getsignal(signal.SIGINT))
        out = subprocess.run(
            ['sh', '-c', 'kill -INT $$; echo survived'],
            stdout=subprocess.PIPE,
            text=True,
        )
        assert out.returncode != 0
        assert 'survived' not in out.stdout
    finally:
        signal.signal(signal.SIGINT, old)
