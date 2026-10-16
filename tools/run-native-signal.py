"""Launch one owned R process with SIGINT enabled, then interrupt its native loop."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

rscript, script, directory, *library = sys.argv[1:]
out = Path(directory)

def enable_interrupt():
    # A noninteractive shell can pass SIG_IGN to background children. Restore
    # the signal before R starts so this tests Rcpp's actual interrupt polling.
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGINT})

with (out / "worker.log").open("wb") as log:
    process = subprocess.Popen([rscript, "--vanilla", script, "--worker", directory, *library],
        stdout=log, stderr=subprocess.STDOUT, start_new_session=True, preexec_fn=enable_interrupt)
    try:
        deadline = time.monotonic() + 30
        ready = out / "ready.pid"
        while not ready.exists():
            if process.poll() is not None:
                raise RuntimeError("R worker exited before reaching the native call")
            if time.monotonic() > deadline:
                raise TimeoutError("R worker did not become ready")
            time.sleep(.01)
        pid = int(ready.read_text().strip())
        if pid != process.pid:
            raise RuntimeError("Ready PID does not belong to the launched worker")
        time.sleep(.03)
        os.kill(pid, signal.SIGINT)
        code = process.wait(timeout=15)
        (out / "launch.json").write_text(json.dumps(dict(sigint_sent=True, exit_code=code)))
        if code != 0:
            raise RuntimeError(f"R worker exited with status {code}")
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
