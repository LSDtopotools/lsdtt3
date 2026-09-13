#!/usr/bin/env python3
"""Run one command and report its peak RSS, CPU time and wall time.

    python rusage_launcher.py RESULT_JSON -- cmd arg ...

Used by run_benchmark_suite.py's run_with_rusage(). Not meant to be run by hand.

WHY THIS FILE EXISTS: THE PEAK-RSS INHERITANCE TRAP

On Linux a child process inherits its parent's LIFETIME PEAK resident set
size, and wait4(2) reports that inherited figure as the child's own peak. At
fork the new mm_struct is a copy of the parent's, high-water mark included; at
execve the old mm's high-water mark is folded into the process's
signal->maxrss, which survives the exec. So a child reports
max(parent's peak at fork, its own real peak).

Measured on WSL2 (kernel 6.6): a Python process touched 1.5 GiB and freed it,
bringing its CURRENT RSS back to 9.5 MiB. `/bin/true` launched from it then
reported a peak of 1509.5 MiB. The same `/bin/true` launched from a fresh
Python process reported 9.0 MiB.

This silently destroyed a whole server campaign's memory column. The benchmark
driver computes fixture statistics (float64 copies of an 8000^2 DEM plus a
RichDEM fill) before measuring anything, reaching a 3.5 GB peak. From then on
441 of 504 rows -- every engine, every fixture below 8000^2 -- reported exactly
3556.98 MiB, and only rows whose true peak exceeded that floor showed a real
number.

The fix is one more process level. This launcher is exec'd fresh, so its own
peak is only what a small Python interpreter needs, and the MEASURED command is
ITS child. The measured command therefore inherits this launcher's small peak,
not the driver's large one. That small peak is still a floor; the driver
measures it (`rss_floor_mb`, by launching /bin/true through here) and flags any
row within reach of it.

KEEP THIS FILE SMALL. Every module imported here raises the floor for every
measurement.

The measured command inherits this process's stdin, stdout and stderr, so the
driver reads its output exactly as if it had launched it directly. The result
goes to a file rather than stdout so that it can never interleave with the
command's own output.

The launcher puts itself in no new session: the driver starts it as a session
leader, so killing that process group on a timeout kills the measured command
too.
"""

import json
import os
import sys
import time


def main() -> int:
    if len(sys.argv) < 4 or sys.argv[2] != "--":
        sys.stderr.write("usage: rusage_launcher.py RESULT_JSON -- cmd ...\n")
        return 2
    result_path = sys.argv[1]
    cmd = sys.argv[3:]

    t0 = time.perf_counter()
    try:
        pid = os.posix_spawnp(cmd[0], cmd, os.environ)
    except OSError as exc:
        sys.stderr.write(f"rusage_launcher: cannot run {cmd[0]!r}: {exc}\n")
        return 127
    _, status, ru = os.wait4(pid, 0)
    wall = time.perf_counter() - t0

    # Imported only after the command has been spawned, so it cannot add to
    # the peak the command inherited.
    import resource

    with open(result_path, "w") as fh:
        json.dump({
            "wall_seconds": wall,
            # kibibytes on Linux
            "peak_rss_kib": ru.ru_maxrss,
            "cpu_seconds": ru.ru_utime + ru.ru_stime,
            "exit_code": os.waitstatus_to_exitcode(status),
            # This launcher's own peak: an upper bound on what the measured
            # command could have inherited from it.
            "launcher_peak_rss_kib":
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        }, fh)

    code = os.waitstatus_to_exitcode(status)
    # A negative code means killed by a signal; report it the way a shell would.
    return code if code >= 0 else 128 - code


if __name__ == "__main__":
    sys.exit(main())
