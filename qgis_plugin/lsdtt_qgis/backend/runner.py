import subprocess


def run_cli(command, feedback=None, timeout=None, env=None, cwd=None):
    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
        cwd=cwd,
        check=False,
    )

    if feedback is not None:
        if process.stdout:
            for line in process.stdout.splitlines():
                feedback.pushInfo(line)
        if process.stderr:
            for line in process.stderr.splitlines():
                feedback.reportError(line)
    return process


def write_log(log_path, command, process):
    lines = []
    lines.append("command: {}".format(" ".join(command)))
    lines.append("returncode: {}".format(process.returncode))
    lines.append("")
    lines.append("[stdout]")
    lines.append(process.stdout or "")
    lines.append("")
    lines.append("[stderr]")
    lines.append(process.stderr or "")
    with open(log_path, "w", encoding="utf-8") as dst:
        dst.write("\n".join(lines))
