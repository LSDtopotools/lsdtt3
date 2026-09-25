"""Starts a process directly, bypassing lsdtt_test_env."""
import subprocess

subprocess.run(["true"], check=False)
