"""Standard library only. The guarded import is optional and must be ignored."""
import json
import sys

try:
    import riversand  # optional: absence is handled below
except ImportError:
    riversand = None
