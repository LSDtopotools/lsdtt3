"""Private helpers for staged Python keyword migrations."""

from __future__ import annotations

import warnings
from typing import Any


MISSING = object()


def resolve_deprecated_keyword(
    canonical_name: str,
    canonical_value: Any,
    legacy_name: str,
    legacy_value: Any,
    *,
    stacklevel: int = 3,
) -> Any:
    """Resolve one canonical/legacy keyword pair without changing its value."""
    if canonical_value is not MISSING and legacy_value is not MISSING:
        raise TypeError(
            f"{canonical_name!r} and deprecated {legacy_name!r} cannot both be supplied"
        )
    if legacy_value is not MISSING:
        warnings.warn(
            f"{legacy_name!r} is deprecated; use {canonical_name!r}",
            DeprecationWarning,
            stacklevel=stacklevel,
        )
        return legacy_value
    return canonical_value
