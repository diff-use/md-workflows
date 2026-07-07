"""Path-resolution helpers shared by step resolvers."""

from __future__ import annotations

from pathlib import Path


def under(base: Path, path: str | Path) -> Path:
    """Resolve ``path`` against ``base`` unless it is already absolute.

    Used so a relative config value (e.g. ``mdp_dir="artifacts"``) is located inside the
    run's ``workdir`` while an absolute path is honored as-is — keeping steps
    cwd-independent.
    """
    p = Path(path)
    return p if p.is_absolute() else base / p
