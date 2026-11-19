"""
style_boot.py — find the shared house style and expose it, loudly.

Why this exists
---------------
Every plotting script in this directory carried the same block:

    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
        from paper_colors import RETAINED, ELIMINATED
    except Exception:
        RETAINED, ELIMINATED = "#2C6FBB", "#C0392B"

Two independent faults, and the `except` hid both:

1. **The path is wrong.** `paper_colors.py` lives in the analysis repo, which
   sits BESIDE this directory rather than above it —
   `.../miniprot/Ascarid-Molecular-Evolution/lib/` against
   `.../miniprot/07_polymorphism/` — so `parent.parent / "lib"` never resolves.

2. **The names are wrong.** `paper_colors` exports `REGION`, a dict. There is no
   `RETAINED` or `ELIMINATED` symbol in it, so even with the path correct the
   import raises ImportError.

The result: every figure these scripts produced used the fallback hexes rather
than the house palette, and nothing ever said so. `#2C6FBB` is not `#4C72B0`.
Side by side with Figures 1-6 the region blue and red would not have matched,
and the cause is invisible from the output.

So this module fails loudly. If the style cannot be found it warns at WARNING
level and falls back to the *correct* house values copied literally, rather than
to a different palette chosen by whoever wrote the fallback.

Usage
-----
    from style_boot import RETAINED, ELIMINATED          # region colours only
    from style_boot import paper_colors                  # the whole module

Set PAPER_COLORS_DIR to skip the search.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

__all__ = ["RETAINED", "ELIMINATED", "REGION", "paper_colors", "STYLE_SOURCE"]

#: The house values, from paper_colors.py. Duplicated here ONLY as the fallback,
#: so that a failed search degrades to the right colours instead of different
#: ones. If paper_colors ever changes these, this copy is stale and the warning
#: below is the signal to go and look.
_HOUSE = {"retained": "#4C72B0", "eliminated": "#C44E52"}


def find_paper_colors() -> Path | None:
    """Return the directory holding paper_colors.py, or None."""
    here = Path(__file__).resolve().parent
    candidates: list[Path] = []
    env = os.environ.get("PAPER_COLORS_DIR")
    if env:
        candidates.append(Path(env).expanduser())
    for base in (here, here.parent, here.parent.parent):
        candidates += [base, base / "lib", base / "scripts", base / "src"]
        # One level of globbing reaches a sibling repo's lib/ without needing to
        # know the repo's name.
        try:
            candidates += sorted(base.glob("*/lib"))
        except OSError:
            pass
    seen: set[Path] = set()
    for c in candidates:
        if c in seen:
            continue
        seen.add(c)
        if (c / "paper_colors.py").is_file():
            return c
    return None


_dir = find_paper_colors()
paper_colors = None
STYLE_SOURCE = "fallback"

if _dir is not None:
    sys.path.insert(0, str(_dir))
    try:
        import paper_colors as _pc          # type: ignore
        paper_colors = _pc
        STYLE_SOURCE = str(_dir / "paper_colors.py")
    except Exception as exc:                # pragma: no cover
        logging.warning("Found %s but could not import it (%s). Falling back to "
                        "the literal house palette.", _dir / "paper_colors.py",
                        exc)

if paper_colors is not None:
    REGION = dict(paper_colors.REGION)
else:
    REGION = dict(_HOUSE)
    logging.warning(
        "paper_colors.py not found — using the literal house palette. Figures "
        "will have the right colours but none of the house sizing, so do NOT "
        "submit them. Set PAPER_COLORS_DIR to the directory holding it "
        "(the lib/ folder of this repository).")

RETAINED = REGION["retained"]
ELIMINATED = REGION["eliminated"]
