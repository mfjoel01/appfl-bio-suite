"""Argonne palette and typography shared by coordinator-side scientific plots."""

from __future__ import annotations

import logging
import os
import re
from functools import lru_cache
from pathlib import Path

PALETTE_PATH = Path(__file__).with_name("plot_assets") / "argonne-palette.css"
# Scientific exports use the light palette, independent of the desktop theme.
_light = PALETTE_PATH.read_text().split("@media", 1)[0]
COLORS = dict(re.findall(r"--([\w-]+):\s*(#[\da-fA-F]{6})", _light))
CATEGORICAL = [COLORS[f"chart-{i}"] for i in range(1, 7)]
SEQUENTIAL = [COLORS[f"seq-{i}"] for i in range(1, 8)]
DIVERGING = [
    COLORS[f"div-{k}"] for k in ("neg-3", "neg-2", "neg-1", "mid", "pos-1", "pos-2", "pos-3")
]
ANCESTRY = dict(zip(("EUR", "AFR", "AMR", "EAS", "CSA", "MID"), CATEGORICAL, strict=True))


@lru_cache(maxsize=4)
def configure_font(font_path: str = "") -> str:
    """Register a supplied Helvetica file/directory, or discover installed Helvetica.

    Set BIOSIM_FONT_PATH for privately licensed fonts; font files are never packaged.
    A fallback is logged explicitly; the return value identifies the actual font.
    """
    from matplotlib import font_manager

    if font_path:
        source = Path(font_path)
        if not source.exists():
            raise FileNotFoundError(source)
        files = sorted(source.rglob("*")) if source.is_dir() else [source]
        for path in files:
            if path.suffix.lower() in {".ttf", ".otf", ".ttc"}:
                font_manager.fontManager.addfont(str(path))
    prop = font_manager.FontProperties(family=["Helvetica"])
    try:
        path = font_manager.findfont(prop, fallback_to_default=False)
    except ValueError:
        if os.environ.get("BIOSIM_REQUIRE_HELVETICA") == "1":
            raise RuntimeError(
                "Helvetica unavailable; set BIOSIM_FONT_PATH to licensed font files"
            ) from None
        path = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans"))
        logging.getLogger(__name__).warning(
            "Helvetica unavailable: preview uses DejaVu Sans. Set BIOSIM_FONT_PATH for Helvetica."
        )
    return font_manager.FontProperties(fname=path).get_name()


def apply_style() -> str:
    """Apply brand colors and Helvetica; return the actual rendering font."""
    import matplotlib as mpl
    from cycler import cycler
    from matplotlib.colors import LinearSegmentedColormap

    font = configure_font(os.environ.get("BIOSIM_FONT_PATH", ""))
    name = "argonne_sequential"
    if name not in mpl.colormaps:
        mpl.colormaps.register(LinearSegmentedColormap.from_list(name, SEQUENTIAL))
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font],
            "mathtext.fontset": "custom",
            "mathtext.rm": font,
            "mathtext.it": f"{font}:italic",
            "mathtext.bf": f"{font}:bold",
            "mathtext.cal": font,
            "mathtext.sf": font,
            "mathtext.tt": font,
            "figure.facecolor": COLORS["anl-white"],
            "axes.facecolor": COLORS["anl-white"],
            "savefig.facecolor": COLORS["anl-white"],
            "text.color": COLORS["anl-ink"],
            "axes.labelcolor": COLORS["anl-gray-700"],
            "axes.edgecolor": COLORS["anl-gray-300"],
            "xtick.color": COLORS["anl-gray-500"],
            "ytick.color": COLORS["anl-gray-500"],
            "grid.color": COLORS["anl-gray-300"],
            "axes.prop_cycle": cycler(color=CATEGORICAL),
            "image.cmap": name,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    return font
