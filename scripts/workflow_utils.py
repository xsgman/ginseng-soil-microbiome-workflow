from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    tomllib = None


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.toml", help="Path to TOML config file.")
    return parser.parse_args()


def load_config(path: str | Path) -> dict:
    if tomllib is None:
        raise RuntimeError("Python 3.11+ is required for TOML config parsing, or install tomli and adapt this script.")
    config_path = Path(path)
    if not config_path.is_absolute():
        config_path = ROOT / config_path
    with config_path.open("rb") as handle:
        return tomllib.load(handle)


def cfg_path(config: dict, key: str) -> Path:
    path = Path(config["paths"][key])
    return path if path.is_absolute() else ROOT / path


def output_dir(config: dict, subdir: str = "") -> Path:
    out = cfg_path(config, "output_dir")
    if subdir:
        out = out / subdir
    out.mkdir(parents=True, exist_ok=True)
    return out


def setup_style(font_size: float = 8.5) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": font_size,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save_figure(fig, out_base: Path, dpi: int = 450) -> None:
    fig.savefig(f"{out_base}.png", dpi=dpi, bbox_inches="tight")
    fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")


def normalize_soil_metric(text: str) -> str:
    text = str(text).strip()
    mapping = {
        "pH": "pH",
        "Electrical conductivity": "EC",
        "Exchangeable total acidity": "Total acidity",
        "Exchangeable H+": "Exch. H+",
        "Exchangeable Al3+": "Exch. Al3+",
        "Soil organic carbon": "SOC",
        "Organic matter": "OM",
        "Readily oxidizable C": "ROC",
        "Alkali-hydrolyzable N": "AN",
        "Available P": "AP",
        "NH4+-N": "NH4-N",
        "NO3--N": "NO3-N",
        "Cation exchange capacity": "CEC",
        "Urease (S-UE)": "Urease",
    }
    if "glucosidase" in text.lower():
        return "BG"
    return mapping.get(text, text)


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    numeric = df.apply(pd.to_numeric, errors="coerce")
    return (numeric - numeric.mean(axis=0)) / numeric.std(axis=0, ddof=1).replace(0, np.nan)


def require_columns(df: pd.DataFrame, columns: list[str], source: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"{source} missing required columns: {missing}")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)
