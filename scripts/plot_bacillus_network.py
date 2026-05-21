from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from workflow_utils import cfg_path, load_config, output_dir, parse_args, save_figure, setup_style


DISPLAY_LABELS = {
    "Bacillus": "Bacillus",
    "Fusarium": "Fusarium",
    "Ilyonectria": "Ilyonectria",
    "Sphingomonas": "Sphingomonas",
    "Stenotrophomonas": "Stenotrophomonas",
    "Cellvibrio": "Cellvibrio",
    "Achromobacter": "Achromobacter",
}
NODE_TYPES = {
    "Bacillus": "center",
    "Fusarium": "pathogen",
    "Ilyonectria": "pathogen",
    "Sphingomonas": "response",
    "Stenotrophomonas": "response",
    "Cellvibrio": "response",
    "Achromobacter": "response",
}
NODE_COLORS = {"center": "#3FA66B", "pathogen": "#C85C5C", "response": "#43A6A0"}
EDGE_COLORS = {"positive": "#D58A2F", "negative": "#3478B8"}
POSITIONS = {
    "Bacillus": (0.0, 0.0),
    "Fusarium": (2.25, 0.95),
    "Ilyonectria": (2.25, -0.95),
    "Sphingomonas": (-2.15, 0.95),
    "Stenotrophomonas": (-2.15, -0.10),
    "Cellvibrio": (-1.20, -1.45),
    "Achromobacter": (0.95, -1.45),
}


def read_genus_abundance(path) -> pd.DataFrame:
    df = pd.read_csv(path)
    lower = {c.lower(): c for c in df.columns}
    if "genus" in lower:
        genus_col = lower["genus"]
        matrix = df.set_index(genus_col)
    else:
        sample_like = [c for c in df.columns if c.lower() in {"sample", "sampleid", "sample_id"}]
        matrix = df.drop(columns=sample_like).T if sample_like else df.T
    matrix.index = matrix.index.astype(str)
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    if x.nunique() < 2 or y.nunique() < 2:
        return math.nan, math.nan
    rho, p = spearmanr(x.to_numpy(dtype=float), y.to_numpy(dtype=float), nan_policy="omit")
    return float(rho), float(p)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    setup_style(9.0)
    p_cutoff = config["plot"].get("p_value_cutoff", 0.05)
    rho_cutoff = config["plot"].get("rho_cutoff", 0.35)
    abundance = read_genus_abundance(cfg_path(config, "genus_abundance"))

    required = [g for g in POSITIONS if g in abundance.index]
    if "Bacillus" not in required:
        raise ValueError("Bacillus must be present in genus_abundance.csv")

    rows = []
    for genus in required:
        if genus == "Bacillus":
            continue
        rho, p = spearman(abundance.loc["Bacillus"], abundance.loc[genus])
        if not math.isnan(rho) and p < p_cutoff and abs(rho) >= rho_cutoff:
            rows.append(
                {
                    "source": "Bacillus",
                    "target": genus,
                    "rho": rho,
                    "p_value": p,
                    "edge_type": "positive" if rho >= 0 else "negative",
                }
            )

    if "Fusarium" in abundance.index and "Ilyonectria" in abundance.index:
        rho, p = spearman(abundance.loc["Fusarium"], abundance.loc["Ilyonectria"])
        if not math.isnan(rho) and p < p_cutoff and abs(rho) >= rho_cutoff:
            rows.append(
                {
                    "source": "Fusarium",
                    "target": "Ilyonectria",
                    "rho": rho,
                    "p_value": p,
                    "edge_type": "positive" if rho >= 0 else "negative",
                }
            )

    edges = pd.DataFrame(rows)
    if edges.empty:
        nodes = ["Bacillus"]
    else:
        nodes = ["Bacillus"] + sorted(set(edges["source"]).union(set(edges["target"])) - {"Bacillus"})

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    ax.axis("off")
    ax.set_xlim(-3.3, 3.5)
    ax.set_ylim(-2.1, 2.1)
    ax.set_title("Bacillus-centered candidate validation network", loc="left", weight="bold", fontsize=14)

    for _, edge in edges.iterrows():
        x1, y1 = POSITIONS[edge["source"]]
        x2, y2 = POSITIONS[edge["target"]]
        width = 1.4 + 4.5 * (abs(edge["rho"]) - rho_cutoff) / max(1e-9, 1 - rho_cutoff)
        ax.plot([x1, x2], [y1, y2], color=EDGE_COLORS[edge["edge_type"]], lw=width, alpha=0.82, solid_capstyle="round")

    node_rows = []
    for node in nodes:
        if node not in POSITIONS:
            continue
        node_type = NODE_TYPES.get(node, "response")
        size = 2100 if node == "Bacillus" else 900
        ax.scatter(*POSITIONS[node], s=size, color=NODE_COLORS[node_type], edgecolor="white", lw=1.8, zorder=3)
        x, y = POSITIONS[node]
        label = DISPLAY_LABELS.get(node, node)
        if node == "Bacillus":
            ax.text(x, y, label, ha="center", va="center", color="white", weight="bold", fontsize=11.5, zorder=4)
        else:
            ha = "left" if x > 0 else "right"
            dx = 0.28 if x > 0 else -0.28
            ax.text(x + dx, y, label, ha=ha, va="center", fontsize=9.8, color="#252A31")
        if edges.empty:
            hit = pd.DataFrame()
        else:
            hit = edges[(edges["source"].eq("Bacillus") & edges["target"].eq(node)) | (edges["target"].eq("Bacillus") & edges["source"].eq(node))]
        node_rows.append(
            {
                "genus": node,
                "display_label": DISPLAY_LABELS.get(node, node),
                "node_type": node_type,
                "mean_abundance": abundance.loc[node].mean(),
                "degree": int((edges["source"].eq(node) | edges["target"].eq(node)).sum()) if not edges.empty else 0,
                "Bacillus_edge_rho": float(hit.iloc[0]["rho"]) if not hit.empty else np.nan,
                "Bacillus_edge_p": float(hit.iloc[0]["p_value"]) if not hit.empty else np.nan,
            }
        )

    ax.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="None", markerfacecolor=NODE_COLORS["center"], markeredgecolor="white", label="Bacillus"),
            Line2D([], [], marker="o", linestyle="None", markerfacecolor=NODE_COLORS["pathogen"], markeredgecolor="white", label="Pathogen-related"),
            Line2D([], [], marker="o", linestyle="None", markerfacecolor=NODE_COLORS["response"], markeredgecolor="white", label="Rhizosphere response"),
            Line2D([], [], color=EDGE_COLORS["positive"], lw=3, label="Positive"),
            Line2D([], [], color=EDGE_COLORS["negative"], lw=3, label="Negative"),
        ],
        loc="upper right",
    )
    ax.text(
        -3.2,
        -1.95,
        "Network is for qPCR/plate-count target prioritization; correlation does not prove causality.",
        fontsize=8.5,
        color="#66707C",
    )

    out = output_dir(config, "metagenomics")
    save_figure(fig, out / "bacillus_validation_network")
    pd.DataFrame(node_rows).to_csv(out / "bacillus_validation_network_nodes.csv", index=False)
    edges.to_csv(out / "bacillus_validation_network_edges.csv", index=False)
    print(f"Saved outputs to {out}")


if __name__ == "__main__":
    main()
