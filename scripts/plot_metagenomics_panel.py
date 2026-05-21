from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

from workflow_utils import cfg_path, load_config, output_dir, parse_args, save_figure, setup_style


TREATMENT_COLORS = {"CK": "#59616C", "KXH": "#2F86A6", "RHB": "#B86B2B", "Other": "#D9DEE3"}
MONTH_MARKERS = {"Aug": "o", "Sep": "^"}
CONTRAST_ORDER = [
    "KXH-Aug vs CK-Aug",
    "RHB-Aug vs CK-Aug",
    "KXH-Sep vs CK-Sep",
    "RHB-Sep vs CK-Sep",
]


def infer_treatment(group: str) -> str:
    text = str(group)
    for treatment in ["CK", "KXH", "RHB"]:
        if text.startswith(treatment):
            return treatment
    return "Other"


def infer_timepoint(group: str) -> str:
    text = str(group)
    if "1" in text or "Aug" in text:
        return "Aug"
    if "2" in text or "Sep" in text:
        return "Sep"
    return "Aug"


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    setup_style(7.6)

    scores = pd.read_csv(cfg_path(config, "pcoa_scores"))
    if "Treatment" not in scores.columns:
        group_col = "Group" if "Group" in scores.columns else "SampleID"
        scores["Treatment"] = scores[group_col].map(infer_treatment)
    if "Timepoint" not in scores.columns:
        group_col = "Group" if "Group" in scores.columns else "SampleID"
        scores["Timepoint"] = scores[group_col].map(infer_timepoint)
    if "SampleID" not in scores.columns:
        scores["SampleID"] = scores.iloc[:, 0]

    centers_path = cfg_path(config, "pcoa_group_centers")
    centers = pd.read_csv(centers_path) if centers_path.exists() else pd.DataFrame()
    if not centers.empty:
        if "Treatment" not in centers.columns:
            group_col = "Group" if "Group" in centers.columns else centers.columns[0]
            centers["Treatment"] = centers[group_col].map(infer_treatment)
        if "Timepoint" not in centers.columns:
            group_col = "Group" if "Group" in centers.columns else centers.columns[0]
            centers["Timepoint"] = centers[group_col].map(infer_timepoint)

    functions = pd.read_csv(cfg_path(config, "function_log2fc"))
    if "feature" not in functions.columns:
        functions = functions.rename(columns={functions.columns[0]: "feature"})
    for col in CONTRAST_ORDER:
        functions[col] = pd.to_numeric(functions[col], errors="coerce")
    functions["max_abs"] = functions[CONTRAST_ORDER].abs().max(axis=1)
    functions = functions.sort_values("max_abs", ascending=False).head(18).reset_index(drop=True)

    fig = plt.figure(figsize=(13.2, 6.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.55], wspace=0.24)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    for treatment in ["Other", "CK", "KXH", "RHB"]:
        sub_t = scores[scores["Treatment"].eq(treatment)]
        if sub_t.empty:
            continue
        for month, marker in MONTH_MARKERS.items():
            sub = sub_t[sub_t["Timepoint"].eq(month)]
            if sub.empty:
                continue
            ax1.scatter(
                sub["PCo1"],
                sub["PCo2"],
                s=28 if treatment == "Other" else 46,
                marker=marker,
                color=TREATMENT_COLORS[treatment],
                edgecolor="white",
                linewidth=0.4,
                alpha=0.30 if treatment == "Other" else 0.92,
            )

    for _, row in centers.iterrows():
        treatment = row["Treatment"]
        if treatment not in {"CK", "KXH", "RHB"}:
            continue
        ax1.scatter(
            row["PCo1"],
            row["PCo2"],
            s=120,
            marker=MONTH_MARKERS.get(row.get("Timepoint", "Aug"), "o"),
            facecolors="none",
            edgecolors=TREATMENT_COLORS[treatment],
            linewidth=1.4,
        )

    ax1.set_title("A   Bray-Curtis PCoA", loc="left", weight="bold")
    ax1.set_xlabel("PCo1")
    ax1.set_ylabel("PCo2")
    ax1.grid(color="#EDF0F2", lw=0.7)
    ax1.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="None", color=TREATMENT_COLORS[t], label=t)
            for t in ["CK", "KXH", "RHB", "Other"]
        ],
        title="Treatment",
        loc="upper left",
    )

    values = functions[CONTRAST_ORDER].to_numpy(dtype=float)
    limit = np.nanmax(np.abs(values))
    cmap = LinearSegmentedColormap.from_list("log2fc", ["#2D6F95", "#F7F5F0", "#B85C2F"])
    im = ax2.imshow(values, cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect="auto")
    ax2.set_title("B   Functional potentials relative to matched CK", loc="left", weight="bold")
    ax2.set_yticks(np.arange(len(functions)))
    ax2.set_yticklabels(functions["feature"])
    ax2.set_xticks(np.arange(len(CONTRAST_ORDER)))
    ax2.set_xticklabels(CONTRST_ORDER if False else CONTRAST_ORDER, rotation=40, ha="right")
    ax2.set_xticks(np.arange(-0.5, len(CONTRAST_ORDER), 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(functions), 1), minor=True)
    ax2.grid(which="minor", color="white", lw=0.8)
    ax2.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax2, fraction=0.035, pad=0.02)
    cbar.set_label("log2FC")

    out = output_dir(config, "metagenomics")
    save_figure(fig, out / "metagenomics_pcoa_function_panel")
    with pd.ExcelWriter(out / "metagenomics_pcoa_function_source_data.xlsx") as writer:
        scores.to_excel(writer, "PCoA_scores", index=False)
        centers.to_excel(writer, "PCoA_centers", index=False)
        functions.to_excel(writer, "Function_log2FC", index=False)
    print(f"Saved outputs to {out}")


if __name__ == "__main__":
    main()
