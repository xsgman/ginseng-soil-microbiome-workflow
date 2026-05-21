from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

from workflow_utils import cfg_path, load_config, output_dir, parse_args, save_figure, setup_style


METRICS = [
    "Disease index",
    "Root diameter",
    "Root fresh weight",
    "Root dry weight",
    "Stem fresh weight",
    "Stem dry weight",
]
COLORS = {"KXH": "#2F86A6", "RHB": "#B86B2B"}
MARKERS = {"Aug": "o", "Sep": "^"}


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    setup_style(8.8)
    workbook = cfg_path(config, "ginseng_workbook")
    months = config["plot"].get("months_ginseng", ["Aug", "Sep"])
    treatments = config["plot"].get("contrast_treatments", ["KXH", "RHB"])

    summary = pd.read_excel(workbook, sheet_name="CK_KXH_RHB_Summary")
    data = summary[summary["metric"].isin(METRICS) & summary["month"].isin(months) & summary["treatment"].isin(treatments)].copy()
    data["benefit_log2_vs_ck"] = pd.to_numeric(data["benefit_log2_vs_ck"], errors="coerce")

    try:
        stats = pd.read_excel(workbook, sheet_name="StatsFocus")
        stats["p_value"] = pd.to_numeric(stats["p_value"], errors="coerce")
        data = data.merge(stats[["metric", "month", "treatment", "p_value"]], on=["metric", "month", "treatment"], how="left")
    except Exception:
        data["p_value"] = pd.NA
    data["significant"] = data["p_value"].lt(config["plot"].get("p_value_cutoff", 0.05))

    y_map = {metric: i for i, metric in enumerate(METRICS[::-1])}
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for metric, y in y_map.items():
        ax.axhline(y, color="#ECEFF1", lw=0.7, zorder=0)
    ax.axvline(0, color="#B8BEC6", lw=1.0)

    for _, row in data.iterrows():
        y = y_map[row["metric"]]
        y += -0.10 if row["treatment"] == treatments[0] else 0.10
        ax.scatter(
            row["benefit_log2_vs_ck"],
            y,
            s=60,
            marker=MARKERS.get(row["month"], "o"),
            color=COLORS.get(row["treatment"], "#888888"),
            edgecolor="#222222" if row["significant"] else "white",
            linewidth=1.2 if row["significant"] else 0.6,
            zorder=3,
        )
        if row["significant"]:
            ax.text(row["benefit_log2_vs_ck"] + 0.04, y + 0.12, "*", fontsize=10, weight="bold")

    ax.set_yticks(list(y_map.values()))
    ax.set_yticklabels(list(y_map.keys()))
    ax.set_xlabel("Benefit log2 response relative to CK")
    ax.set_title("Ginseng phenotype responses relative to CK", loc="left", weight="bold")
    ax.grid(axis="x", color="#EDF0F2", lw=0.7)

    treatment_handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=COLORS[t], markeredgecolor="none", label=t)
        for t in treatments
    ]
    marker_handles = [
        Line2D([], [], marker=MARKERS[m], color="#666666", linestyle="None", label=m) for m in months if m in MARKERS
    ]
    sig_handle = Line2D([], [], marker="o", color="none", markerfacecolor="white", markeredgecolor="#222222", label="p < 0.05")
    leg1 = ax.legend(handles=treatment_handles, title="Treatment", loc="upper right")
    ax.add_artist(leg1)
    ax.legend(handles=marker_handles + [sig_handle], title="Month/stat", loc="lower right")

    out = output_dir(config, "ginseng")
    save_figure(fig, out / "ginseng_phenotype_effect")
    data.to_excel(out / "ginseng_phenotype_effect_source_data.xlsx", index=False)
    print(f"Saved outputs to {out}")


if __name__ == "__main__":
    main()
