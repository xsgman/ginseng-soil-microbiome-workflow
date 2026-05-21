from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

from workflow_utils import cfg_path, load_config, normalize_soil_metric, output_dir, parse_args, save_figure, setup_style


METRIC_GROUPS = {
    "Acid-related": ["pH", "EC", "Total acidity", "Exch. H+", "Exch. Al3+"],
    "Carbon": ["SOC", "OM", "ROC"],
    "Nutrients": ["AN", "AP", "NH4-N", "NO3-N", "CEC"],
    "Enzymes": ["BG", "Urease"],
}
GROUP_COLORS = {
    "Acid-related": "#8EA6B4",
    "Carbon": "#C9A05F",
    "Nutrients": "#8FB27E",
    "Enzymes": "#B58EBE",
}


def metric_order() -> list[str]:
    ordered: list[str] = []
    for metrics in METRIC_GROUPS.values():
        ordered.extend(metrics)
    return ordered


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    setup_style(8.8)

    months = config["plot"].get("months_soil", ["Jun", "Jul", "Aug", "Sep"])
    treatments = config["plot"].get("treatments_focus", ["CK", "KXH", "RHB"])
    contrasts = [t for t in config["plot"].get("contrast_treatments", ["KXH", "RHB"]) if t != "CK"]

    workbook = cfg_path(config, "soil_workbook")
    summary = pd.read_excel(workbook, sheet_name="Summary")
    reps = pd.read_excel(workbook, sheet_name="ReplicateMeans")
    summary["metric_label"] = summary["metric_en"].map(normalize_soil_metric)
    reps["metric_label"] = reps["metric_en"].map(normalize_soil_metric)

    keep_summary = summary[
        summary["phase"].eq("treatment")
        & summary["month"].isin(months)
        & summary["treatment"].isin(treatments)
        & summary["metric_label"].isin(metric_order())
    ].copy()
    keep_reps = reps[
        reps["phase"].eq("treatment")
        & reps["month"].isin(months)
        & reps["treatment"].isin(treatments)
        & reps["metric_label"].isin(metric_order())
    ].copy()

    pooled = {}
    skipped = []
    for metric in metric_order():
        values = pd.to_numeric(keep_reps.loc[keep_reps["metric_label"].eq(metric), "value"], errors="coerce").dropna()
        if len(values) < 2:
            values = pd.to_numeric(keep_summary.loc[keep_summary["metric_label"].eq(metric), "mean"], errors="coerce").dropna()
        sd = float(values.std(ddof=1)) if len(values) >= 2 else np.nan
        pooled[metric] = sd
        if not np.isfinite(sd) or sd == 0:
            skipped.append({"metric": metric, "reason": "pooled SD is zero or missing"})

    rows = []
    for metric in metric_order():
        sd = pooled.get(metric, np.nan)
        if not np.isfinite(sd) or sd == 0:
            continue
        for month in months:
            ck = keep_summary[
                keep_summary["metric_label"].eq(metric)
                & keep_summary["month"].eq(month)
                & keep_summary["treatment"].eq("CK")
            ]
            if ck.empty:
                continue
            ck_mean = float(ck.iloc[0]["mean"])
            for treatment in contrasts:
                hit = keep_summary[
                    keep_summary["metric_label"].eq(metric)
                    & keep_summary["month"].eq(month)
                    & keep_summary["treatment"].eq(treatment)
                ]
                if hit.empty:
                    continue
                mean = float(hit.iloc[0]["mean"])
                rows.append(
                    {
                        "row_label": f"{treatment}-CK {month}",
                        "metric": metric,
                        "month": month,
                        "treatment": treatment,
                        "mean_treatment": mean,
                        "mean_CK": ck_mean,
                        "pooled_sd": sd,
                        "standardized_difference": (mean - ck_mean) / sd,
                    }
                )

    long_df = pd.DataFrame(rows)
    row_order = [f"{t}-CK {m}" for t in contrasts for m in months]
    matrix = long_df.pivot(index="row_label", columns="metric", values="standardized_difference").reindex(row_order)[metric_order()]

    fig, ax = plt.subplots(figsize=(11.2, 5.6))
    cmap = LinearSegmentedColormap.from_list("effect", ["#2D6F95", "#F7F5F0", "#B85C2F"], N=256)
    limit = np.nanmax(np.abs(matrix.to_numpy()))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    im = ax.imshow(matrix.to_numpy(), cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    ax.set_title("Soil multi-metric responses relative to CK", loc="left", weight="bold")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if pd.notna(value) and abs(value) >= 1.5:
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", color="white", weight="bold", fontsize=7)

    start = 0
    for group, metrics in METRIC_GROUPS.items():
        width = len(metrics)
        ax.add_patch(Rectangle((start - 0.5, -1.28), width, 0.18, color=GROUP_COLORS[group], clip_on=False))
        ax.text(start + width / 2 - 0.5, -1.45, group, ha="center", va="bottom", fontsize=8, weight="bold")
        start += width

    ax.set_xticks(np.arange(-0.5, len(matrix.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("Standardized difference from CK")

    out = output_dir(config, "soil")
    save_figure(fig, out / "soil_multi_metric_responses_vs_ck")
    with pd.ExcelWriter(out / "soil_multi_metric_responses_vs_ck_source_data.xlsx") as writer:
        long_df.to_excel(writer, "LongEffects", index=False)
        matrix.to_excel(writer, "HeatmapMatrix")
        pd.DataFrame({"metric": list(pooled), "pooled_sd": list(pooled.values())}).to_excel(writer, "PooledSD", index=False)
        pd.DataFrame(skipped).to_excel(writer, "Skipped", index=False)
    print(f"Saved outputs to {out}")


if __name__ == "__main__":
    main()
