from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from workflow_utils import cfg_path, load_config, normalize_soil_metric, output_dir, parse_args, save_figure, setup_style, zscore


COLORS = {"CK": "#59616C", "KXH": "#2F86A6", "RHB": "#B86B2B"}
MARKERS = {"Jun": "o", "Jul": "s", "Aug": "^", "Sep": "D"}
LOADINGS_KEEP = ["pH", "Al3+", "H+", "SOC", "AN", "NH4-N", "Total acidity", "CEC"]


def pca_svd(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    x = matrix.to_numpy(dtype=float)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    loadings = vt[:2, :].T
    explained = (s**2) / np.sum(s**2)
    scores_df = pd.DataFrame(scores, columns=["PC1", "PC2"], index=matrix.index)
    loadings_df = pd.DataFrame(loadings, columns=["PC1", "PC2"], index=matrix.columns).reset_index(names="metric")
    explained_df = pd.DataFrame({"PC": ["PC1", "PC2"], "explained_ratio": explained[:2]})
    return scores_df, loadings_df, explained_df


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    setup_style(8.5)

    months = config["plot"].get("months_soil", ["Jun", "Jul", "Aug", "Sep"])
    treatments = config["plot"].get("treatments_focus", ["CK", "KXH", "RHB"])
    workbook = cfg_path(config, "soil_workbook")
    reps = pd.read_excel(workbook, sheet_name="ReplicateMeans")
    reps["metric_label"] = reps["metric_en"].map(normalize_soil_metric)
    keep = reps[reps["phase"].eq("treatment") & reps["month"].isin(months) & reps["treatment"].isin(treatments)].copy()
    keep["sample_id"] = keep["treatment"].astype(str) + "_" + keep["month"].astype(str) + "_R" + keep["replicate"].astype(str)

    wide = keep.pivot_table(
        index=["sample_id", "treatment", "month", "replicate"],
        columns="metric_label",
        values="value",
        aggfunc="mean",
    ).reset_index()

    metric_cols = [c for c in wide.columns if c not in {"sample_id", "treatment", "month", "replicate"}]
    imputed = wide.copy()
    for metric in metric_cols:
        imputed[metric] = pd.to_numeric(imputed[metric], errors="coerce")
        imputed[metric] = imputed.groupby(["treatment", "month"])[metric].transform(lambda s: s.fillna(s.mean()))
        imputed[metric] = imputed.groupby(["month"])[metric].transform(lambda s: s.fillna(s.mean()))
        imputed[metric] = imputed[metric].fillna(imputed[metric].mean())

    z = zscore(imputed[metric_cols]).dropna(axis=1, how="any")
    scores, loadings, explained = pca_svd(z)
    scores = pd.concat([imputed[["sample_id", "treatment", "month", "replicate"]], scores.reset_index(drop=True)], axis=1)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for treatment in treatments:
        sub_t = scores[scores["treatment"].eq(treatment)]
        for month, marker in MARKERS.items():
            sub = sub_t[sub_t["month"].eq(month)]
            if sub.empty:
                continue
            ax.scatter(sub["PC1"], sub["PC2"], s=42, color=COLORS.get(treatment, "#888888"), marker=marker, alpha=0.9)

    top_loadings = loadings[loadings["metric"].isin(LOADINGS_KEEP)].copy()
    top_loadings["strength"] = np.sqrt(top_loadings["PC1"] ** 2 + top_loadings["PC2"] ** 2)
    top_loadings = top_loadings.sort_values("strength", ascending=False).head(6)
    scale = min(scores["PC1"].abs().max(), scores["PC2"].abs().max()) * 1.25
    for _, row in top_loadings.iterrows():
        x = row["PC1"] * scale
        y = row["PC2"] * scale
        ax.arrow(0, 0, x, y, color="#4B5159", width=0.006, head_width=0.08, length_includes_head=True)
        ax.text(x * 1.08, y * 1.08, row["metric"], fontsize=8, color="#333333")

    ax.axhline(0, color="#D9DEE3", lw=0.8)
    ax.axvline(0, color="#D9DEE3", lw=0.8)
    ax.grid(color="#EDF0F2", lw=0.6)
    ax.set_title("Integrated soil environmental states", loc="left", weight="bold")
    ax.set_xlabel(f"PC1 ({explained.loc[0, 'explained_ratio'] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained.loc[1, 'explained_ratio'] * 100:.1f}%)")

    treatment_handles = [
        Line2D([], [], marker="o", color="none", markerfacecolor=COLORS[t], markeredgecolor="none", label=t)
        for t in treatments
    ]
    month_handles = [
        Line2D([], [], marker=m, color="#666666", linestyle="None", label=month) for month, m in MARKERS.items() if month in months
    ]
    leg1 = ax.legend(handles=treatment_handles, title="Treatment", loc="upper right")
    ax.add_artist(leg1)
    ax.legend(handles=month_handles, title="Month", loc="lower right")

    out = output_dir(config, "soil")
    save_figure(fig, out / "soil_integrated_pca")
    with pd.ExcelWriter(out / "soil_integrated_pca_source_data.xlsx") as writer:
        wide.to_excel(writer, "WideRaw", index=False)
        imputed.to_excel(writer, "WideImputed", index=False)
        z.to_excel(writer, "ZscoreMatrix")
        scores.to_excel(writer, "PCAScores", index=False)
        loadings.to_excel(writer, "Loadings", index=False)
        explained.to_excel(writer, "ExplainedVariance", index=False)
    print(f"Saved outputs to {out}")


if __name__ == "__main__":
    main()
