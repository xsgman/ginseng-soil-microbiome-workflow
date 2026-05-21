# Figure Explanation

## 1. Soil Multi-Metric Responses Relative to CK

Purpose: show how KXH and RHB shift soil indicators relative to the matched CK
baseline.

Calculation:

```text
standardized difference = (treatment mean - matched CK mean) / pooled SD
```

`pooled SD` is the standard deviation of the same metric across CK, KXH, RHB
and selected months. It puts pH, EC, carbon, nutrients, and enzyme activities
on a common scale.

Interpretation: red/blue indicate increase/decrease relative to CK. They do not
automatically mean good/bad; interpretation depends on the metric.

## 2. Integrated Soil PCA

Purpose: test whether all soil metrics together separate CK, KXH, and RHB.

Each sample is represented by many soil metrics. Metrics are z-score
standardized and compressed into PC1 and PC2. Points close together have similar
integrated soil states; points far apart differ more strongly.

## 3. Ginseng Phenotype Effect Plot

Purpose: summarize phenotype and disease responses relative to CK on one
improvement-oriented axis.

The plotted value is `benefit_log2_vs_ck`. Positive values mean the direction
is favorable after harmonizing metric direction: lower disease index is better,
whereas higher biomass and root diameter are better.

## 4. Metagenomics PCoA + Functional Heatmap

Purpose: connect treatment effects with root-zone microbiome structure and
functional potentials.

Panel A: Bray-Curtis PCoA of metagenomic community composition.  
Panel B: functional log2 fold change relative to matched CK.

Functional potentials are gene-level potential, not direct proof of expression
or activity.

## 5. Bacillus-Centered Validation Network

Purpose: prioritize qPCR or plate-count validation targets, not to present a
full co-occurrence network.

Edges are Spearman correlations. The default thresholds are:

```text
P < 0.05 and |rho| >= 0.35
```

Correlation does not prove inhibition, mechanism, or ecological causality.
