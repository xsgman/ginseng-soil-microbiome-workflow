# Data Layout

This repository does not include private experimental data.

Prepare your own files under `data/processed/` or another local path, then copy
`config.example.toml` to `config.toml` and update the paths.

Expected inputs:

- Soil workbook with sheets `Summary` and `ReplicateMeans`.
- Ginseng workbook with sheet `CK_KXH_RHB_Summary`; optional sheet `StatsFocus`.
- Metagenomics files:
  - `pcoa_scores.csv`
  - optional `pcoa_group_centers.csv`
  - `function_log2FC_vs_CK.csv`
  - `genus_abundance.csv`

`genus_abundance.csv` can be either:

1. rows as genera and columns as samples, with the first column named `genus`; or
2. columns as genera and rows as samples, with one sample ID column.

Do not commit raw sequencing files, Excel workbooks, private results, or generated
figures to this repository.
