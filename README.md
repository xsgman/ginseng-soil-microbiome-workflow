# Ginseng Soil-Microbiome Figure Workflow

Reusable Python workflow for preparing presentation-ready figures from ginseng
soil, phenotype, disease, and metagenomics summary tables.

This repository is a workflow template. It does not contain private experimental
data, generated figures, local disk paths, raw sequencing files, or Excel
workbooks.

## Figures

The workflow covers five figure types:

1. Soil multi-metric responses relative to CK.
2. Integrated soil environmental PCA.
3. Ginseng phenotype and disease response effect plot.
4. Metagenomics PCoA plus functional response heatmap.
5. Bacillus-centered candidate validation network.

## Quick Start

```bash
python -m venv .venv
.venv\\Scripts\\activate
pip install -r requirements.txt
copy config.example.toml config.toml
```

Edit `config.toml` so that it points to your own processed data files.

Run scripts individually:

```bash
python scripts/plot_soil_heatmap.py --config config.toml
python scripts/plot_soil_pca.py --config config.toml
python scripts/plot_ginseng_effect.py --config config.toml
python scripts/plot_metagenomics_panel.py --config config.toml
python scripts/plot_bacillus_network.py --config config.toml
```

Outputs are written to `outputs/figures/` by default.

## Data Privacy

The `.gitignore` excludes common private and large file types, including:

- Excel workbooks
- CSV/TSV/TXT data tables
- FASTQ/GZ/BAM/SAM sequencing files
- ZIP/RAR archives
- generated figures and output folders

Before publishing, run a manual check for private paths and sensitive files.

## Citation / Use

This is an analysis and plotting workflow template for internal research and
presentation use. Validate all statistics and figure interpretations against
your own experimental design before publication.
