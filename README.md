# Metabolic Model Project

This repository contains a workflow for building a gene-expression-constrained metabolic model from the Human-GEM reconstruction. The project is implemented in both Python and Julia and is designed to map transcriptomics data onto reaction constraints before running FBA or FVA.

## Overview

The general workflow is:

1. Load the Human-GEM metabolic model.
2. Filter genes that are present in the model.
3. Normalise expression values.
4. Map gene expression to reaction-level constraints using GPR rules.
5. Adjust reaction bounds according to expression values.
6. Run constrained FBA or FVA to assess feasible metabolic fluxes.

## Repository structure

- `Julia/` — Julia implementation
  - `FVA_julia_notebook_.ipynb` — notebook for flux variability analysis in Julia
  - `map_gene_reaction.jl` — gene-to-reaction mapping and bound constraint logic
  - `data/normalised_cpm_counts.csv` — normalised gene expression data
  - `model/Human-GEM.mat` — Human-GEM model file
- `python/` — Python implementation
  - `Constrained_FBA.ipynb` — constrained FBA notebook
  - `Constrained_FVA.ipynb` — constrained FVA notebook
  - `util.py` — utility functions for reaction-bound adjustment
  - `data/normalised_cpm_counts.csv` — normalised gene expression data
  - `model/Human-GEM.mat` — Human-GEM model file
  - `results/` — output flux files

## Python workflow

The Python workflow uses COBRApy with the Human-GEM model and expression matrix. The utilities in `python/util.py` map gene-level expression onto reaction constraints and update lower and upper bounds before running the analysis.

Typical execution steps:

1. Open the notebook in VS Code or Jupyter.
2. Load the model and expression data.
3. Run the preprocessing and constraint steps.
4. Execute constrained FBA or FVA.
5. Inspect outputs saved in `python/results/`.

Expected output files include:

- `python/results/FBA_Fluxes.csv`
- `python/results/maxfluxes.csv`
- `python/results/minfluxes.csv`

## Julia workflow

The Julia workflow uses COBREXA-style reaction constraint logic with a gene-expression-to-reaction mapping function. It follows the same conceptual approach as the Python implementation but in Julia syntax and data structures.

Typical execution steps:

1. Open the Julia notebook or script.
2. Load the Human-GEM model.
3. Map each reaction to gene-expression-derived values using the GPR rules.
4. Apply reaction-bound adjustments.
5. Run FVA or relevant flux analysis.

## Requirements

### Python

Install the main dependencies with:

```bash
pip install pandas numpy scikit-learn tqdm cobra
```

A COBRA-compatible solver such as GLPK must also be available.

### Julia

The Julia workflow depends on packages used to read the model and work with COBREXA-style objects, including support for MATLAB `.mat` models. Use the Julia package manager to install the required libraries before running the notebook.

## Notes

- Gene-protein-reaction (GPR) rules are used to translate gene expression into reaction-level effects.
- Model outputs depend on the objective function, the selected solver, and the scaling parameter applied during bound adjustment.
