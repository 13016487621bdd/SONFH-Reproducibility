# SONFH QC03 upstream reproduction

This folder provides a portable verification of the frozen QC03 gene-expression matrix. The computation and audit assertions are preserved from the supplied, previously successful script; only input paths and an upfront missing-file check have been updated.

## Folder layout

Place this directory inside the SONFH project's `03_QC` folder:

```text
03_QC/
|-- SONFH_QC03_gene_expression_median_log2.csv.gz   # Frozen reference: DO NOT OVERWRITE
`-- QC03_reproduction/
    |-- GSE123568_series_matrix.txt.gz
    |-- GPL15207_annotation.zip
    |-- SONFH_QC02_probe_to_gene_unambiguous.csv
    |-- SONFH_QC03_probe_decision_audit.csv.gz
    |-- SONFH_QC03_reproduction_audit.py
    `-- README.md
```

Important: the four input files inside `QC03_reproduction/` and the frozen reference matrix at `03_QC/` must be copied into place. They are **not included** in this script-only archive. File names must exactly match the layout above. Preserve the original file bytes when renaming.

## Python environment

Requires Python 3 and packages `numpy` and `pandas` (also Python standard-library `gzip`, `pathlib`, and `zipfile`). Example install:

```text
python -m pip install numpy pandas
```

## Run on Windows (PowerShell / command prompt)

Open a terminal in `03_QC/QC03_reproduction/`, then run:

```text
python SONFH_QC03_reproduction_audit.py
```

Paths are resolved relative to the script, not the current working directory, so it also works when launched using an absolute script path.

## Expected audited findings (from the prior successful run)

- Source expression matrix: 49,293 probes and 40 samples.
- Missing/unannotated probe mappings: 430.
- Ambiguous mappings (`///`): 1,994.
- Two additional slash-containing nonstandard gene symbols removed.
- Retained probes: 46,867; final genes: 18,833.
- Frozen and reproduced matrices: exactly matching gene and sample names/order.
- Elementwise comparison: 753,320 cells; maximum prior observed difference 1.78e-15, below the script's 1e-12 tolerance.

The script checks archived intermediate mapping and decision files before comparing the reconstructed matrix with the frozen reference. It raises an error if an assertion fails.

## Outputs (created only after passing checks)

- `SONFH_QC03_reproduction_audit_report.txt`: audit report.
- `SONFH_QC03_reproduced_temporary.csv.gz`: reconstructed **temporary** gene-expression matrix.

Neither output overwrites the frozen reference matrix.

## Verification status

The original script reportedly passed on the user's prior dataset. This repackaged version was syntax-checked and its missing-input handling was tested, but **has not been rerun end-to-end here**, since the five required input files are not attached in this conversation. After putting those inputs in the folder, rerun once and archive the new report. Do not describe this repackaged script as independently rerun until that check succeeds.
