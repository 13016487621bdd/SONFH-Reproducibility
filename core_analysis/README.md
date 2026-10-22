# SONFH JTM figures and reproducibility package

Updated 21 October 2026. This version includes both blood source files and the native R evidence returned by the user. Original input matrices and legacy scripts are preserved under `inputs/`; current entry points are under `code/`.

## Verified outcome

Reproduction is verified stage by stage from the supplied inputs to the reported analyses and figure data. The evidence spans Python in the preparation environment and native R on the user's Windows computer; a single fresh run of the integrated launcher was not performed.

| Stage | Executed evidence | Outcome |
| --- | --- | --- |
| Blood deposited RMA probe matrix + GPL15207 to gene matrix | `results/blood_upstream_python/` | 49,293 probes to 18,833 genes x 40 samples; identifiers/order exact; maximum error 1.78e-15 |
| Frozen scores to principal statistics, sensitivity analyses and figure data | `results/verification_checks.csv` | All 76 Python numerical checks passed |
| Native GSVA/ssGSEA, CR2 inference, gene-removal scoring and tissue scoring sensitivity | `results/native_R/` | All 16 checks passed in user-run R |
| Tissue integer counts to VST, frozen annotation, GSVA and corrected Pair11 | `results/tissue_native_R/` | VST and mapped matrices exactly matched; Z/contrast error 4.22e-15; effect 1.2937192585681867; exact P=5/126 |
| Figures, tables, input checksums | `results/artifact_validation.csv`, `results/input_SHA256.csv` | Five main figures, three supplementary figures, ten supplementary tables |

The principal findings remain 8/85 individual gene-set associations, 0/21 functional-group comparisons and 1/11 opposing-process contrasts after their separate BH corrections. The tissue comparison uses five SONFH and four fracture-control samples.

## Scope

Blood starts at the deposited RMA-normalized probe expression matrix; CEL-to-RMA processing is outside this package. Tissue starts at supplied integer gene counts; FASTQ alignment and count generation are outside this package. The frozen 85-set panel and frozen annotation tables are analysis inputs. The package does not depend on querying a current annotation database.

The original tissue preprocessing source was incomplete. Four fixed transformation recipes were compared with the archived full expression matrix, without using a P-value criterion. Exactly one matched: an intercept-only design with DESeq2 `vst`, `blind=TRUE`, `nsub=1000`, `fitType="parametric"`, and ratio size factors. This establishes an equivalent reconstructed workflow, not uniquely recoverable historical source settings. The default tissue entry point now fixes that matched recipe.

## Current entry points

From the extracted `SONFH_JTM` project directory:

```bash
python code/reproduce_blood_upstream.py
python code/reproduce.py
python code/make_figures.py
python code/package_checks.py
```

For native R, set the working directory to `SONFH_JTM` and run:

```r
source("code/reproduce_tissue_upstream.R", encoding="UTF-8")
source("code/verify_scoring_and_CR2.R", encoding="UTF-8")
```

`code/reproduce_blood_upstream.R` provides a base-R alternative to the executed Python probe-to-gene path. It has not itself been executed in R. `code/run_all_R.R` is a convenience launcher for the three native R entries; it was syntax-checked, not executed as one integrated job. The user does not need to rerun these stages to validate the evidence already returned.

The canonical tissue entry point was derived from the user-executed recovery script by selecting its matching recipe and adapting project paths. Its unchanged numerical operations have native execution evidence; this integrated entry point was syntax-checked only. The actual recovery script is preserved as `results/tissue_native_R/RUN_TISSUE_CHECK_executed.R`. `code/verify_tissue_upstream.R` is the original standalone four-candidate recovery script, requiring the standalone tissue-check layout.

Use Python 3.12 and the dependencies in `requirements.txt` or `requirements-tested.txt`. Native execution used R 4.5.0, GSVA 2.4.9, clubSandwich 0.7.0, DESeq2 1.50.2, readxl 1.5.0.1 and SummarizedExperiment 1.40.0. Full returned environment information is retained with the evidence.

## Provenance and corrected definitions

- Blood QC02/QC03 audit tables were regenerated from the uploaded series matrix and platform annotation. They are newly reconstructed records, not claims of equality with unavailable historical audit tables. The frozen reference expression matrix was used only for final validation.
- The erythroid expression proxy uses ALAS2, AHSP, GYPA and SLC4A1. HBB is added only in a separate five-gene sensitivity analysis. Supplementary tables use that four-gene primary definition consistently.
- Tissue annotation removes rows with missing symbols, then simultaneously excludes all rows whose ENSEMBL identifier or symbol is duplicated in the valid table. This reproduces 18,572 unique mappings. ENSG00000112096 is appended as SOD2, producing 18,573 genes.
- Tissue figures use SOD2-corrected scores. The older uncorrected effect 1.4678777 and P=0.015873 are superseded.
- The old shared-gene sensitivity table reports the corrected baseline effect approximately 1.293720, 7.4e-7 above the sample-level value. New Table S10 uses the verified sample-level effect and retains the archived value in a separate provenance column. Shared-gene exclusion rescoring was independently verified in user-run R.
- User-run R evidence is identified as such. The pasted tissue report is saved as a normalized excerpt; its full intermediate reconstructed RDS remains on the user's computer. All four input MD5 values from that report match the packaged inputs.

## Submission files

- `figures/Figure1.pdf` through `Figure5.pdf`: main figures, with SVG and 600-dpi LZW TIFF alternatives.
- `figures/Supplementary_Figure_S1.pdf` through `S3.pdf`: supplementary figures; S1 has three pages.
- `SONFH_JTM_Figure_Review.pdf`: ten-page visual review file.
- `documentation/Figure_legends.txt`: English figure titles and legends.
- `documentation/Methods_Results_insertions.txt`: manuscript additions, including verified preprocessing settings and figure/table citation map.
- `documentation/Submission_audit_CN.txt`: current Chinese audit summary and manuscript assembly items.
- `results/Supplementary_Table_S1_*.csv` through `S10_*.csv`: regenerated supplementary tables.
- `results/Table1_available_sample_characteristics.csv`: available sample characteristics.
- `results/execution_evidence.json`: stage-level execution provenance and scope.

Main figures are 170 mm wide; vector PDF/SVG and lossless 600-dpi TIFF are included. The figure preparation reference was the JTM submission guideline accessed 20 October 2026: https://link.springer.com/journal/12967/submission-guidelines

This package completes the scoped computational checks and figure preparation. Final manuscript assembly still requires inserting the updated methods/results and figure citations, resolving Introduction placeholders `[REF-A]` and `[REF-B]`, and completing references, author information, declarations and the cover letter. It does not establish that journal submission has occurred.
