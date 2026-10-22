# SONFH functional transcriptomic analyses and six publication figures

Prepared 22 October 2026. Public repository: https://github.com/13016487621bdd/SONFH-Reproducibility . This repository provides a versioned Git commit; no DOI has been minted.

This archive supports a cross-sectional secondary analysis of GSE123568 (30 SONFH cases and 10 glucocorticoid-exposed participants without SONFH) and exploratory femoral head tissue analyses of GSE345158 (5 SONFH samples and 4 fracture controls).

## Repository layout

- `core_analysis/`: archived inputs, code, results, and stage-level execution evidence. Its original README and legacy plotting scripts are retained as historical provenance; generated legacy figure exports are omitted. Current manuscript illustration is under `publication_figures/`.
- `native_R/Figure1` through `Figure6`: self-contained R packages for the study diagram or full figure-specific recomputation. Open `RUN_FIGURE1.R` through `RUN_FIGURE6.R` in the respective directory and use Source in RStudio. Each package has its own inputs and instructions.
- `publication_figures/`: the current six-figure design, exact plot inputs, vector PDF/SVG outputs; raster exports are reproducible from the drawing code. `REDRAW_ALL_FIGURES.R` redraws the layout from checked plot data; it does not replace model recomputation.
- `provenance/GEO_20261022/`: official metadata and successful count-file download/checksum verification.
- `download_tissue_counts.R`: optional repeat of the public count-file download with MD5 verification.
- `SHA256SUMS.tsv`: hashes of all packaged files except this manifest itself. Duplicate historical image exports and large raster exports are omitted from this code repository; the independent manuscript submission figure package retains them.

## R workflows

The native environment recorded R 4.5.0, GSVA 2.4.9, clubSandwich 0.7.0, DESeq2 1.50.2, readxl 1.5.0.1, and SummarizedExperiment 1.40.0. Preserve these archived versions when checking numerical equality. Required packages are checked by scripts; they are not installed automatically.

For the archived integrated core workflow, set the R working directory to `core_analysis` and run:

```r
source("code/run_all_R.R", encoding = "UTF-8")
```

The convenience launcher was syntax-checked but not executed as a single fresh job. Native R execution evidence exists for its component stages. The original blood base-R upstream wrapper itself was not run; the equivalent deposited-probe-to-gene path was executed and independently checked in Python. Figure-specific R packages have their own execution scopes and status reports. Figure 3 full recomputation was successfully returned by the user; do not infer fresh complete execution of every package from the existence of a plot PDF. The current polished layout-only R script was prepared but not executed in native R in the preparation environment.

For an individual figure, open its standalone script, for example `native_R/Figure3/RUN_FIGURE3.R`, and select Source. Statistical recomputation and figure redraws are distinct operations. Inspect generated `FINAL_STATUS.txt`, numerical checks and session information.

## Scope and inference

The three primary, parallel testing families contain 85 individual gene sets, 21 functional groups, and 11 opposing-process contrasts, with separate BH corrections. They are not a sequential significance filter. Layer 2 compares the disease effects on standardized constituent-score differences, without separately declaring either cohort physiologically abnormal. Fixed-score leave-one-out analyses do not recalculate full-cohort scoring. Proxy analyses are exploratory expression adjustments rather than measured blood-cell fractions.

Blood starts from deposited RMA expression, not CEL-to-RMA reconstruction. Tissue starts from deposited integer counts, not FASTQ alignment/count generation. The fixed tissue reconstruction uses ratio size factors and DESeq2 `vst` with `design = ~1`, `blind = TRUE`, `nsub = 1000`, and `fitType = "parametric"`. Frozen annotation and the SOD2 correction are supplied inputs. No current annotation query or P-value-based recipe choice is required.

## Public tissue source verified on 22 October 2026

- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345158
- Official file directory: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345158/suppl/
- Counts used here: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345158/suppl/GSE345158_gene_count_exp.annotation.xlsx
- Official metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345158/soft/GSE345158_family.soft.gz

The count file was downloaded successfully from NCBI and matched the archived analysis input byte for byte (MD5 `5ae538bfa2cc1c5aee87f539e8be2768`; SHA256 `f565ea527d67e3ba98a2f1e98939cea939ef982dc26134e42bc2920745f86de4`). It contains 20,345 gene rows and nine sample columns. The TPM file is a different supplementary file and was not used for count-to-VST reconstruction.

The GEO title uses “healthy controls”, while its overall design explicitly identifies these as femoral head fracture patients and reports age/sex matching. This repository uses the more specific term “fracture controls”. Individual-level age/sex matching was not independently verified. These are not controls matched to the blood cohort for glucocorticoid exposure; tissue evidence is exploratory, not external clinical validation of a blood risk marker.

## Data and software attribution

GEO source accession GSE123568: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123568 ; platform GPL15207. Tissue platform GPL24676 is an Illumina NovaSeq 6000 platform, with BioProject PRJNA1518914. Source annotations and archived source records are preserved with their origins. Read `THIRD_PARTY_NOTICES.md` for MSigDB/GO attribution. No blanket license is asserted over third-party data or archived author-provided source scripts.

No manuscripts, author correspondence, credentials or account tokens are intentionally included. This release records supplied archived results and stage-specific verification; it does not claim prospective prediction, diagnostic utility or a completed journal submission.

## Distributed blood annotation extract

The original GPL15207 full annotation archive was 40,068,950 bytes. To keep the public repository practical, `core_analysis/inputs/01_Blood_Source/GPL15207_analysis_fields.zip` distributes exactly the three archived columns consumed by the upstream scripts: ID, Gene Symbol, and Annotation Date. It is an explicitly documented derivative, not a byte-identical copy of the complete source download. No annotation values were updated from a current database. The archived SHA256 identifies the original full archive; the root SHA256 manifest identifies the distributed files. See `provenance/GPL15207_analysis_fields.json` and the reproduced expression-matrix check. Historical input hash reports under core_analysis remain provenance of the original full input and are not a manifest of this revised distribution.

The distributed annotation extract was executed through the deposited-probe-to-gene reconstruction on 22 October 2026. All 18,833 gene identifiers, all 40 sample identifiers and their orders matched; 753,320 expression values agreed within 1e-12 (maximum absolute difference 1.7763568394002505e-15). This is a new upstream verification of the distribution change, not a new execution of all native R statistical analyses.

## Audited archived objects

The tissue analysis object contains only nine generic GEO sample labels and group categories, frozen gene annotation, gene-set members, size factors and statistical summaries. Legacy tissue effect/P values in that object are superseded as explained in the core README; they are not reconstruction targets for the corrected signal. The original master workbook records definitions, overlap checks and planning. Its historical angiogenesis–osteoblast cross-domain planning worksheet was not implemented in this manuscript and is not one of the 11 opposing-process contrasts. No result for that planned comparison is claimed. See `provenance/Public_disclosure_content_audit_20261022.json`.

## Figure 6 intermediate vector output

The redundant archived `Figure6_native_R_source.pdf` intermediate is not distributed. Full R recomputation and drawing require no such PDF. Run `native_R/Figure6/RUN_FIGURE6.R` to produce `results/<run>/figure/Figure6_R.pdf`. The optional Python layout exporter locates that locally generated native R output when redrawing Figure 6. The final six publication PDFs/SVGs are included. This is a distribution choice, not a change to the analysis or the plotted values.
