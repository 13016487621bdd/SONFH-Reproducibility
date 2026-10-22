# SONFH E1 consistency closure audit

Date: 2026-10-18

## Scope
Recomputed the E1 erythroid/globin-associated expression proxy and its exploratory sensitivity models from the archived 40-sample gene-level expression matrix, sample manifest, and frozen GSVA score matrix. This audit is limited to E1 consistency and does not replace the prespecified primary analyses.

## Frozen E1 proxy definition
- Requested marker panel in the manuscript/table: HBB, ALAS2, AHSP, GYPA, SLC4A1.
- Measurable genes actually used in the frozen E1 proxy: ALAS2, AHSP, GYPA, SLC4A1 (four genes).
- HBB is present in the final expression matrix but is not part of the frozen four-gene E1 proxy. It is handled separately in E4 as a sensitivity definition.
- For each gene, expression was standardized across the 40 participants using the sample standard deviation; the proxy is the arithmetic mean of the four standardized gene values.

## Independent recomputation results
- N = 40 (30 SONFH; 10 glucocorticoid-exposed non-SONFH controls).
- Four-gene proxy mean: controls = 1.1305; SONFH = -0.3768.
- Spearman correlation between proxy and hydrogen peroxide catabolism GSVA score: rho = 0.7694, P = 6.58e-09.
- Primary sex-adjusted catabolism model: beta = -0.3855, 95% CI [-0.5180, -0.2530], P = 8.76e-07.
- Exploratory model additionally adjusted for E1 proxy: beta = -0.0978, 95% CI [-0.3335, 0.1379], P = 0.4055.
- Primary biosynthesis-minus-catabolism contrast: beta = 1.5203, P = 0.00358.
- Exploratory proxy-adjusted contrast: beta = 0.8887, 95% CI [-0.4714, 2.2488], P = 0.1935.
- Recomputed primary FDR results across 85 gene sets: 8 significant, matching the frozen primary result file.
- Recomputed exploratory proxy-adjusted models across 85 gene sets: 0 significant at BH-FDR < 0.05.
- The recalculated all-85 exploratory results match the prior E3 output to floating-point precision (maximum differences < 2e-15 for reported p-values and < 5e-16 for coefficients).
- HBB-inclusive five-gene definition is retained as a separate sensitivity analysis in E4; it does not replace the frozen four-gene E1 proxy.

## Interpretation boundary
The proxy is an unvalidated erythroid/globin-associated expression score, not a measured cell fraction. Attenuation after adjustment does not prove cell-composition confounding: the proxy may reflect cell abundance, cell state, disease-related expression, treatment, or specimen-related effects, and adjustment may remove disease-related signal. These exploratory models must not replace the prespecified primary models or be interpreted causally.

## Manuscript alignment
The current E1-harmonized manuscript states the exact four-gene proxy, reports the association and attenuation results, explicitly labels proxy-adjusted models exploratory, and treats HBB inclusion separately. Sections to keep aligned are Abstract, Results 3.5, Methods 2.7, Discussion, and Supplementary Tables E1-E4.

## Remaining checks
1. Ensure the manuscript's Supplementary Information actually includes the included E1-E4 tables and the gene list used for each proxy.
2. Keep E1 column names distinguishing `genes_requested` from `genes_measured_and_used`.
3. Do not describe this proxy as a cell fraction or assert that it proves cell-composition confounding.
4. Keep primary 85-set findings unchanged; proxy-adjusted models are exploratory only.
