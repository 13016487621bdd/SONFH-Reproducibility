#!/usr/bin/env python3
"""QC03 upstream reproducibility audit; does not modify frozen QC03."""
from pathlib import Path
import gzip
import zipfile
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
series=ROOT/'GSE123568_series_matrix.txt.gz'
platform=ROOT/'GPL15207_annotation.zip'
reference=ROOT.parent/'SONFH_QC03_gene_expression_median_log2.csv.gz'
map_audit=ROOT/'SONFH_QC02_probe_to_gene_unambiguous.csv'
decisions=ROOT/'SONFH_QC03_probe_decision_audit.csv.gz'

# Fail early with a complete list of missing required input files.
required = {
    'GSE123568 Series Matrix': series,
    'GPL15207 platform annotation': platform,
    'Frozen QC03 reference matrix': reference,
    'Archived QC02 probe mapping': map_audit,
    'Archived QC03 probe decisions': decisions,
}
missing_inputs = {name: path for name, path in required.items() if not path.is_file()}
if missing_inputs:
    details = "\n".join(f"  - {name}: {path}" for name, path in missing_inputs.items())
    raise FileNotFoundError(
        "QC03 reproduction cannot start; required input files are missing:\n"
        + details
        + "\nSee README.md for the expected folder layout."
    )

with zipfile.ZipFile(platform) as z:
    with z.open(z.namelist()[0]) as fh:
        gpl=pd.read_csv(fh,sep='\t',skiprows=36,usecols=['ID','Gene Symbol'],dtype=str)
with gzip.open(series,'rt') as fh:
    for line in fh:
        if line.startswith('!series_matrix_table_begin'):
            break
    probe_matrix=pd.read_csv(fh,sep='\t',index_col=0)
probe_matrix=probe_matrix.loc[~probe_matrix.index.str.startswith('!series_matrix_table_end')]
probe_matrix=probe_matrix.apply(pd.to_numeric)

assert gpl.ID.is_unique and probe_matrix.index.is_unique
gpl_index=gpl.set_index('ID')['Gene Symbol'].reindex(probe_matrix.index)
text=gpl_index.fillna('').str.strip()
missing=text.eq('') | text.eq('---')
ambiguous=text.str.contains('///',regex=False)
unambig=~(missing|ambiguous)
raw_unique=pd.DataFrame({'ID_REF':probe_matrix.index[unambig], 'Gene_Symbol':text[unambig].to_numpy()})
legacy=pd.read_csv(map_audit,dtype=str)
assert len(raw_unique)==46869, (len(raw_unique),'unexpected unambiguous count')
assert raw_unique.equals(legacy), 'Platform mapping differs from archived QC02 mapping'

nonstandard=raw_unique.Gene_Symbol.str.contains('/',regex=False)
exceptional=raw_unique[nonstandard]
assert set(exceptional.ID_REF)=={'11761930_at','11762411_at'}
retained=raw_unique[~nonstandard].copy()
assert len(retained)==46867 and retained.Gene_Symbol.nunique()==18833

qc_decisions=pd.read_csv(decisions,dtype=str)
actual_kept=qc_decisions.loc[qc_decisions.decision=='RETAIN_UNIQUE_GENE',['ID_REF','Gene_Symbol']]
assert retained.reset_index(drop=True).equals(actual_kept.reset_index(drop=True)), 'Independent annotation filtering disagrees with archived decisions'

matrix=probe_matrix.loc[retained.ID_REF].copy()
matrix.index=retained.Gene_Symbol.to_numpy()
reproduced=matrix.groupby(level=0,sort=True).median()
frozen=pd.read_csv(reference,index_col=0)
assert reproduced.shape==(18833,40) and frozen.shape==(18833,40)
assert reproduced.index.equals(frozen.index) and reproduced.columns.equals(frozen.columns)
diff=np.abs(reproduced.to_numpy()-frozen.to_numpy())
assert np.isfinite(diff).all()
max_abs=float(diff.max())
assert max_abs<1e-12

report=f'''SONFH QC03 reproduction audit
==============================
Series Matrix probes: {len(probe_matrix):,}
Samples: {probe_matrix.shape[1]}
GPL annotation probe rows: {len(gpl):,}
Unmapped/blank: {int(missing.sum()):,}
Ambiguous (///): {int(ambiguous.sum()):,}
Unambiguous mapping: {len(raw_unique):,}
Initial unique gene symbols: {raw_unique.Gene_Symbol.nunique():,}
Nonstandard excluded probes: {exceptional.ID_REF.to_list()}
Their symbols: {exceptional.Gene_Symbol.to_list()}
Final retained probes: {len(retained):,}
Reproduced genes: {reproduced.shape[0]:,}
Samples: {reproduced.shape[1]:,}
Compared expression cells: {diff.size:,}
Gene names & order: EXACT
Sample names & order: EXACT
Maximum absolute numeric difference: {max_abs:.17g}
All numerical differences <= 1e-12: {bool((diff<=1e-12).all())}
Note: float64 conversion/CSV parsing explains machine-precision differences.
Source GPL annotation date: 30-Mar-2016, as present in GPL metadata.
''' 
print(report)
(ROOT/'SONFH_QC03_reproduction_audit_report.txt').write_text(report,encoding='utf-8')
reproduced.to_csv(ROOT/'SONFH_QC03_reproduced_temporary.csv.gz',compression='gzip',float_format='%.16g')
