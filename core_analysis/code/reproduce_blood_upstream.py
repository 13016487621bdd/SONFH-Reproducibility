"""Rebuild gene expression from deposited RMA probe expression and GPL annotation.

The two historical audit tables are OUTPUTS here, not independent validation inputs.
No additional normalization, log transformation, or gene-symbol modernization is used.
Run from any directory; defaults resolve from this file's package root.
"""
from pathlib import Path
import argparse, gzip, hashlib, io, json, platform as sysplatform, sys, zipfile
import numpy as np
import pandas as pd

def unique_file(folder,pattern):
    matches=list(folder.glob(pattern))
    if len(matches)!=1:
        raise FileNotFoundError(f'Expected exactly one {pattern} in {folder}; found {len(matches)}. Supply an explicit command-line path.')
    return matches[0]

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--series',type=Path)
    ap.add_argument('--platform',type=Path)
    ap.add_argument('--reference',type=Path,default=root/'inputs/03_QC/SONFH_QC03_gene_expression_median_log2.csv.gz')
    ap.add_argument('--out',type=Path,default=root/'results/blood_upstream_python')
    a=ap.parse_args()
    a.series=a.series or unique_file(root/'inputs/01_Blood_Source','GSE123568_series_matrix*.gz')
    a.platform=a.platform or unique_file(root/'inputs/01_Blood_Source','GPL15207*.zip')
    for p in [a.series,a.platform,a.reference]:
        if not p.is_file():raise FileNotFoundError(p)

    metadata={}
    with gzip.open(a.series,'rt',encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            if line.startswith(('!Series_geo_accession','!Sample_data_processing','!Sample_geo_accession')):
                key=line.split('\t',1)[0];metadata.setdefault(key,[]).append(line.rstrip())
        else:raise ValueError('Series matrix start marker missing')
        lines=[];found_end=False
        for line in f:
            if line.startswith('!series_matrix_table_end'):
                found_end=True;break
            lines.append(line)
        assert found_end,'Series matrix end marker missing'
    probes=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0)
    assert probes.shape==(49293,40),probes.shape
    assert probes.index.is_unique and probes.columns.is_unique
    assert np.isfinite(probes.to_numpy()).all()
    with zipfile.ZipFile(a.platform) as z:
        candidates=[n for n in z.namelist() if not n.endswith('/') and Path(n).name.startswith('GPL15207') and n.lower().endswith('.txt')]
        assert len(candidates)==1,candidates
        member=candidates[0]
        with z.open(member) as f:
            gpl=pd.read_csv(f,sep='\t',comment='#',usecols=['ID','Gene Symbol','Annotation Date'],dtype=str)
    assert gpl.ID.is_unique,'Duplicate platform probe IDs'
    sym=gpl.set_index('ID')['Gene Symbol'].reindex(probes.index).fillna('').str.strip()
    missing=sym.isin(['','---'])
    ambiguous=~missing & sym.str.contains('///',regex=False)
    unique=~missing & ~ambiguous
    nonstandard=unique & sym.str.contains('/',regex=False)
    retained=unique & ~nonstandard
    assert int(missing.sum())==430
    assert int(ambiguous.sum())==1994
    assert int(unique.sum())==46869
    assert set(sym[nonstandard].index)=={'11761930_at','11762411_at'}
    assert int(nonstandard.sum())==2 and int(retained.sum())==46867
    rebuilt=probes.loc[retained].copy()
    rebuilt.index=sym[retained].to_numpy()
    rebuilt=rebuilt.groupby(level=0,sort=True).median()
    frozen=pd.read_csv(a.reference,index_col=0)
    assert rebuilt.shape==frozen.shape==(18833,40)
    assert rebuilt.index.equals(frozen.index),'Gene identifiers/order differ'
    assert rebuilt.columns.equals(frozen.columns),'Sample identifiers/order differ'
    diff=np.abs(rebuilt.to_numpy()-frozen.to_numpy())
    maxdiff=float(diff.max())
    assert np.isfinite(diff).all() and maxdiff<1e-12,maxdiff
    a.out.mkdir(parents=True,exist_ok=True)
    qc2=pd.DataFrame({'ID_REF':probes.index[unique],'Gene_Symbol':sym[unique].to_numpy()})
    qc2.to_csv(a.out/'QC02_probe_mapping_RECONSTRUCTED.csv',index=False)
    decisions=np.select([missing,ambiguous,nonstandard],['EXCLUDE_UNANNOTATED','EXCLUDE_AMBIGUOUS_MULTIGENE','EXCLUDE_NONSTANDARD_SLASH_SYMBOL'],default='RETAIN_UNIQUE_GENE')
    pd.DataFrame({'ID_REF':probes.index,'Gene_Symbol':sym.to_numpy(),'decision':decisions}).to_csv(a.out/'QC03_probe_decisions_RECONSTRUCTED.csv.gz',index=False)
    pd.DataFrame({'ID_REF':sym[nonstandard].index,'Gene_Symbol':sym[nonstandard].to_numpy()}).to_csv(a.out/'nonstandard_excluded_probes.csv',index=False)
    # The frozen reference is never overwritten; retain a separate regenerated matrix.
    rebuilt.to_csv(a.out/'QC03_gene_expression_REPRODUCED.csv.gz',float_format='%.16g')
    source_rows=[]
    for role,p in [('deposited_probe_expression',a.series),('platform_annotation',a.platform),('frozen_reference',a.reference)]:
        source_rows.append({'role':role,'filename':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    pd.DataFrame(source_rows).to_csv(a.out/'source_SHA256.csv',index=False)
    summary=dict(status='PASS',scope='Deposited RMA-normalized probe matrix to frozen gene-level matrix; not CEL-to-RMA reproduction',
        input_probes=49293,samples=40,platform_rows=int(len(gpl)),annotation_dates=sorted(gpl['Annotation Date'].dropna().unique().tolist()),
        unannotated_probes=int(missing.sum()),ambiguous_probes=int(ambiguous.sum()),unambiguous_probes_before_symbol_filter=int(unique.sum()),
        nonstandard_excluded=int(nonstandard.sum()),retained_probes=int(retained.sum()),final_genes=18833,
        expression_values_compared=int(diff.size),gene_identifiers_and_order='EXACT',sample_identifiers_and_order='EXACT',
        maximum_absolute_difference=maxdiff,numeric_tolerance=1e-12,
        audit_provenance='QC02 mapping and QC03 decisions newly reconstructed from the two sources; not compared to unavailable historical audit tables',
        reference_role='Validation only; not used to select probes, construct annotations, or calculate gene medians',
        python=sys.version,numpy=np.__version__,pandas=pd.__version__,system=sysplatform.platform())
    (a.out/'blood_upstream_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    pd.DataFrame({'sample':rebuilt.columns,'max_absolute_difference':diff.max(axis=0)}).to_csv(a.out/'samplewise_matrix_comparison.csv',index=False)
    (a.out/'source_preprocessing_metadata.txt').write_text('\n'.join(v for values in metadata.values() for v in values)+'\n',encoding='utf-8')
    report='BLOOD UPSTREAM REPRODUCTION PASSED\n\n'+'\n'.join(f'{k}: {v}' for k,v in summary.items())+'\n'
    (a.out/'blood_upstream_report.txt').write_text(report,encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
