"""Executed structural audit; does not independently estimate DESeq2 dispersions."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    paths = [args.inputs/n for n in ('GSE345158_gene_count_exp.annotation.xlsx',
             'GSE345158_Pair11_analysis_objects.rds',
             'GSE345158_reproducibility_inputs.rds')]
    raw = pd.read_excel(paths[0], engine='openpyxl')
    archive = rdata.read_rds(paths[1]); ref = rdata.read_rds(paths[2])
    x = raw.iloc[:, 1:].to_numpy(float)
    expected = [f'Case-{i}' for i in range(1,6)] + ['Ctrl-1','Ctrl-2','Ctrl-4','Ctrl-5']
    assert raw.shape == (20345,10) and list(raw.columns[1:]) == expected
    assert np.isfinite(x).all() and (x>=0).all() and (x==np.round(x)).all()
    vst = ref['vst_matrix']
    assert np.array_equal(vst.coords['dim_0'].values,raw.iloc[:,0].values)
    assert np.array_equal(vst.coords['dim_1'].values,expected)
    positive = (x>0).all(axis=1)
    log_x = np.log(x[positive])
    sf = np.exp(np.median(log_x-log_x.mean(axis=1,keepdims=True),axis=0))
    sferr = float(np.max(np.abs(sf-np.array(archive['size_factors']))))
    assert sferr < 1e-12
    pd.DataFrame({'Sample':expected,'recalculated':sf,
                  'archived':np.array(archive['size_factors'])}).to_csv(args.out/'size_factors.csv',index=False)
    gm = archive['gene_mapping'].dropna(subset=['SYMBOL'])
    gm = gm[gm.SYMBOL.ne('')]
    unique = gm[~gm.duplicated('ENSEMBL',keep=False) & ~gm.duplicated('SYMBOL',keep=False)]
    assert unique.shape == (18572,2)
    for target in (archive['unique_mapping'],ref['unique_mapping']):
        assert np.array_equal(unique.to_numpy(),target.to_numpy())
    unique.to_csv(args.out/'mapping_reconstructed.csv',index=False)
    v = pd.DataFrame(np.array(vst),index=vst.coords['dim_0'].values,columns=expected)
    mapped = v.loc[unique.ENSEMBL].copy(); mapped.index = unique.SYMBOL.to_numpy()
    sod2 = v.loc[['ENSG00000112096']].copy(); sod2.index=['SOD2']
    mapped = pd.concat([mapped,sod2])
    target = ref['gsva_matrix_sod2']
    assert np.array_equal(mapped.index.values,target.coords['dim_0'].values)
    maperr = float(np.max(np.abs(mapped.to_numpy()-np.array(target))))
    assert maperr == 0
    # Diagnostic only: infer a one-parameter transform from the archived zeros.
    # This checks algebraic consistency; it does NOT estimate dispersions from counts.
    y = np.array(vst); q = x/sf
    assert np.max(y[x==0])-np.min(y[x==0]) == 0
    k = float(4*2**y[x==0][0])
    fitted_to_archive = np.log2((np.sqrt(q)+np.sqrt(q+k))**2/4)
    diagnostic_err = float(np.max(np.abs(fitted_to_archive-y)))
    summary = dict(status='STRUCTURAL_CHECKS_PASS_NATIVE_VST_PENDING',genes=20345,samples=9,
        size_factor_max_absolute_error=sferr,mapping_rows=18572,
        mapping_rule='Drop missing/empty symbols; simultaneously remove all duplicated ENSEMBL or SYMBOL rows',
        corrected_genes=18573,SOD2_ENSEMBL='ENSG00000112096',
        frozen_vst_to_symbol_matrix_max_error=maperr,
        parametric_diagnostic_inferred_k=k,parametric_diagnostic_max_error=diagnostic_err,
        diagnostic_limit='k was inferred from the archived VST. Not independent dispersion fitting or VST reproduction.',
        native_R_status='Not executed in this preparation environment',
        annotation_scope='Rebuilt from full frozen annotation; not a new annotation-database query',
        sources=[dict(filename=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths])
    (args.out/'tissue_structure_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
