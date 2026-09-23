#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild CDD–HTN shared SMR screen and audit archived 412-probe coloc.

No GWAS, SMR, coloc or SuSiE is rerun. Python 3 standard library only.
Run inside the extracted package: python 06_code/reproduce_SMR_and_coloc.py
Outputs are written to 05_results_recalculated (can override with --out).
"""
import argparse, csv, hashlib, json, math, sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROBES = {"ENSG00000076685": "NT5C2", "ENSG00000141452": "RMC1", "ENSG00000213277": "MARCKSL1P1", "ENSG00000148842": "CNNM2"}
PRIOR_ORDER = [1e-6,1e-5,5e-5]
PAIRS = ('EQTL_HTN','EQTL_CDD','HTN_CDD')

def load_smr(path):
    with path.open('r', encoding='utf-8-sig') as f:
        rdr = csv.DictReader(f, delimiter=' ', skipinitialspace=True)
        # csv with delimiter=' ' does not handle arbitrary tabs / multiple spaces.
        # The original .smr is whitespace-separated; parse each line with split().
    with path.open('r', encoding='utf-8-sig') as f:
        header = f.readline().split()
        required = {'probeID','p_SMR','p_HEIDI','b_SMR','se_SMR','nsnp_HEIDI'}
        if not required.issubset(header):
            raise ValueError(f'Missing SMR fields in {path}: {required-set(header)}')
        records = {}
        for lineno, line in enumerate(f,2):
            if not line.strip(): continue
            cells = line.split()
            if len(cells) != len(header):
                raise ValueError(f'{path}:{lineno} unexpected field count {len(cells)} vs {len(header)}')
            row = dict(zip(header, cells))
            pid = row['probeID']
            if pid in records: raise ValueError(f'Duplicate probeID {pid} in {path}')
            try: p = float(row['p_SMR'])
            except Exception as ex: raise ValueError(f'{path}:{lineno} invalid p_SMR: {ex}')
            if not math.isfinite(p) or not 0 <= p <= 1:
                raise ValueError(f'{path}:{lineno} p_SMR outside [0,1]')
            records[pid] = row
        return records

def load_tsv(path):
    with path.open('r',encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader();w.writerows(rows)

def bh(pvals):
    n=len(pvals); order=sorted(range(n),key=lambda i:(pvals[i],i)); out=[None]*n; minimum=1.
    for k in range(n-1,-1,-1):
        i=order[k]; minimum=min(minimum, pvals[i]*n/(k+1));out[i]=minimum
    return out

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--out',type=Path,default=None)
    a=ap.parse_args()
    root=a.root.resolve(); out=(a.out or root/'05_results_recalculated').resolve();out.mkdir(parents=True,exist_ok=True)
    cdd_path=root/'01_primary_inputs/CDD_eQTLGen_all.smr'
    htn_path=root/'01_primary_inputs/HTN_eQTLGen_all.smr'
    batch_path=root/'03_historical_selection/batch_coloc_412_all_priors.tsv'
    frozen_path=root/'03_historical_selection/candidate_freeze_2026-09-13.tsv'
    cdd=load_smr(cdd_path);htn=load_smr(htn_path)
    common=sorted(set(cdd)&set(htn)); pvals=[max(float(cdd[p]['p_SMR']),float(htn[p]['p_SMR'])) for p in common]
    qvals=bh(pvals)
    allrows=[]
    for pid,p,q in zip(common,pvals,qvals):
        cr,hr=cdd[pid],htn[pid]
        r={'probeID':pid,'candidate_symbol_if_final':PROBES.get(pid,''),'ProbeChr_CDD':cr['ProbeChr'],'Probe_bp_CDD':cr['Probe_bp'],
           'topSNP_CDD':cr['topSNP'],'topSNP_HTN':hr['topSNP'],
           'b_SMR_CDD':cr['b_SMR'],'se_SMR_CDD':cr['se_SMR'],'p_SMR_CDD':cr['p_SMR'],
           'p_HEIDI_CDD':cr['p_HEIDI'],'nsnp_HEIDI_CDD':cr['nsnp_HEIDI'],
           'b_SMR_HTN':hr['b_SMR'],'se_SMR_HTN':hr['se_SMR'],'p_SMR_HTN':hr['p_SMR'],
           'p_HEIDI_HTN':hr['p_HEIDI'],'nsnp_HEIDI_HTN':hr['nsnp_HEIDI'],
           'Pshared':f'{p:.15g}','q_BH_shared':f'{q:.15g}',
           'nominal_both':str(p<.05),
           'HEIDI_both_gt_001':str(all(r['p_HEIDI']!='NA' and float(r['p_HEIDI'])>.01 for r in (cr,hr)))}
        allrows.append(r)
    allrows.sort(key=lambda r:(float(r['Pshared']),r['probeID']))
    for i,r in enumerate(allrows,1):r['p_shared_order']=i
    cols=list(allrows[0])
    nominal=[r for r in allrows if r['nominal_both']=='True']
    final=[r for r in allrows if r['probeID'] in PROBES]
    for name,rs in [('shared_SMR_all.tsv',allrows),('shared_SMR_nominal.tsv',nominal),('four_candidates_SMR.tsv',final)]:
        write_tsv(out/name,rs,cols)
    if len(final)!=4:raise AssertionError('A named final probe was not found in SMR source')
    batch=load_tsv(batch_path)
    batch_by={}
    for r in batch:
        key=(r['probeID'],r['pair'],float(r['p12']))
        if key in batch_by: raise ValueError(f'Duplicate coloc key {key}')
        batch_by[key]=r
    probe_set={r['probeID'] for r in batch}
    prior_rows=[]; selections={}
    for prior in PRIOR_ORDER:
        selected=[]
        for probe in sorted(probe_set):
            vals=[batch_by.get((probe,pair,prior)) for pair in PAIRS]
            if any(v is None for v in vals): raise AssertionError(f'Missing coloc triple: {probe} {prior}')
            if all(float(v['H4'])>=.5 for v in vals):selected.append(probe)
        selections[prior]=selected
        prior_rows.append({'p12':f'{prior:.0e}','n_probes_total':len(probe_set),
                           'number_meeting_all_three_H4_ge_0.50':len(selected),
                           'probeIDs':';'.join(selected)})
    write_tsv(out/'coloc_prior_sensitivity_counts.tsv',prior_rows,list(prior_rows[0]))
    main_set=set(selections[1e-5]); frozen=load_tsv(frozen_path);freeze_ids={r['probeID'] for r in frozen}
    if main_set!=freeze_ids: raise AssertionError(f'Primary-prior selection differs from freeze table: {main_set^freeze_ids}')
    details=[]
    for pid in sorted(freeze_ids):
        row={'probeID':pid,'candidate_symbol_if_final':PROBES.get(pid,'')}
        for pair in PAIRS:
            row[f'H4_{pair}_primary']=batch_by[(pid,pair,1e-5)]['H4']
        row['min_H4_primary']=f"{min(float(row[f'H4_{p}_primary']) for p in PAIRS):.15g}"
        details.append(row)
    write_tsv(out/'four_candidates_screen_coloc.tsv',details,list(details[0]))
    summary={
      'scope':'reconstructed_from_saved_input_files; no SMR/coloc/SuSiE rerun',
      'CDD_rows':len(cdd),'HTN_rows':len(htn),'joint_probe_count':len(common),
      'CDD_only_count':len(set(cdd)-set(htn)),'HTN_only_count':len(set(htn)-set(cdd)),
      'nominal_shared_Pshared_lt_0.05':len(nominal),
      'BH_FDR_lt_0.05':sum(float(r['q_BH_shared'])<.05 for r in allrows),
      'min_Pshared':min(pvals),'min_BH_q':min(qvals),
      'all_3_HEIDI_gt_0.01_not_a_selection_rule':None,
      'nominal_HEIDI_both_gt_0.01_descriptive':sum(r['HEIDI_both_gt_001']=='True' for r in nominal),
      'batch_rows':len(batch),'batch_probe_count':len(probe_set),
      'prior_pass_counts':{f'{p:.0e}':len(ids) for p,ids in selections.items()},
      'selected_primary_probeIDs':sorted(main_set),
      'freeze_probeIDs_match_selected_primary':main_set==freeze_ids,
      'source_SHA256':{'CDD_eQTLGen_all.smr':sha256(cdd_path),'HTN_eQTLGen_all.smr':sha256(htn_path),
                       'batch_coloc_412_all_priors.tsv':sha256(batch_path),
                       'candidate_freeze_2026-09-13.tsv':sha256(frozen_path)},
      'limitations':['Cannot prove that threshold was defined before seeing data','Cannot establish chronology of freeze, formal colocalization and MVP from these records alone','Three pairwise H4 values do not prove one common causal variant across all three traits','Formal coloc of selected candidates is not an independent replication']}
    (out/'reconstruction_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS source_rows CDD=%s HTN=%s joint=%s nominal=%s FDR_significant=%s' % (len(cdd),len(htn),len(common),len(nominal),summary['BH_FDR_lt_0.05']))
    print('PASS coloc batch_rows=%s unique_probes=%s prior_counts=%s' % (len(batch),len(probe_set),summary['prior_pass_counts']))
    print('PASS freeze_primary_match=%s final_probes=%s' % (summary['freeze_probeIDs_match_selected_primary'],','.join(sorted(main_set))))
    print('OUTPUT',out)

if __name__=='__main__':main()
