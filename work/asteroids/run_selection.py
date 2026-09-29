"""Run audited equal-weight selection from completed, independent candidate pools."""
import argparse,csv,hashlib,json,time
from pathlib import Path
import numpy as np
from selection import Objective,select
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'outputs/Tau5_joint_refinement'
OUT=ROOT/'outputs/Tau5_ASTEROIDS'
SCALES={'CA':.5,'C':.5,'CB':.5,'N':2.45,'NH':.49,'HA':.25}
AA3=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
def save(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def table(path,header,rows):
    with path.open('w') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)
def load(sample,rep):
    d=BASE/sample;names=[f'pool_{rep}',f'pool_{rep+2}',f'helix_lengths_{rep}'];counts=[250,1750,110 if sample=='WT' else 90]
    paths=[]
    for name,count in zip(names,counts):
        source=sorted((d/name).glob('conformer_*.pdb'))
        found=sorted((d/'backcalc'/name).glob('conformer_*/result.json'))
        if len(source)!=count or len(found)!=count:raise ValueError(f'{sample}/{name}: requires {count} complete candidates; found {len(source)} structures, {len(found)} predictions')
        if {p.stem for p in source}!={p.parent.name for p in found}:raise ValueError('Candidate ID mismatch')
        paths.extend(found)
    records=[json.loads(p.read_text()) for p in paths]
    sequence=''.join(l.strip() for l in (ROOT/f'outputs/ASTEROIDS_setup/{sample}/{sample}_120aa.fasta').read_text().splitlines() if not l.startswith('>'))
    for r in records:
        pdb=ROOT/r['pdb'];seq=''.join(AA3[l[17:20].strip()] for l in pdb.read_text().splitlines() if l.startswith('ATOM') and l[12:16].strip()=='CA')
        if seq!=sequence:raise ValueError(f'Sequence mismatch: {pdb}')
        if len(r['dssp'])!=len(sequence):raise ValueError('DSSP length mismatch')
    restraints=[r for r in json.loads((ROOT/f'outputs/ASTEROIDS_setup/{sample}/chemical_shift_restraints.json').read_text()) if 3<=r['model_residue']<=119]
    keys=[f'{r["model_residue"]}:{r["atom"]}' for r in restraints]
    # Missing predictions fail explicitly; never silently discard a measured nucleus.
    cs=np.array([[r['shifts'][k] for k in keys] for r in records]);rc=np.array([[r['random_coil_shifts'][k] for k in keys] for r in records])
    if np.max(np.ptp(rc,axis=0))>1e-6:raise ValueError('Inconsistent sequence-specific random-coil reference')
    observed=np.array([r['shift_ppm'] for r in restraints]);atoms=np.array([r['atom'] for r in restraints]);scales=np.array([SCALES[a] for a in atoms])
    q,y,err=np.loadtxt(ROOT/f'outputs/ASTEROIDS_setup/{sample}/{sample}_EOM_fit_range.dat',unpack=True)
    if not np.all(err>0):raise ValueError('Nonpositive SAXS errors')
    for r in records:
        if min(q)<min(r['q'])-1e-10 or max(q)>max(r['q'])+1e-10:raise ValueError('SAXS extrapolation prohibited')
    curves=np.array([np.interp(q,r['q'],r['intensity']) for r in records])
    fingerprint=hashlib.sha256(b''.join(p.read_bytes() for p in paths)+observed.tobytes()+y.tobytes()+err.tobytes()).hexdigest()
    return records,restraints,cs,rc[0],observed,atoms,scales,q,y,err,curves,fingerprint

def run(sample,rep,generations=1000,population=100,seeds=(20260929,20261029,20261129)):
    records,restraints,cs,rc,observed,atoms,scales,q,y,err,curves,fingerprint=load(sample,rep)
    out=OUT/sample/f'replicate_{rep}';out.mkdir(parents=True,exist_ok=True)
    obj=Objective(cs,observed,scales,curves,y,err)
    manifest=dict(sample=sample,replicate=rep,n_candidates=len(records),n_selected=500,input_sha256=fingerprint,
        shift_scales_ppm=SCALES,objective='sum squared scaled shift residuals + sum squared SAXS residuals / experimental variance',
        note='Carbon 0.5 ppm is an optimization scale, not an acceptance tolerance or a calibrated model error. No shift offsets, no additive SAXS background, no continuous conformer weights.',
        population=population,generations=generations,seeds=list(seeds),independent_validation=False)
    save(out/'run_manifest.json',manifest)
    table(out/'candidates.csv',['candidate_index','source_pdb','protonated_pdb'],[[i,r['pdb'],r['protonated_pdb']] for i,r in enumerate(records)])
    dssp=np.array([r['dssp'] for r in records]);profiles={'DSSP_H':dssp=='H','DSSP_HGI':np.isin(dssp,['H','G','I'])}
    def metrics(ids):
        pred=cs[ids].mean(axis=0);curve=curves[ids].mean(axis=0);scale=max(0,float(np.dot(curve/err,y/err)/np.dot(curve/err,curve/err)))
        residual=pred-observed
        return dict(RMSD_ppm={a:float(np.sqrt(np.mean(residual[atoms==a]**2))) for a in sorted(set(atoms))},SAXS_chi2_mean=float(np.mean(((scale*curve-y)/err)**2)),SAXS_scale=scale)
    raw=metrics(np.arange(len(records)))
    for seed in seeds:
        dest=out/f'seed_{seed}';dest.mkdir(exist_ok=True)
        done=dest/'summary.json'
        if done.exists():
            old=json.loads(done.read_text())
            if old['input_sha256']==fingerprint and old['generations']==generations and old['population']==population:continue
            raise ValueError('Existing run has different inputs/settings; use a new output location')
        start=time.time()
        def checkpoint(gen,ids,score,trace):
            save(dest/'checkpoint.json',dict(generation=gen,selected_indices=ids.tolist(),objective=float(score),trace=trace,input_sha256=fingerprint,seed=seed,updated_unix=time.time(),restart_note='Best-so-far checkpoint; restarting reruns this deterministic seed from generation zero'))
            print(sample,rep,seed,gen,float(score),flush=True)
        ids,trace=select(obj,500,seed=seed,population=population,generations=generations,callback=checkpoint)
        fitted=metrics(ids)
        table(dest/'selected.csv',['candidate_index','source_pdb','protonated_pdb','weight'],[[int(i),records[i]['pdb'],records[i]['protonated_pdb'],1/500] for i in ids])
        rawcs=cs.mean(axis=0);fitcs=cs[ids].mean(axis=0)
        table(dest/'shifts.csv',['AR_residue','atom','experimental_secondary_ppm','pool_secondary_ppm','selected_secondary_ppm','Exp_minus_pool_ppm','Exp_minus_selected_ppm'],[[r['author_residue'],r['atom'],observed[j]-rc[j],rawcs[j]-rc[j],fitcs[j]-rc[j],observed[j]-rawcs[j],observed[j]-fitcs[j]] for j,r in enumerate(restraints)])
        table(dest/'saxs.csv',['q_Ainv','I_exp','sigma_exp','I_pool','I_selected'],zip(q,y,err,raw['SAXS_scale']*curves.mean(axis=0),fitted['SAXS_scale']*curves[ids].mean(axis=0)))
        table(dest/'helicity.csv',['AR_residue','definition','pool','selected'],[[j+328,name,float(matrix[:,j].mean()),float(matrix[ids,j].mean())] for name,matrix in profiles.items() for j in range(2,119)])
        save(done,dict(**manifest,seed=seed,method='Independent ASTEROIDS-style equal-weight genetic subset selection',equal_weight=.002,unique_selected=len(set(ids)),pool=raw,selected=fitted,final_objective=float(trace[-1]),elapsed_seconds=time.time()-start,complete=True,convergence_established=False))
    save(out/'status.json',dict(state='selection_runs_complete',updated_unix=time.time(),convergence_assessment_pending=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('sample',choices=['WT','AA']);ap.add_argument('replicate',type=int,choices=[1,2]);ap.add_argument('--generations',type=int,default=1000);ap.add_argument('--population',type=int,default=100);a=ap.parse_args();run(a.sample,a.replicate,a.generations,a.population)
