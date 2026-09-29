"""Wait for full forward-model caches, then run equal-weight selection once per pool."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'outputs/Tau5_ASTEROIDS';OUT.mkdir(exist_ok=True)
PENDING=[(s,r) for s in ['WT','AA'] for r in [1,2]]
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENMM_CPU_THREADS='1',MPLCONFIGDIR=str(ROOT/'work/mpl_cache'))
def write(state,**kwargs):
    dest=OUT/'selection_pipeline_status.json';tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(dict(state=state,updated_unix=time.time(),pending=PENDING,**kwargs),indent=2));tmp.replace(dest)
try:
    while PENDING:
        counts={};ready=None
        for sample,rep in PENDING:
            base=ROOT/'outputs/Tau5_joint_refinement'/sample/'backcalc'
            specs=[(f'pool_{rep}',250),(f'pool_{rep+2}',1750),(f'helix_lengths_{rep}',110 if sample=='WT' else 90)]
            ns=[len(list((base/name).glob('conformer_*/result.json'))) for name,n in specs]
            counts[f'{sample}_{rep}']=dict(zip([n for n,k in specs],ns))
            if ns==[n for name,n in specs] and ready is None:ready=(sample,rep)
        write('waiting_for_complete_predictions',counts=counts)
        if ready:
            sample,rep=ready;write('selecting',sample=sample,replicate=rep,counts=counts)
            with (OUT/f'{sample}_{rep}_selection.log').open('a') as f:
                subprocess.run([sys.executable,ROOT/'work/asteroids/run_selection.py',sample,str(rep)],cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT,check=True)
            PENDING.remove(ready)
        else:time.sleep(30)
    write('selection_runs_complete',convergence_assessment_pending=True)
except BaseException as e:
    write('failed',error=repr(e));raise
