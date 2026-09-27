"""Run frozen primary 12 MR-PRESSO jobs with two seeds and auditable precision."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import json
import os
import subprocess
import time

ROUND=Path(__file__).resolve().parents[1]
RSCRIPT=Path('C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe')
EXPOSURES={'GCST004131':'IBD','GCST004132':'CD','GCST004133':'UC'}
OUTCOMES=['M13_ANKYLOSPON','H7_IRIDOCYCLITIS','L12_PSORIASIS','M13_PSORIARTH']
SEEDS=[271828,314159]
N_SIM=3000


def run(job):
    ident,exposure,outcome,seed=job
    stem=f'{ident}__{outcome}__seed{seed}__n{N_SIM}'
    log=ROUND/'logs'/f'mrpresso__{stem}.log'
    cmd=[str(RSCRIPT),'--vanilla','scripts/run_one_mrpresso.R',ident,exposure,outcome,str(seed),str(N_SIM)]
    env=dict(os.environ,LC_ALL='C',LANG='C')
    start=time.time()
    proc=subprocess.run(cmd,cwd=ROUND,env=env,capture_output=True,encoding='utf-8',errors='replace',timeout=3600)
    log.write_text(proc.stdout+'\nSTDERR:\n'+proc.stderr,encoding='utf-8')
    if proc.returncode: raise RuntimeError(stem+' failed: '+proc.stderr[-1000:])
    return dict(gwas_id=ident,exposure=exposure,outcome=outcome,seed=seed,n_sim=N_SIM,
                elapsed_seconds=time.time()-start,log=str(log),status='PASS')


if __name__=='__main__':
    jobs=[(i,e,o,s) for i,e in EXPOSURES.items() for o in OUTCOMES for s in SEEDS]
    records=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(run,j):j for j in jobs}
        for future in as_completed(futures):
            record=future.result();records.append(record);print(json.dumps(record),flush=True)
    records.sort(key=lambda x:(x['gwas_id'],x['outcome'],x['seed']))
    (ROUND/'reports/mrpresso_execution_status.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    assert len(records)==24
    print('PASS: primary 12 pairs, two MR-PRESSO seeds each',flush=True)
