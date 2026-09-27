"""Increase MR-PRESSO precision only for seed-unstable primary pairs."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json
import os
import subprocess
import time

import pandas as pd


ROUND = Path(__file__).resolve().parents[1]
RSCRIPT = Path("C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe")
IDS = {"IBD": "GCST004131", "CD": "GCST004132", "UC": "GCST004133"}
SEEDS = [271828, 314159]
N_SIM = 10000


def run(job):
    exposure, outcome, seed = job
    ident = IDS[exposure]
    stem = f"{ident}__{outcome}__seed{seed}__n{N_SIM}"
    log = ROUND / "logs" / f"mrpresso_refinement__{stem}.log"
    cmd = [str(RSCRIPT), "--vanilla", "scripts/run_one_mrpresso.R", ident,
           exposure, outcome, str(seed), str(N_SIM)]
    env = dict(os.environ, LC_ALL="C", LANG="C")
    start = time.time()
    proc = subprocess.run(cmd, cwd=ROUND, env=env, capture_output=True,
                          encoding="utf-8", errors="replace", timeout=7200)
    log.write_text(proc.stdout + "\nSTDERR:\n" + proc.stderr, encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(stem + " failed: " + proc.stderr[-1000:])
    return dict(gwas_id=ident, exposure=exposure, outcome=outcome, seed=seed,
                n_sim=N_SIM, elapsed_seconds=time.time()-start, log=str(log), status="PASS")


def main():
    stability = pd.read_csv(ROUND / "tables/mrpresso_seed_stability.csv")
    selected = stability.loc[stability.needs_higher_precision.astype(bool), ["exposure", "outcome"]]
    jobs = [(row.exposure, row.outcome, seed) for row in selected.itertuples(index=False)
            for seed in SEEDS]
    records = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(run, job): job for job in jobs}
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(json.dumps(record), flush=True)
    records.sort(key=lambda x: (x["gwas_id"], x["outcome"], x["seed"]))
    report = dict(selection_rule="Only 3000-simulation pairs with nonidentical outlier sets across the two locked seeds",
                  n_selected_pairs=len(selected), n_jobs=len(jobs), n_sim_each=N_SIM,
                  records=records)
    (ROUND / "reports/mrpresso_refinement_execution_status.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    assert len(records) == len(jobs)
    print(json.dumps({"status": "PASS", "pairs": len(selected), "jobs": len(jobs)}), flush=True)


if __name__ == "__main__":
    main()
