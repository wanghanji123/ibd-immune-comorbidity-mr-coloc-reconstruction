"""Fail closed on local inconsistencies; never equate test PASS with submission readiness."""
from pathlib import Path
import json
import os
import platform
import subprocess
import sys
import pandas as pd
from build_instruments import ROUND, PREV, sha


def call(cmd, filename, env=None):
    result = subprocess.run([str(x) for x in cmd], cwd=ROUND, env=env,
                            capture_output=True, encoding='utf-8', errors='replace', timeout=600)
    (ROUND/'logs'/filename).write_text(result.stdout+'\nSTDERR:\n'+result.stderr,encoding='utf-8')
    if result.returncode:
        raise RuntimeError(filename+': '+result.stderr[-1000:])
    print(result.stdout.strip() or result.stderr.strip(),flush=True)


def main():
    call([sys.executable,'-X','utf8','-m','unittest','discover','-s','scripts','-p','test_*.py','-v'],'unit_tests.txt')
    env = dict(os.environ,LC_ALL='C',LANG='C')
    call(['C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe','--vanilla','scripts/check_MR_runtime.R'], 'MR_runtime_checks.txt',env)
    call([sys.executable,'-X','utf8','scripts/check_MR_crosslanguage.py'],'MR_crosslanguage.txt')
    lock = json.loads((ROUND/'metadata/execution_lock.json').read_text())
    assert sha(ROUND/'config/analysis_lock.json') == lock['analysis_lock_sha256']
    source_checks=[]
    for ident in ['GCST004131','GCST004132','GCST004133']:
        qc=json.loads((PREV/'reports'/f'{ident}_full_qc.json').read_text())
        actual=sha(qc['path'])
        assert actual == qc['sha256'], 'Original exposure changed: '+ident
        source_checks.append(dict(gwas_id=ident,sha256=actual,match_previous_full_source_audit=True))
    report = json.loads((ROUND/'reports/clumping_validation.json').read_text())
    assert report['status'] == 'PASS' and len(report['checks']) == 12
    coord = pd.read_csv(ROUND/'tables/instrument_coordinate_audit.csv')
    assert len(coord) == report['n_unique_instrument_variants']
    assert not coord.MarkerName.duplicated().any()
    mapped=coord.loc[coord.reciprocal_exact]
    keys=[(r.chrom38,int(r.pos38),tuple(sorted([r.effect_allele38,r.other_allele38]))) for r in mapped.itertuples()]
    assert len(keys)==len(set(keys)), 'Different source markers collide at same target variant'
    # A later mapping failure must remain visible; it must not disappear from the report.
    failed = coord.loc[~coord.reciprocal_exact]
    failed.to_csv(ROUND/'tables/coordinate_mapping_holds.csv',index=False)
    summary = pd.read_csv(ROUND/'tables/instrument_rebuild_summary.csv')
    assert len(summary) == 12
    comparisons = []
    for ident in summary.gwas_id.unique():
        def ids(branch):
            return set(pd.read_csv(ROUND/'tables'/f'{ident}_{branch}_instruments.csv').reference_rsid)
        primary, nomhc, relaxed = ids('primary_all'),ids('primary_nonMHC'),ids('relaxed_cohort_p_all')
        for mode,other in [('nonMHC',nomhc),('relaxed_cohort_p',relaxed)]:
            for snp in sorted(primary|other):
                comparisons.append(dict(gwas_id=ident,comparison=mode,SNP=snp,in_primary=snp in primary,
                                        in_sensitivity=snp in other))
    pd.DataFrame(comparisons).to_csv(ROUND/'tables/instrument_branch_comparison.csv',index=False)
    regions = pd.read_csv(ROUND/'tables/systematic_coloc_region_plan_GRCh37.csv')
    assert not regions.region_id.duplicated().any()
    validation = dict(status='LOCAL_INPUT_PREPARATION_PASS_NOT_SUBMISSION_READY',
        python=platform.python_version(),unit_tests=13,synthetic_R_Python_IVW_fixtures=4,
        instrument_branches=12,unique_variants=len(coord),coordinate_holds=len(failed),
        planned_coloc_regions=len(regions),actual_MR_run=False,actual_coloc_run=False,
        outcome_alleles_verified=False,finngen_registration='NOT_COMPLETED_USER_CONFIRMED',
        original_round02_manifest_unchanged=True,original_source_checks=source_checks)
    call([sys.executable,'-X','utf8',PREV/'scripts/verify_round.py','--manifest-only'],'round02_readonly_verification.txt')
    (ROUND/'reports/round03_validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    print(json.dumps(validation),flush=True)


if __name__ == '__main__':
    main()
