"""Frozen-rule local instrument rebuild; never reads outcome associations."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
import pandas as pd

ROUND = Path(__file__).resolve().parents[1]
ROOT = ROUND.parent
PREV = ROOT / 'round02_2026-09-18'
RUNTIME = Path('C:/Users/Lenovo/IBD_rebuild_runtime/round03_2026-09-27')
RELAXABLE = {'PAPER_COHORT_P_QC_RULE', 'COHORT_P_ORDER_ROUNDING_AMBIGUOUS'}


def sha(path):
    with Path(path).open('rb') as inp:
        return hashlib.file_digest(inp, 'sha256').hexdigest()


def eligible(reasons, relaxed=False):
    reasons = set(filter(None, str(reasons).split(';')))
    return not (reasons - (RELAXABLE if relaxed else set()))


def run_checked(cmd, cwd, output):
    result = subprocess.run([str(x) for x in cmd], cwd=cwd, capture_output=True,
                            encoding='utf-8', errors='replace', timeout=1800)
    output.write_text(result.stdout + '\nSTDERR:\n' + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'Command failed {result.returncode}: {output}')
    return result.stdout


def merge_windows(records, radius=500000):
    merged = []
    for chrom, pos, ident in sorted(records):
        start, end = max(1, pos-radius), pos+radius
        if merged and merged[-1]['chrom37'] == chrom and start <= merged[-1]['end37']:
            merged[-1]['end37'] = max(end, merged[-1]['end37'])
            merged[-1]['index_sources'].append(ident)
        else:
            merged.append(dict(chrom37=chrom, start37=start, end37=end, index_sources=[ident]))
    for i, item in enumerate(merged, 1):
        item['region_id'] = f'REGION_{i:03d}'
        item['index_sources'] = ';'.join(sorted(set(item['index_sources'])))
    return merged


def main():
    for folder in ['logs', 'tables', 'reports', 'metadata']:
        (ROUND / folder).mkdir(exist_ok=True)
    RUNTIME.mkdir(parents=True, exist_ok=True)
    lockpath = ROUND / 'config/analysis_lock.json'
    lock = json.loads(lockpath.read_text(encoding='utf-8'))
    existing = ROUND / 'metadata/execution_lock.json'
    if existing.exists():
        assert json.loads(existing.read_text())['analysis_lock_sha256'] == sha(lockpath), 'Lock changed; create new version'
    else:
        existing.write_text(json.dumps(dict(analysis_lock_sha256=sha(lockpath),
            frozen_before_execution_utc=datetime.now(timezone.utc).isoformat()), indent=2), encoding='utf-8')
    subprocess.run([sys.executable, '-X', 'utf8', str(PREV/'scripts/verify_round.py'), '--manifest-only'], check=True)
    ref = json.loads((PREV/'reports/public_ld_cache_audit.json').read_text())
    prefix = Path(ref['cache_prefix'])
    for ext in ['bed', 'bim', 'fam']:
        assert sha(str(prefix)+'.'+ext) == ref['local_sha256'][ext], f'Reference changed: {ext}'
    plink = prefix.parents[3] / 'tools/plink.exe'
    version = run_checked([plink, '--version'], prefix.parent, ROUND/'logs/plink_version.txt').strip()
    metadata = dict(plink_path=str(plink), plink_sha256=sha(plink), version=version,
                    public_reference_read_only=True, reference_sha256=ref['local_sha256'],
                    analysis_lock_sha256=sha(lockpath), runtime=str(RUNTIME), commands=[])
    all_summary, region_records = [], []
    for trait, ident in lock['exposures'].items():
        auditpath = PREV / 'tables' / f'{ident}_preclump_candidate_audit.csv.gz'
        data = pd.read_csv(auditpath, keep_default_na=False)
        assert not data.MarkerName.duplicated().any()
        assert (data['P.value'] < lock['primary_p_threshold']).all()
        data['mhc'] = (data.chrom37 == 6) & data.pos37.between(25000000, 34000000)
        for relaxed in [False, True]:
            for no_mhc in [False, True]:
                branch = ('relaxed_cohort_p' if relaxed else 'primary') + ('_nonMHC' if no_mhc else '_all')
                allowed = data.preflight_hold_reasons.map(lambda x: eligible(x, relaxed))
                if no_mhc:
                    allowed &= ~data.mhc
                selected = data.loc[allowed].copy().sort_values(['P.value', 'chrom37', 'pos37', 'reference_rsid'], kind='stable')
                assert not selected.reference_rsid.duplicated().any()
                stem = f'{ident}_{branch}'
                infile = RUNTIME / f'{stem}_clump_input.tsv'
                selected[['reference_rsid', 'P.value']].rename(columns={'reference_rsid':'SNP', 'P.value':'P'}).to_csv(infile, sep='\t', index=False)
                out = RUNTIME/stem
                cmd = [plink, '--bfile', prefix.name, '--clump', infile, '--clump-p1', str(lock['primary_p_threshold']),
                       '--clump-p2', str(lock['clump_p2']), '--clump-r2', str(lock['clump_r2']), '--clump-kb', str(lock['clump_kb']),
                       '--memory', str(lock['clump_memory_mb']), '--threads', str(lock['clump_threads']), '--out', out]
                run_checked(cmd, prefix.parent, ROUND/'logs'/f'{stem}_stdout.txt')
                clumped = pd.read_csv(str(out)+'.clumped', sep=r'\s+')
                assert not clumped.SNP.duplicated().any()
                assert set(clumped.SNP).issubset(set(selected.reference_rsid))
                instruments = selected.loc[selected.reference_rsid.isin(clumped.SNP)].copy()
                instruments['trait'] = trait
                instruments['gwas_id'] = ident
                instruments['branch'] = branch
                instruments['analysis_status'] = 'LD_CLUMPED_EXPOSURE_ONLY_NOT_HARMONISED_WITH_OUTCOME'
                instruments.to_csv(ROUND/'tables'/f'{stem}_instruments.csv', index=False)
                trace = data.copy()
                trace['eligible_for_branch'] = allowed
                trace['selected_index'] = trace.reference_rsid.isin(clumped.SNP) & allowed
                trace['branch'] = branch
                trace['branch_status'] = 'EXCLUDED_BY_QC'
                trace.loc[allowed, 'branch_status'] = 'CLUMPED_NON_INDEX'
                if no_mhc:
                    trace.loc[trace.mhc, 'branch_status'] = 'MHC_EXCLUDED_BEFORE_CLUMP'
                trace.loc[trace.selected_index, 'branch_status'] = 'SELECTED_INDEX'
                trace.to_csv(ROUND/'tables'/f'{stem}_candidate_trace.csv.gz', index=False, compression='gzip')
                for ext in ['clumped', 'log']:
                    (ROUND/'logs'/f'{stem}.{ext}').write_bytes(Path(str(out)+'.'+ext).read_bytes())
                (ROUND/'tables'/infile.name).write_bytes(infile.read_bytes())
                record = dict(trait=trait, gwas_id=ident, branch=branch, nominal_gws=len(data), eligible=len(selected),
                              clumped_index=len(instruments), mhc_index=int(instruments.mhc.sum()),
                              min_F=float(instruments.F_beta_over_SE_squared.min()), median_F=float(instruments.F_beta_over_SE_squared.median()))
                all_summary.append(record)
                metadata['commands'].append(dict(name=stem, command=[str(x) for x in cmd], cwd=str(prefix.parent),
                    input_sha256=sha(infile), source_candidate_sha256=sha(auditpath)))
                print(json.dumps(record), flush=True)
                if branch == 'primary_all':
                    region_records.extend((int(r.chrom37), int(r.pos37), f'{ident}:{r.reference_rsid}') for r in instruments.itertuples())
    pd.DataFrame(all_summary).to_csv(ROUND/'tables/instrument_rebuild_summary.csv', index=False)
    pd.DataFrame(merge_windows(region_records)).to_csv(ROUND/'tables/systematic_coloc_region_plan_GRCh37.csv', index=False)
    (ROUND/'metadata/instrument_execution.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print('PASS: exposure instruments built; no outcome estimates or MR conclusions generated.', flush=True)


if __name__ == '__main__':
    main()
