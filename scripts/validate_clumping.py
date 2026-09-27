"""Reconcile every eligible SNP and independently check residual lead-pair LD."""
from collections import Counter
import json
import re
import pandas as pd
from build_instruments import ROUND, PREV, RUNTIME, run_checked, sha
from pathlib import Path


def membership(clumped):
    records = []
    for row in clumped.itertuples():
        records.append(dict(SNP=row.SNP, index_SNP=row.SNP, role='INDEX'))
        secondary = [] if row.SP2 == 'NONE' else row.SP2.split(',')
        assert len(secondary) == row.TOTAL, (row.SNP, row.TOTAL, len(secondary))
        for cell in secondary:
            match = re.fullmatch(r'(.+)\(1\)', cell)
            assert match, cell
            records.append(dict(SNP=match.group(1), index_SNP=row.SNP, role='SECONDARY'))
    assert max(Counter(r['SNP'] for r in records).values()) == 1
    return pd.DataFrame(records)


def main():
    execution = json.loads((ROUND/'metadata/instrument_execution.json').read_text())
    reference = json.loads((PREV/'reports/public_ld_cache_audit.json').read_text())
    prefix = Path(reference['cache_prefix'])
    plink = Path(execution['plink_path'])
    assert sha(plink) == execution['plink_sha256']
    checks, union = [], []
    for row in pd.read_csv(ROUND/'tables/instrument_rebuild_summary.csv').itertuples():
        stem = f'{row.gwas_id}_{row.branch}'
        input_data = pd.read_csv(ROUND/'tables'/f'{stem}_clump_input.tsv', sep='\t')
        clumped = pd.read_csv(ROUND/'logs'/f'{stem}.clumped', sep=r'\s+', keep_default_na=False)
        member = membership(clumped)
        assert set(input_data.SNP) == set(member.SNP), f'Unaccounted candidates: {stem}'
        member.to_csv(ROUND/'tables'/f'{stem}_membership.csv.gz', index=False, compression='gzip')
        leads = pd.read_csv(ROUND/'tables'/f'{stem}_instruments.csv', keep_default_na=False)
        union.append(leads[['MarkerName','chrom37','pos37','effect_allele','other_allele','reference_rsid']])
        ids = RUNTIME/f'{stem}_lead_ids.txt'
        ids.write_text('\n'.join(leads.reference_rsid)+'\n', encoding='ascii')
        out = RUNTIME/f'{stem}_lead_ld'
        cmd = [plink, '--bfile', prefix.name, '--extract', ids, '--r2', 'dprime',
               '--ld-window', '999999', '--ld-window-kb', '10000', '--ld-window-r2', '0',
               '--memory', '2048', '--threads', '1', '--out', out]
        run_checked(cmd, prefix.parent, ROUND/'logs'/f'{stem}_lead_ld_stdout.txt')
        ld = pd.read_csv(str(out)+'.ld', sep=r'\s+')
        assert not ld.R2.isna().any(), f'Unestimable pair LD: {stem}'
        # PLINK table output is rounded. Tolerance is less than its last reported digit near the cutoff.
        violations = ld.loc[ld.R2 > .00100001]
        assert violations.empty, f'Residual LD violates locked clump rule: {stem}'
        positions = list(leads[['chrom37','pos37']].itertuples(index=False, name=None))
        expected_pairs = sum(c1 == c2 and abs(p1-p2) <= 10000000
                             for i, (c1,p1) in enumerate(positions) for c2,p2 in positions[i+1:])
        assert len(ld) == expected_pairs, f'LD pair coverage mismatch: {stem}'
        ld.to_csv(ROUND/'tables'/f'{stem}_lead_pair_ld.csv', index=False)
        result = dict(gwas_id=row.gwas_id, branch=row.branch, eligible_reconciled=len(member),
                      leads=len(leads), within_10Mb_pairs_checked=len(ld), residual_LD_violations=len(violations),
                      max_R2=float(ld.R2.max()) if len(ld) else None)
        checks.append(result)
        print(json.dumps(result), flush=True)
    variants = pd.concat(union).drop_duplicates()
    assert not variants.MarkerName.duplicated().any(), 'Conflicting mapping of same source marker'
    variants.sort_values(['chrom37','pos37']).to_csv(ROUND/'tables/unique_instrument_variants_GRCh37.csv', index=False)
    (ROUND/'reports/clumping_validation.json').write_text(json.dumps(dict(checks=checks,
        n_unique_instrument_variants=len(variants), status='PASS',
        note='Pairwise LD verified within the prespecified 10Mb window; this is not proof of zero long-range LD.'), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
