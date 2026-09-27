from pathlib import Path
import json
import math
import pandas as pd
from harmonise_core import ivw_crosscheck,bh_with_missing

ROUND=Path(__file__).resolve().parents[1]
results=pd.read_csv(ROUND/'tables/mr_all_locked_branches.csv')
ivw=results[results.method=='IVW multiplicative random effects']
assert len(ivw)==60
checks=[]
for row in ivw.itertuples():
    d=pd.read_csv(ROUND/'tables'/f'harmonised_{row.gwas_id}_{row.branch}_to_{row.outcome}.csv')
    x=ivw_crosscheck(d.beta_exposure,d.beta_outcome,d.se_outcome)
    for field,rfield in [('b','beta'),('se','se'),('pval','p'),('Q','Q'),('Q_df','Q_df'),('Q_pval','Q_p')]:
        assert math.isclose(x[field],getattr(row,rfield),rel_tol=2e-10,abs_tol=1e-12),(row.exposure,row.branch,row.outcome,field)
    checks.append(dict(exposure=row.exposure,branch=row.branch,outcome=row.outcome,status='PASS'))
primary=pd.read_csv(ROUND/'tables/mr_primary12.csv')
assert len(primary)==12 and len(primary[['exposure','outcome']].drop_duplicates())==12
expected=bh_with_missing(primary.p.to_numpy())
assert all(math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15) for a,b in zip(expected,primary.FDR_BH_primary12))
loo=pd.read_csv(ROUND/'tables/mr_primary12_leave_one_out.csv')
assert len(loo)==sum(primary.nsnp)
assert loo[['gwas_id','outcome','omitted']].duplicated().sum()==0
report=dict(status='PASS',ivw_pairs_checked=len(checks),primary_family=12,leave_one_out_rows=len(loo),
    statement='Independent Python formula matches installed TwoSampleMR IVW beta, SE, p, Q and Q p for all 60 branches. This validates arithmetic, not MR assumptions.')
(ROUND/'reports/actual_MR_crosslanguage_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report),flush=True)
