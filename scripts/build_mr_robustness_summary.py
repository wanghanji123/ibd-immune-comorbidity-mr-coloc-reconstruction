from pathlib import Path
import json
import math
import pandas as pd

ROUND=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROUND/'tables/mr_all_locked_branches.csv')
loo=pd.read_csv(ROUND/'tables/mr_primary12_leave_one_out.csv')
outcomes=['M13_ANKYLOSPON','H7_IRIDOCYCLITIS','L12_PSORIASIS','M13_PSORIARTH']


def one(exposure,outcome,branch,method):
    x=d[(d.exposure==exposure)&(d.outcome==outcome)&(d.branch==branch)&(d.method==method)]
    assert len(x)==1
    return x.iloc[0]


rows=[]
for exposure in ['IBD','CD','UC']:
  for outcome in outcomes:
    ivw=one(exposure,outcome,'primary_all','IVW multiplicative random effects')
    med=one(exposure,outcome,'primary_all','Weighted median')
    egger=one(exposure,outcome,'primary_all','MR-Egger')
    nomhc=one(exposure,outcome,'primary_nonMHC','IVW multiplicative random effects')
    relaxed=one(exposure,outcome,'relaxed_cohort_p_all','IVW multiplicative random effects')
    ld=loo[(loo.exposure==exposure)&(loo.outcome==outcome)]
    assert len(ld)==ivw.nsnp
    strict=one(exposure,'M13_ANKYLOSPON_STRICT','primary_all','IVW multiplicative random effects') if outcome=='M13_ANKYLOSPON' else None
    most=ld.iloc[(ld.beta-ivw.beta).abs().argmax()]
    sign=lambda x: 0 if x==0 else (1 if x>0 else -1)
    rows.append(dict(exposure=exposure,outcome=outcome,nsnp=int(ivw.nsnp),primary_OR=ivw.OR,primary_CI_low=ivw.CI_low,
      primary_CI_high=ivw.CI_high,primary_p=ivw.p,primary_FDR=ivw.FDR_BH_primary12,heterogeneity_Q=ivw.Q,
      heterogeneity_Q_df=int(ivw.Q_df),heterogeneity_Q_p=ivw.Q_p,
      weighted_median_OR=med.OR,weighted_median_p=med.p,weighted_median_same_direction=sign(med.beta)==sign(ivw.beta),
      egger_OR=egger.OR,egger_p=egger.p,egger_same_direction=sign(egger.beta)==sign(ivw.beta),
      egger_intercept=egger.egger_intercept,egger_intercept_p=egger.egger_intercept_p,
      nonMHC_nsnp=int(nomhc.nsnp),nonMHC_OR=nomhc.OR,nonMHC_p=nomhc.p,nonMHC_same_direction=sign(nomhc.beta)==sign(ivw.beta),
      relaxed_cohort_p_OR=relaxed.OR,relaxed_cohort_p_p=relaxed.p,relaxed_cohort_same_direction=sign(relaxed.beta)==sign(ivw.beta),
      loo_beta_min=ld.beta.min(),loo_beta_max=ld.beta.max(),loo_all_same_direction=bool((ld.beta.map(sign)==sign(ivw.beta)).all()),
      loo_max_deviation_marker=most.omitted,loo_max_abs_beta_deviation=abs(most.beta-ivw.beta),
      strict_AS_OR=(strict.OR if strict is not None else math.nan),strict_AS_p=(strict.p if strict is not None else math.nan),
      strict_AS_same_direction=(sign(strict.beta)==sign(ivw.beta) if strict is not None else None)))
frame=pd.DataFrame(rows)
frame['primary_FDR_positive']=frame.primary_FDR<.05
frame['heterogeneity_present']=frame.heterogeneity_Q_p<.05
frame['weighted_median_supports']=frame.weighted_median_same_direction&(frame.weighted_median_p<.05)
frame['nonMHC_supports']=frame.nonMHC_same_direction&(frame.nonMHC_p<.05)
frame['directional_pleiotropy_signal']=frame.egger_intercept_p<.05
frame['robustness_note']='PRIMARY_NULL_OR_INCONCLUSIVE'
mask=frame.primary_FDR_positive&frame.weighted_median_supports&frame.nonMHC_supports
frame.loc[mask,'robustness_note']='POSITIVE_ACROSS_IVW_WEIGHTED_MEDIAN_NONMHC_WITH_HETEROGENEITY'
frame.loc[mask&frame.directional_pleiotropy_signal,'robustness_note']='POSITIVE_BUT_DIRECTIONAL_PLEIOTROPY_SIGNAL'
frame.loc[frame.primary_FDR_positive&~frame.weighted_median_supports,'robustness_note']='PRIMARY_POSITIVE_METHOD_SENSITIVE'
frame.loc[~frame.primary_FDR_positive&frame.nonMHC_supports,'robustness_note']='PRIMARY_NULL_BUT_NONMHC_POSITIVE_MHC_SENSITIVE'
frame.to_csv(ROUND/'tables/mr_robustness_summary_prePRESSO.csv',index=False)
report=dict(primary_FDR_positive=int(frame.primary_FDR_positive.sum()),heterogeneity_present=int(frame.heterogeneity_present.sum()),
            weighted_median_supports=int(frame.weighted_median_supports.sum()),nonMHC_supports=int(frame.nonMHC_supports.sum()),
            directional_pleiotropy_signals=int(frame.directional_pleiotropy_signal.sum()),
            note_counts=frame.robustness_note.value_counts().to_dict(),
            caution='These are genetic-liability association diagnostics. Labels do not establish causality, clinical comorbidity risk, mechanism, or absence of pleiotropy.')
(ROUND/'reports/mr_robustness_prePRESSO.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report),flush=True)
