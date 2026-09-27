"""Build all frozen analysis branches from exact point queries with full trace."""
from collections import Counter
from pathlib import Path
import json
import math
import pandas as pd
from harmonise_core import align

ROUND=Path(__file__).resolve().parents[1]
PRIMARY=['M13_ANKYLOSPON','H7_IRIDOCYCLITIS','L12_PSORIASIS','M13_PSORIARTH']
ALL=PRIMARY+['M13_ANKYLOSPON_STRICT']
TRAITS={'GCST004131':'IBD','GCST004132':'CD','GCST004133':'UC'}
SAMPLES={
 'M13_ANKYLOSPON':(3838,353224),'H7_IRIDOCYCLITIS':(9958,473095),
 'L12_PSORIASIS':(12760,482181),'M13_PSORIARTH':(4290,315115),
 'M13_ANKYLOSPON_STRICT':(1705,496755)}


def number(row,name):
    try: value=float(row[name])
    except Exception: return math.nan
    return value


def main():
    raw=pd.read_csv(ROUND/'tables/outcome_instrument_point_query_raw.csv',dtype={'chrom':str},keep_default_na=False)
    coords=pd.read_csv(ROUND/'tables/instrument_coordinate_audit.csv',keep_default_na=False)
    coord={r.MarkerName:r for r in coords.itertuples()}
    assert set(raw.endpoint)==set(ALL)
    traces=[]; summary=[]
    for ident,trait in TRAITS.items():
      for branch in ['primary_all','primary_nonMHC','relaxed_cohort_p_all','relaxed_cohort_p_nonMHC']:
        instruments=pd.read_csv(ROUND/'tables'/f'{ident}_{branch}_instruments.csv',keep_default_na=False)
        assert not instruments.MarkerName.duplicated().any()
        for endpoint in ALL:
          counts=Counter(); aligned=[]
          for inst in instruments.itertuples():
            c=coord[inst.MarkerName]
            candidates=raw[(raw.endpoint==endpoint)&(raw.query_marker==inst.MarkerName)]
            accepted=[]; rejection=[]
            for out in candidates.to_dict('records'):
              result=align(c.effect_allele38,c.other_allele38,out['ref'],out['alt'],number(out,'beta'),number(out,'sebeta'),number(out,'af_alt'))
              result.update(out)
              if result['status']=='ALIGNED': accepted.append(result)
              else: rejection.append(result['status'])
            if len(accepted)==1:
              out=accepted[0]
              if not (0 < number(out,'pval') <= 1): status='EXCLUDE_INVALID_OUTCOME_P'
              elif str(out['chrom']) != str(c.chrom38).removeprefix('chr') or int(out['pos']) != int(c.pos38): status='EXCLUDE_QUERY_COORDINATE_MISMATCH'
              else: status='ALIGNED'
            elif len(accepted)>1: status='EXCLUDE_MULTIPLE_EXACT_ALLELE_ROWS'
            elif candidates.empty: status='EXCLUDE_OUTCOME_POSITION_ABSENT'
            elif rejection: status=';'.join(sorted(set(rejection)))
            else: status='EXCLUDE_UNCLASSIFIED'
            row=dict(exposure=trait,gwas_id=ident,branch=branch,outcome=endpoint,MarkerName=inst.MarkerName,
              reference_rsid=inst.reference_rsid,chrom37=inst.chrom37,pos37=inst.pos37,chrom38=c.chrom38,pos38=c.pos38,
              effect_allele_exposure=c.effect_allele38,other_allele_exposure=c.other_allele38,
              beta_exposure=inst.Effect_beta_unchanged,se_exposure=inst.SE,F_exposure=inst.F_beta_over_SE_squared,
              outcome_rows_at_position=len(candidates),matching_outcome_rows=len(accepted),harmonisation_status=status)
            if status=='ALIGNED':
              ncase,ncontrol=SAMPLES[endpoint]
              row.update(effect_allele_outcome=c.effect_allele38,other_allele_outcome=c.other_allele38,
                beta_outcome=out['beta_outcome'],se_outcome=out['se_outcome'],eaf_outcome=out['eaf_outcome'],
                p_outcome=number(out,'pval'),outcome_ref=out['ref'],outcome_alt=out['alt'],outcome_rsids=out['rsids'],
                ncase_outcome=ncase,ncontrol_outcome=ncontrol,n_outcome=ncase+ncontrol)
              aligned.append(row)
            traces.append(row); counts[status]+=1
          pair=pd.DataFrame(aligned)
          pair.to_csv(ROUND/'tables'/f'harmonised_{ident}_{branch}_to_{endpoint}.csv',index=False)
          summary.append(dict(exposure=trait,gwas_id=ident,branch=branch,outcome=endpoint,n_instruments=len(instruments),
             n_aligned=len(pair),n_excluded=len(instruments)-len(pair),status_counts=json.dumps(counts,sort_keys=True)))
          print(json.dumps(summary[-1]),flush=True)
    trace=pd.DataFrame(traces)
    trace.to_csv(ROUND/'tables/harmonisation_trace_all_branches.csv.gz',index=False,compression='gzip')
    pd.DataFrame(summary).to_csv(ROUND/'tables/harmonisation_summary.csv',index=False)
    assert len(summary)==3*4*5
    problematic=trace[trace.harmonisation_status=='EXCLUDE_MULTIPLE_EXACT_ALLELE_ROWS']
    assert problematic.empty,'Duplicate exact outcome rows require source investigation'
    print('PASS: 60 frozen exposure-branch/outcome pairs traced; no keep-first duplicate handling',flush=True)


if __name__=='__main__': main()
