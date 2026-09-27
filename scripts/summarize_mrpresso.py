from pathlib import Path
import json
import math
import re
import pandas as pd

ROUND=Path(__file__).resolve().parents[1]
ROOT=ROUND/'results/mrpresso'


def parse_p(value):
    text=str(value).strip()
    bound=text.startswith('<')
    number=float(text[1:] if bound else text)
    return number,bound


def main():
    seeds=[]; outliers=[]
    dirs=sorted(p for p in ROOT.iterdir() if p.is_dir())
    assert len(dirs)==24, f'Expected 24 completed jobs, found {len(dirs)}'
    for folder in dirs:
        global_row=pd.read_csv(folder/'global_test.csv').iloc[0].to_dict()
        p,p_bound=parse_p(global_row['Pvalue'])
        main=pd.read_csv(folder/'main_results.csv')
        warnings=(folder/'warnings.txt').read_text(encoding='utf-8').strip()
        raw=main[main['MR Analysis']=='Raw'].iloc[0]
        corrected=main[main['MR Analysis']=='Outlier-corrected']
        row=dict(exposure=global_row['exposure'],gwas_id=global_row['gwas_id'],outcome=global_row['outcome'],
            seed=int(global_row['seed']),n_sim=int(global_row['n_sim']),nsnp_input=int(global_row['nsnp_input']),
            global_RSS=float(global_row['RSSobs']),global_p_reported=str(global_row['Pvalue']),global_p_numeric_bound=p,
            global_p_is_upper_bound=p_bound,global_significant=(p<.05),raw_beta=float(raw['Causal Estimate']),
            raw_se=float(raw['Sd']),raw_p=float(raw['P-value']),warnings=warnings)
        if len(corrected):
            c=corrected.iloc[0]
            row.update(corrected_beta=float(c['Causal Estimate']),corrected_se=float(c['Sd']),corrected_p=float(c['P-value']))
        distortion=folder/'distortion_test.csv'
        if distortion.exists():
            d=pd.read_csv(distortion).iloc[0]
            row.update(outlier_markers=str(d['outlier_markers']),distortion_coefficient=float(d['distortion_coefficient']),
                       distortion_p_reported=str(d['distortion_p']))
        else:
            row.update(outlier_markers='',distortion_coefficient=math.nan,distortion_p_reported='')
        test=folder/'outlier_test.csv'
        if test.exists():
            t=pd.read_csv(test,dtype={'Pvalue':str})
            for r in t.itertuples():
                op,ob=parse_p(r.Pvalue)
                outliers.append(dict(**{k:row[k] for k in ['exposure','gwas_id','outcome','seed','n_sim']},
                    MarkerName=r.MarkerName,RSSobs=r.RSSobs,p_reported=r.Pvalue,p_numeric_bound=op,p_is_upper_bound=ob,
                    outlier_at_0_05=op<.05))
        seeds.append(row)
    seed_frame=pd.DataFrame(seeds)
    seed_frame.to_csv(ROUND/'tables/mrpresso_seed_results.csv',index=False)
    pd.DataFrame(outliers).to_csv(ROUND/'tables/mrpresso_outlier_tests.csv.gz',index=False,compression='gzip')
    pairs=[]
    for keys,d in seed_frame.groupby(['exposure','gwas_id','outcome'],sort=True):
        assert len(d)==2 and set(d.seed)=={271828,314159}
        sets=[set(filter(None,str(x).split(';'))) for x in d.outlier_markers]
        pair=dict(exposure=keys[0],gwas_id=keys[1],outcome=keys[2],n_sim_each=int(d.n_sim.iloc[0]),
            global_significant_both=bool(d.global_significant.all()),global_classification_consistent=d.global_significant.nunique()==1,
            outlier_sets_identical=sets[0]==sets[1],outlier_union=';'.join(sorted(sets[0]|sets[1])),
            outlier_intersection=';'.join(sorted(sets[0]&sets[1])),n_outlier_union=len(sets[0]|sets[1]),
            corrected_beta_range=float(d.corrected_beta.max()-d.corrected_beta.min()) if d.corrected_beta.notna().all() else math.nan,
            any_warning=bool(d.warnings.astype(bool).any()))
        pair['needs_higher_precision']=not (pair['global_classification_consistent'] and pair['outlier_sets_identical'])
        pairs.append(pair)
    pair_frame=pd.DataFrame(pairs)
    pair_frame.to_csv(ROUND/'tables/mrpresso_seed_stability.csv',index=False)
    report=dict(jobs=len(seed_frame),pairs=len(pair_frame),all_global_classification_consistent=bool(pair_frame.global_classification_consistent.all()),
        all_outlier_sets_identical=bool(pair_frame.outlier_sets_identical.all()),pairs_needing_higher_precision=pair_frame.loc[pair_frame.needs_higher_precision,['exposure','outcome']].to_dict('records'),
        caution='MR-PRESSO simulation p values retain inequality bounds and finite resolution; non-significance is not proof of no horizontal pleiotropy.')
    (ROUND/'reports/mrpresso_seed_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report),flush=True)


if __name__=='__main__': main()
