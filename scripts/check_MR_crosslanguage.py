import json
import math
import pandas as pd
from harmonise_core import ivw_crosscheck
from build_instruments import ROUND

fixtures = {
    'underdispersed': ([.1,.2,.3,.4],[.0301,.0599,.0901,.1199],[.05]*4),
    'heterogeneous': ([.1,.2,.3,.4],[.3,-.1,.4,.05],[.02,.03,.04,.02]),
    'two_snps': ([.1,.3],[.03,.05],[.03,.04]),
    'signed': ([-.1,.2,-.3,.4],[-.0301,.0599,-.0901,.1199],[.05]*4),
}
results = []
for row in pd.read_csv(ROUND/'tables/synthetic_MR_R_crosscheck.csv').to_dict('records'):
    expected = ivw_crosscheck(*fixtures[row['fixture']])
    for col,value in expected.items():
        assert math.isclose(row[col],value,rel_tol=1e-10,abs_tol=1e-12), (row['fixture'],col,row[col],value)
    results.append(dict(fixture=row['fixture'],status='PASS'))
(ROUND/'reports/MR_crosslanguage_validation.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print('PASS: four synthetic datasets, six IVW statistics each; Python independent formula agrees with installed TwoSampleMR.')
