"""Strict, non-palindromic, GRCh38 SNV alignment after reciprocal liftover.

Outcome alleles are REF/ALT on the reference assembly. Do not rescue a mismatch
by arbitrary strand complementation; the chain conversion already did that.
"""
import math
import numpy as np
from scipy.stats import norm, chi2


def align(effect_allele38, other_allele38, ref, alt, beta, se, af_alt):
    alleles = [effect_allele38, other_allele38, ref, alt]
    if any(not isinstance(a, str) or a not in {'A','C','G','T'} for a in alleles):
        return {'status': 'EXCLUDE_NON_ACGT_SNV'}
    if effect_allele38 == other_allele38 or ref == alt:
        return {'status': 'EXCLUDE_IDENTICAL_ALLELES'}
    if set([effect_allele38, other_allele38]) in [{'A','T'}, {'C','G'}]:
        return {'status': 'EXCLUDE_PALINDROMIC_NO_EXPOSURE_EAF'}
    if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
        return {'status': 'EXCLUDE_INVALID_EFFECT_OR_SE'}
    frequency_missing = af_alt is None or math.isnan(af_alt)
    if not frequency_missing and (not math.isfinite(af_alt) or not 0 < af_alt < 1):
        return {'status': 'EXCLUDE_INVALID_OUTCOME_FREQUENCY'}
    if (effect_allele38, other_allele38) == (alt, ref):
        sign = 1
    elif (effect_allele38, other_allele38) == (ref, alt):
        sign = -1
    else:
        return {'status': 'EXCLUDE_REFERENCE_ALLELE_MISMATCH'}
    eaf = math.nan if frequency_missing else (af_alt if sign == 1 else 1-af_alt)
    return {'status': 'ALIGNED', 'beta_outcome': sign*beta, 'se_outcome': se,
            'eaf_outcome': eaf, 'beta_sign': sign,
            'frequency_status': 'MISSING_NOT_REQUIRED_FOR_NONPALINDROMIC_MR' if frequency_missing else 'OBSERVED_OUTCOME_FREQUENCY'}


def checked_unique(records, key):
    result = {}
    for row in records:
        ident = key(row)
        if ident in result:
            raise ValueError('Duplicate variant requires explicit source review; never keep first')
        result[ident] = row
    return result


def ivw_crosscheck(bx, by, sy):
    """Independent test oracle, not the production MR package."""
    bx, by, sy = (np.asarray(x, dtype=float) for x in [bx,by,sy])
    if len(bx) < 2 or len(bx) != len(by) or len(bx) != len(sy):
        raise ValueError('IVW requires >=2 aligned instruments')
    if not (np.isfinite(bx).all() and np.isfinite(by).all() and np.isfinite(sy).all() and (sy > 0).all()):
        raise ValueError('Nonfinite or invalid input')
    information = np.sum(bx**2/sy**2)
    if information <= 0:
        raise ValueError('No exposure information')
    beta = np.sum(bx*by/sy**2)/information
    q = np.sum((by-beta*bx)**2/sy**2)
    se = math.sqrt(max(1, q/(len(bx)-1))/information)
    return dict(b=float(beta), se=se, pval=float(2*norm.sf(abs(beta/se))),
                Q=float(q), Q_df=len(bx)-1, Q_pval=float(chi2.sf(q,len(bx)-1)))


def bh_with_missing(pvalues, family_size=12):
    if len(pvalues) != family_size:
        raise ValueError('Family must contain all predeclared hypotheses')
    raw = np.asarray(pvalues, dtype=float)
    observed = np.isfinite(raw)
    if np.any((raw[observed] < 0) | (raw[observed] > 1)):
        raise ValueError('Invalid p value')
    vals = np.where(observed, raw, 1.)
    order = np.argsort(vals, kind='stable')
    adjusted = np.minimum.accumulate((vals[order]*family_size/np.arange(1,family_size+1))[::-1])[::-1]
    out = np.empty(family_size)
    out[order] = np.minimum(1.,adjusted)
    out[~observed] = np.nan
    return out
