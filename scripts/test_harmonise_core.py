import unittest
import math
import numpy as np
from harmonise_core import align, checked_unique, ivw_crosscheck, bh_with_missing


class HarmoniseTests(unittest.TestCase):
    def test_matching_alt(self):
        r = align('A','C','C','A',.2,.05,.1)
        self.assertEqual(r['beta_outcome'],.2)
        self.assertEqual(r['eaf_outcome'],.1)

    def test_matching_ref(self):
        r = align('A','C','A','C',.2,.05,.1)
        self.assertEqual(r['beta_outcome'],-.2)
        self.assertEqual(r['eaf_outcome'],.9)
        self.assertEqual(r['se_outcome'],.05)

    def test_no_unjustified_complement_rescue(self):
        self.assertEqual(align('A','C','T','G',.2,.05,.1)['status'], 'EXCLUDE_REFERENCE_ALLELE_MISMATCH')

    def test_palindromes_not_rescued_with_outcome_frequency(self):
        self.assertTrue(align('A','T','T','A',.2,.05,.01)['status'].startswith('EXCLUDE_PALINDROMIC'))

    def test_invalid_inputs(self):
        for beta,se,af in [(math.nan,.1,.1),(.1,0,.1),(.1,-.1,.1),(.1,.1,1),(.1,.1,math.inf)]:
            self.assertNotEqual(align('A','C','C','A',beta,se,af)['status'],'ALIGNED')
        self.assertNotEqual(align('A','AC','AC','A',.1,.1,.1)['status'],'ALIGNED')

    def test_missing_frequency_not_blanket_barrier_for_nonpalindromic_MR(self):
        for af in [None, math.nan]:
            result = align('A','C','C','A',.1,.1,af)
            self.assertEqual(result['status'],'ALIGNED')
            self.assertTrue(math.isnan(result['eaf_outcome']))

    def test_duplicates_fail_even_if_identical(self):
        with self.assertRaises(ValueError):
            checked_unique([{'id':'a'}, {'id':'a'}], lambda r:r['id'])

    def test_ivw_does_not_apply_to_one_snp(self):
        with self.assertRaises(ValueError):
            ivw_crosscheck([.1],[.2],[.1])

    def test_ivw_sign_invariance_for_both_effects(self):
        a = ivw_crosscheck([.1,.2,.3],[.03,.11,.09],[.01,.02,.03])
        b = ivw_crosscheck([-.1,.2,-.3],[-.03,.11,-.09],[.01,.02,.03])
        self.assertAlmostEqual(a['b'],b['b'])
        self.assertAlmostEqual(a['se'],b['se'])

    def test_missing_hypotheses_not_dropped_from_fdr(self):
        out = bh_with_missing([.01]+[math.nan]*11)
        self.assertAlmostEqual(out[0],.12)
        self.assertTrue(np.isnan(out[1:]).all())
        with self.assertRaises(ValueError):
            bh_with_missing([.01])


if __name__ == '__main__':
    unittest.main()
