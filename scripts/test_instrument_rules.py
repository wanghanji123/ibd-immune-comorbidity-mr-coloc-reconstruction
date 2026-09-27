import unittest
from build_instruments import eligible, merge_windows


class InstrumentRules(unittest.TestCase):
    def test_primary_does_not_relax(self):
        self.assertTrue(eligible(''))
        self.assertFalse(eligible('PAPER_COHORT_P_QC_RULE'))
        self.assertFalse(eligible('PALINDROMIC_NO_SOURCE_EAF'))

    def test_sensitivity_does_not_hide_other_failures(self):
        self.assertTrue(eligible('PAPER_COHORT_P_QC_RULE', True))
        self.assertFalse(eligible('PAPER_COHORT_P_QC_RULE;PALINDROMIC_NO_SOURCE_EAF', True))
        self.assertFalse(eligible('LD_ID_NOT_UNIQUE_OR_ABSENT', True))

    def test_regions_merge_transitively_and_not_across_chromosomes(self):
        rows = merge_windows([(1, 1, 'a'), (1, 900000, 'b'), (1, 1800000, 'c'), (2, 1, 'd')])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['start37'], 1)
        self.assertEqual(rows[0]['end37'], 2300000)
        self.assertEqual(rows[0]['index_sources'], 'a;b;c')


if __name__ == '__main__':
    unittest.main()
