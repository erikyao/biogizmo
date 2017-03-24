import unittest
from util.ucsc_gb_client import UcscGbClient


class UcscGbClientCase(unittest.TestCase):
    def test_query_chia_pet_cluster(self):
        with UcscGbClient() as gb:
            df = gb.query_chia_pet_cluster('K562CtcfRep1', scope_type='within',
                                           chrom='chr1', chrom_start=830000, chrom_end=860000).set_index('id')

            self.assertEqual(df.shape[0], 2)

            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b1_chrom'], 'chr1')
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b1_chromStart'], 839841)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b1_chromEnd'], 840733)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b2_chrom'], 'chr1')
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b2_chromStart'], 855600)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839841_856780', 'b2_chromEnd'], 856780)

            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b1_chrom'], 'chr1')
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b1_chromStart'], 839974)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b1_chromEnd'], 840594)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b2_chrom'], 'chr1')
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b2_chromStart'], 847929)
            self.assertEqual(df.loc['K562CtcfRep1_chr1_839974_848569', 'b2_chromEnd'], 848569)


if __name__ == '__main__':
    unittest.main()
