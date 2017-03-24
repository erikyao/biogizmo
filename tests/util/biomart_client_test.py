import unittest
from util.biomart_client import BiomartClient


class BiomartClientCase(unittest.TestCase):
    def test_query_dfm(self):
        with BiomartClient() as bm_client:
            result = bm_client.query_snp(['rs10005089', 'rs1000237'], out_fmt='dfm').set_index('name')

            self.assertEqual(result.loc['rs1000237', 'chrom'], 'chr19')
            self.assertEqual(result.loc['rs1000237', 'chromStart'], 19518315)
            self.assertEqual(result.loc['rs1000237', 'chromEnd'], 19518316)
            self.assertEqual(result.loc['rs10005089', 'chrom'], 'chr4')
            self.assertEqual(result.loc['rs10005089', 'chromStart'], 75694605)
            self.assertEqual(result.loc['rs10005089', 'chromEnd'], 75694606)

    def test_query_bed(self):
        with BiomartClient() as bm_client:
            result = bm_client.query_snp(['rs10005089', 'rs1000237'], out_fmt='BED')

            expected = "chr19	19518315	19518316	rs1000237	0	+\n" \
                       "chr4	75694605	75694606	rs10005089	0	+\n"

            self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()
