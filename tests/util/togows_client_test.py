import unittest
from util.togows_client import TogowsClient


class TogowsClientTestCase(unittest.TestCase):
    def test_query(self):
        with TogowsClient() as tg_client:
            result = tg_client.query('chr1', 872838, 872848)
            expected = ">hg19:chr1:872838-872848\n" \
                       "GGAGAGGGTAG\n"
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
