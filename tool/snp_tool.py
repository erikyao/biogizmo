from util.biomart_client import BiomartClient


def query(rsid_lst):
    with BiomartClient() as bm_client:
        result = bm_client.query_snp(rsid_lst, out_fmt='BED')
        if result:
            print(result)
        else:
            print("No such SNP found.")
