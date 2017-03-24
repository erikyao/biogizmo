import os
from util.togows_client import TogowsClient
from util.ucsc_gb_client import UcscGbClient


def query(table_key, scope_type, chrom, chrom_start, chrom_end):
    with UcscGbClient() as gb_client:
        return gb_client.query_chia_pet_cluster(table_key, scope_type, chrom, chrom_start, chrom_end)


def to_fasta(chia_pet_df, output_dir):
    if chia_pet_df.empty:
        print("No ChIA-PET cluster found.")
        pass

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def __to_fasta(row, _tg_client):
        b1_fasta = _tg_client.query(row['b1_chrom'], row['b1_chromStart'], row['b1_chromEnd'], row['id'] + '_b1|')
        b2_fasta = _tg_client.query(row['b2_chrom'], row['b2_chromStart'], row['b2_chromEnd'], row['id'] + '_b2|')

        fasta = b1_fasta + b2_fasta
        filename = os.path.join(output_dir, "cluster{}_{}.fasta".format(row.name, row['id']))
        print(filename)
        with open(filename, 'w') as f:
            f.write(fasta)

    with TogowsClient() as tg_client:
        chia_pet_df.apply(__to_fasta, axis=1, _tg_client=tg_client)

    print("{} cluster fasta files written to {}".format(chia_pet_df.shape[0], os.path.abspath(output_dir)))


# if __name__ == '__main__':
#     with UcscGbClient() as gb:
#         df = gb.query_chia_pet_cluster('K562CtcfRep1', scope_type='within',
#                                         chrom='chr1', chrom_start=830000, chrom_end=860000)
#         to_fasta(df, '.')


