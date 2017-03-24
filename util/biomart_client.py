from biomart import BiomartServer
from StringIO import StringIO
import pandas as pd


class BiomartClient:
    def __init__(self):
        # GRCh37 is also known as hg19

        # server = BiomartServer("http://useast.ensembl.org/biomart")
        self.server = BiomartServer("http://grch37.ensembl.org/biomart")

        # set verbose to True to get some messages
        # server.verbose = True

        # server.show_databases()

        self.database = self.server.databases["ENSEMBL_MART_SNP"]

        # db.show_datasets()

        self.dataset = self.database.datasets["hsapiens_snp"]

        # dataset.show_filters()
        # dataset.show_attributes()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self.dataset
        del self.database
        del self.server

    def query_snp(self, rsid_seq, out_fmt='dfm'):
        if out_fmt.lower() not in ['dfm', 'bed']:
            raise ValueError('{} format is not supported. Use "dfm" or "BED".'.format(out_fmt))

        response = self.dataset.search({
            'filters': {
                'snp_filter': rsid_seq
            },
            'attributes': [
                'chr_name', 'chrom_start', 'chrom_end', 'refsnp_id', 'chrom_strand'
            ]
        }, header=1)

        lines = [line.decode('utf-8') for line in response.iter_lines()]

        tsv_str = StringIO("\r\n".join(lines))

        dfm = pd.read_csv(tsv_str, header=0, sep="\t",
                              names=["chrom", 'chromStart', 'chromEnd', "name", "strand"])

        # Change to BED style
        dfm.loc[:, 'chrom'] = 'chr' + dfm.loc[:, 'chrom'].astype(str)
        dfm.loc[:, 'strand'] = dfm.loc[:, 'strand'].replace({1: '+', -1: '-'})
        # Change from 0-based to 1-based
        dfm.loc[:, 'chromStart'] = dfm.loc[:, 'chromStart'] - 1

        if out_fmt.lower() == 'dfm':
            return dfm.loc[:, ["name", "chrom", 'chromStart', 'chromEnd', "strand"]]
        else:  # out_fmt.lower() == 'bed'
            dfm = dfm.assign(score='0')
            dfm = dfm.loc[:, ["chrom", 'chromStart', 'chromEnd', "name", "score", "strand"]]

            bed_str = StringIO()
            dfm.to_csv(bed_str, sep='\t', header=None, index=None)

            return bed_str.getvalue()
