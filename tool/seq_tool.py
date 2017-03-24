from util.togows_client import TogowsClient


def hg_refseq(build, chrom, chrom_start, chrom_end):
    if build not in ['hg19', 'hg38']:
        raise ValueError("Genome build '{}' is not supported. Use 'hg19' or 'hg38'")

    with TogowsClient(build) as tg_client:
        print(tg_client.query(chrom, chrom_start, chrom_end))
