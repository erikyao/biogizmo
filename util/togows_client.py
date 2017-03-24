from urllib2 import Request, urlopen


class TogowsClient:
    def __init__(self, build="hg19"):
        self.build = build

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self

    def query(self, chrom, chrom_start, chrom_end, defline_prefix=None):
        url_template = 'http://togows.org/api/ucsc/{build}/{chrom}:{chromStart}-{chromEnd}.fasta'
        url = url_template.format(build=self.build, chrom=chrom, chromStart=chrom_start, chromEnd=chrom_end)

        req = Request(url)

        response = urlopen(req)
        page = response.read().decode("utf8")
        if defline_prefix:
            idx = page.find('>') + len('>')
            page = page[:idx] + defline_prefix + page[idx:]

        response.close()

        return page

if __name__ == '__main__':
    with TogowsClient() as tg_client:
        print(tg_client.query('chr1', 872838, 872848))