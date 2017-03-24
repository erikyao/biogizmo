import argparse
from util.ucsc_gb_client import chia_pet_tables
from tool import snp_tool, chia_pet_tool, alignment_tool, seq_tool


def add_subparser_hg_refseq(_subparsers):
    # Human Genome Sequence sub-command
    parser_hg_refseq = _subparsers.add_parser("hg-refseq", add_help=True,
                                              description="Fetch Human Genome Reference Sequence.")
    parser_hg_refseq.add_argument('--build', dest='build', type=str, required=True, choices=['hg19', 'hg38'],
                                  help="human genome build")
    parser_hg_refseq.add_argument('--chrom', dest='chrom', type=str, required=True,
                                  help="chrom of interest (Must start with 'chr'. E.g. 'chr1')")
    parser_hg_refseq.add_argument('--chrom-start', dest='chrom_start', type=int, required=True,
                                  help="start position on chromosome (0-based)")
    parser_hg_refseq.add_argument('--chrom-end', dest='chrom_end', type=int, required=True,
                                  help="end position on chromosome (0-based)")


def add_subparser_snp(_subparsers):
    # SNP sub-command
    parser_snp = _subparsers.add_parser("snp", add_help=True,
                                        description="Fetch SNP information.")
    parser_snp.add_argument('-i', '--rsid', dest='rsid_lst', type=str, nargs='+', required=True,
                            help="rsid of interest; separate by space if there are multiple")


def add_subparser_chia_pet(_subparsers):
    # ChIA-PET sub-command
    parser_chia_pet = _subparsers.add_parser("chia-pet", add_help=True,
                                             description="Fetch ChIA-PET cluster sequences.")
    parser_chia_pet.add_argument('--table-key', dest='table_key', type=str, required=True,
                                 choices=chia_pet_tables.keys(),
                                 help="short name of the UCSC Genome Browser ChIA-PET tables")
    parser_chia_pet.add_argument('--scope-type', dest='scope_type', type=str, required=True,
                                 choices=['within', 'overlap'],
                                 help="when 'within', search ChIA-PET clusters within range ['start', 'end']; "
                                      "when 'overlap', search ChIA-PET clusters overlapping range ['start', 'end']")
    parser_chia_pet.add_argument('--chrom', dest='chrom', type=str, required=True,
                                 help="chrom of interest (Must start with 'chr'. E.g. 'chr1')")
    parser_chia_pet.add_argument('--chrom-start', dest='chrom_start', type=int, required=True,
                                 help="start position to search on chromosome (0-based)")
    parser_chia_pet.add_argument('--chrom-end', dest='chrom_end', type=int, required=True,
                                 help="end position to search on chromosome (0-based)")
    parser_chia_pet.add_argument('--out-dir', dest='out_dir', type=str, required=True,
                                 help="directory to for the output fasta files")


def add_subparser_pw_aln(_subparsers):
    def __penalty_type(x):
        x = float(x)
        if x > 0:
            raise argparse.ArgumentTypeError("A penalty cannot be positive. Your input is {}".format(x))
        return x

    def __score_type(x):
        x = float(x)
        if x <= 0:
            raise argparse.ArgumentTypeError("A score must be positive. Your input is {}".format(x))
        return x

    # Pairwise Alignment sub-command
    parser_pw_aln = _subparsers.add_parser("pw-aln", add_help=True,
                                           description="Perform pairwise alignment.")
    parser_pw_aln.add_argument('--alg', dest='alg', type=str, required=True,
                               choices=['global', 'local'],
                               help="whether to perform global alignment (Needleman-Wunsch) "
                                    "or local alignment (Smith-Waterman)")
    parser_pw_aln.add_argument('--input-fasta', dest='input_fasta', type=str, required=True,
                               help="the fasta file of input sequences (must contain exactly 2 sequences)")
    parser_pw_aln.add_argument('--match', dest='match', type=__score_type, required=True,
                               help="match score (must be positive)")
    parser_pw_aln.add_argument('--mismatch', dest='mismatch', type=__penalty_type, required=True,
                               help="mismatch penalty (cannot be positive)")
    parser_pw_aln.add_argument('--gap-open', dest='gap_open', type=__penalty_type, required=True,
                               help="penalty when opening a gap (cannot be positive)")
    parser_pw_aln.add_argument('--gap-extend', dest='gap_extend', type=__penalty_type, required=True,
                               help="penalty when extending an existing gap (cannot be positive)")
    parser_pw_aln.add_argument('--score-only', dest='score_only', action="store_true",
                               help="If specified, only the best alignment score will be reported.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Welcome to use biogizmo toolkit!", add_help=True)
    subparsers = parser.add_subparsers(dest="subparser_name")

    add_subparser_hg_refseq(subparsers)
    add_subparser_snp(subparsers)
    add_subparser_chia_pet(subparsers)
    add_subparser_pw_aln(subparsers)

    args = parser.parse_args()

    if args.subparser_name == 'hg-refseq':
        seq_tool.hg_refseq(args.build, args.chrom, args.chrom_start, args.chrom_end)
    elif args.subparser_name == 'snp':
        snp_tool.query(args.rsid_lst)
    elif args.subparser_name == 'chia-pet':
        chia_pet_df = chia_pet_tool.query(args.table_key, args.scope_type,
                                          args.chrom, args.chrom_start, args.chrom_end)
        chia_pet_tool.to_fasta(chia_pet_df, args.out_dir)
    elif args.subparser_name == 'pw-aln':
        seq_a, seq_b = alignment_tool.read_pw_seq(args.input_fasta)
        alignment_tool.pw_aln(seq_a, seq_b, args.alg, args.match, args.mismatch,
                              args.gap_open, args.gap_extend, args.score_only)
    else:
        # This should never happen!
        raise ValueError("'{}' cannot be recognized as a sub command.".format(args.subparser_name))

