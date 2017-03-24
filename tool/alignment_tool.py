from Bio import SeqIO
from Bio.pairwise2 import align, format_alignment


def read_pw_seq(fasta_file):
    fasta = SeqIO.parse(fasta_file, 'fasta')
    sequences = [record.seq for record in fasta]

    if len(sequences) != 2:
        raise ValueError("Input fasta file must contain exactly 2 sequences!")

    return sequences


def pw_aln(seq_a, seq_b, alg, match, mismatch, gap_open, gap_extend, score_only):
    if alg not in ['global', 'local']:
        raise ValueError("'{}' algorithm is not supported. Choose between 'global' and 'local'.")

    if alg == 'global':
        aln_func = align.globalms
    else:
        aln_func = align.localms

    if score_only:
        best_score = aln_func(seq_a, seq_b, match, mismatch, gap_open, gap_extend, score_only=True)
        print(best_score)
    else:
        alignments = aln_func(seq_a, seq_b, match, mismatch, gap_open, gap_extend, score_only=False)
        for aln in alignments:
            print(format_alignment(*aln))

if __name__ == '__main__':
    seq = read_pw_seq('cluster0_K562CtcfRep1_chr1_839841_856780.fasta')
    pw_aln(*seq, alg='local', match=1, mismatch=2, gap_open=1, gap_extend=-1, score_only=False)
