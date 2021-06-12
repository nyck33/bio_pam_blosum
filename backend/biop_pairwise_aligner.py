from Bio import Align
from Bio.Align import substitution_matrices
from rq import Queue, Retry
from worker2 import conn
import time


# sample seqs
x = ("MATLKDQLIVNLLKEEQAPQNKITVVGVGAVGMACAISILMKDLADELALVDVMEDKLKGEMMDLQHGSL"
+"FLKTPKIVSSKDYCVTANSKLVIITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKYSPHCKLLIVSNPV"
+"DILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHALSCHGWVLGEHGDSSVPVWSGVNVAGVS"
+"LKSLNPELGTDADKEQWKEVHKQVVDSAYEVIKLKGYTSWAIGLSVADLAESIMKNLRRVHPISTMIKGL"
+"YGINEDVFLSVPCILGQNGISDVVKVTLTPEEEARLKKSADTLWGIQKELQF")

y = ("MATLKDQLIYNLLKEEQTPQNKITVVGVGAVGMACAISILMKDLADELALVDVIEDKLKGEMMDLQHGSL"
+"FLRTPKIVSGKVDILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHPLSCHGWVLGEHGDSSV"
+"PVWSGMNVAGVSLKTLHPDLGTDKDKEQWKEVHKQVVESAYEVIKLKGYTSWAIGLSVADLAESIMKNLR"
+"RVHPVSTMIKGLYGIKDDVFLSVPCILGQNGISDLVKVTLTSEEEARLKKSADTLWGIQKELQF")

x = "ACGTACGT"

y = "TACGTACG"

#aligner = Align.PairwiseAligner()


def global_align(aligner, seq1, seq2, matrix, gap_open=-10, gap_extend=-0.5):
    """

    :param aligner: aligner obj
    :param seq1:
    :param seq2:
    :param matrix:
    :return: alignments an iterator of alignment objects
    """
    # set params
    aligner.match_score = 1.0
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    aligner.mode = 'global'
    assert aligner.algorithm == 'Needleman-Wunsch', f'{aligner.algorithm}'
    #print('here')
    aligner.substitution_matrix = matrix
    alignments = aligner.align(seq1, seq2)
    # all alignments have same score
    score = alignments.score
    print(score)

    return alignments

def global_align2(seq1, seq2, matrix, gap_open=-10, gap_extend=-0.5):
    aligner = Align.PairwiseAligner()
    # set params
    aligner.match_score = 1.0
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    aligner.mode = 'global'
    #assert aligner.algorithm == 'Needleman-Wunsch', f'{aligner.algorithm}'
    aligner.substitution_matrix = matrix
    alignments = aligner.align(seq1, seq2)
    # all alignments have same score
    score = alignments.score
    #print(score)
    return alignments

def local_align(aligner, seq1, seq2, matrix, gap_open=-10, gap_extend=-0.5):
    """

    :param aligner: aligner obj
    :param seq1:
    :param seq2:
    :param matrix:
    :return: alignments an iterator of alignment objects
    """
    # set params
    aligner.match_score = 1.0
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    aligner.mode = 'local'
    assert aligner.algorithm == 'Smith-Waterman', f'{aligner.algorithm}'
    aligner.substitution_matrix = matrix
    alignments = aligner.align(seq1, seq2)
    # all alignments have same score
    score = alignments.score

    return alignments

def get_protein_alignment(alignments, idx):
    """

    :param alignments: iterator of alignment objects
    :param idx: choose alignments[idx]
    :return: seqA and seqB alignments to make aligned FASTA
    """
    target = str(alignments[idx])
    seq_arr = target.split()
    seqA, _, seqB = seq_arr

    return seqA, seqB

def main():
    q = Queue(connection=conn)
    mat_name = "BLOSUM62"
    matrix = substitution_matrices.load(mat_name)
    aligner = Align.PairwiseAligner()
    # aligner.substitution_matrix = matrix
    job = q.enqueue(global_align, args=(aligner, x, y, matrix))
    # alignments = global_align()
    count = 0
    while True:
        if job.result != None or count > 100000:
            break
        time.sleep(2)
        count += 1
        print(f'job.get_id(): {job.get_id()}, '
              f'job.result:{job.result}')
    alignments = job.result
    print(f'alignments[0]:{alignments[0]}\n score: {alignments[0].score}')

def tester():
    main()

if __name__ == "__main__":
    tester()