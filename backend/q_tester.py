from Bio import Align
from Bio.Align import substitution_matrices
from rq import Queue, Retry
#from worker import conn
from worker2 import conn2
import time
from biop_pairwise_aligner import (global_align2, x, y,
                                   local_align)


def tester():
    q = Queue(connection=conn2)
    mat_name = "BLOSUM62"
    matrix = substitution_matrices.load(mat_name)
    #aligner = Align.PairwiseAligner()
    # aligner.substitution_matrix = matrix
    job = q.enqueue(local_align, args=(x, y, matrix))
    # alignments = global_align()
    count = 0

    while True:
        if job.result is not None or count > 100:
            #print(f'result: {job.result}')
            break
        time.sleep(1)
        count += 1
        print(f'job.get_id(): {job.get_id()}\n')
              #f'job result: {job.result}')

    #only returning one
    alignment = job.result
    seqA, connector, seqB = get_protein_alignment(alignment)
    #todo: this is just returning seqA
    seqA_adj = make_single_seq(seqA, connector)
    seqB_adj = make_single_seq(seqB, connector)

    #assert len(seqA_adj) == len(seqB_adj) == len(connector)

    print(f'strA:\n{str(seqA)}\n'
          f'strB:\n{str(seqB)}')

    print(f'adjA:\n {seqA_adj}\n, adjB:\n {seqB_adj}')
    #print(type(alignment))
    #print(f'alignments[0]:{alignment}\n score: {alignment.score}')
    # practice pretty printing
    count = 0
    while True:
        print(f'line: {count}')
        print(seqA[count*50: (count*50)+50])
        #print('\n')
        print(connector[count * 50: (count * 50) + 50])
        #print('\n')
        print(seqB[count * 50: (count * 50) + 50])
        #print('\n')
        if (count*50) + 50 > len(connector):
            break
        count += 1


def get_protein_alignment(alignment):
    """

    :param alignments: iterator of alignment objects
    :param idx: choose alignments[idx]
    :return: seqA and seqB alignments to make aligned FASTA
    """
    target = str(alignment)
    seq_arr = target.split()
    seqA, connectors, seqB = seq_arr

    return seqA, connectors, seqB


def make_single_seq(seq, connector):
    """
    If connector is "|" want the protein else "-" for gap
    :param seq:
    :param connector:
    :return:
    """
    seq_with_gaps = ""
    for i in range(len(connector)):
        if connector[i] != "|":
            seq_with_gaps += "-"
        else:
            seq_with_gaps += seq[i]

    return seq


def tester2():
    q = Queue(connection=conn2)
    mat_name = "BLOSUM62"
    matrix = substitution_matrices.load(mat_name)
    aligner = Align.PairwiseAligner()
    Align.PairwiseAligner.__module__ ="__init__"
    aligner.substitution_matrix = matrix
    job = q.enqueue(aligner.align, args=(x, y))
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

def tester3():
    pass

if __name__=="__main__":
    tester()