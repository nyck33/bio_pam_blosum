from Bio import Align
from Bio.Align import substitution_matrices
from rq import Queue, Retry
from worker2 import conn
import time
from biop_pairwise_aligner import (global_align, global_align2, x, y)


def tester():
    q = Queue(connection=conn)
    mat_name = "BLOSUM62"
    matrix = substitution_matrices.load(mat_name)
    #aligner = Align.PairwiseAligner()
    # aligner.substitution_matrix = matrix
    job = q.enqueue(global_align2, args=(x, y, matrix))
    # alignments = global_align()
    count = 0

    while True:
        if job.result is not None or count > 1000:
            print(f'result: {job.result}')
            break
        time.sleep(1)
        count += 1
        print(f'job.get_id(): {job.get_id()}\n'
              f'job result: {job.result}')

    alignments = job.result
    print(f'alignments[0]:{alignments[0]}\n score: {alignments[0].score}')

def tester2():
    q = Queue(connection=conn)
    mat_name = "BLOSUM62"
    matrix = substitution_matrices.load(mat_name)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    job = q.enqueue(aligner.align, args=(x, y))
    # alignments = global_align()
    count = 0
    '''
    while True:
        if job.result != None or count > 100000:
            break
        time.sleep(2)
        count += 1
        print(f'job.get_id(): {job.get_id()}, '
              f'job.result:{job.result}')
    '''
    alignments = job.result
    print(f'alignments[0]:{alignments[0]}\n score: {alignments[0].score}')


if __name__=="__main__":
    tester()