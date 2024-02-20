"""
The match parameters are:

CODE  DESCRIPTION & OPTIONAL KEYWORDS
x     No parameters. Identical characters have score of 1, otherwise 0.
m     A match score is the score of identical chars, otherwise mismatch
      score. Keywords ``match``, ``mismatch``.
d     A dictionary returns the score of any pair of characters.
      Keyword ``match_dict``.
c     A callback function returns scores. Keyword ``match_fn``.

The gap penalty parameters are:

CODE  DESCRIPTION & OPTIONAL KEYWORDS
x     No gap penalties.
s     Same open and extend gap penalties for both sequences.
      Keywords ``open``, ``extend``.
d     The sequences have different open and extend gap penalties.
      Keywords ``openA``, ``extendA``, ``openB``, ``extendB``.
c     A callback function returns the gap penalties.
      Keywords ``gap_A_fn``, ``gap_B_fn``.

matrices: 
https://github.com/biopython/biopython/tree/master/Bio/Align/substitution_matrices/data

"""

from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices
#local alignment using Aligner:
from Bio import Align
import sys


matrix=substitution_matrices.load("BLOSUM62")
#seq1 = SeqIO.read("example_fasta_files/homo_sapiens_lactate.fasta", "fasta")
#seq2 = SeqIO.read("example_fasta_files/mus_musculus_lactate.fasta", "fasta")
# very long sequences of nucleotides
#seq1 = SeqIO.read("example_fasta_files/corvyllus.fasta", "fasta")
#seq2 = SeqIO.read("example_fasta_files/drosophilaAE014297.fasta", "fasta")
# very long sequences of proteins
seq1 = SeqIO.read("example_fasta_files/titinHomoSapiens.fasta", "fasta")
seq2 = SeqIO.read("example_fasta_files/unnamedHermetica.fasta", "fasta")
#params: match, mismatch, gap open, gap extend

seq_str1 = str(seq1.seq)
seq_str2 = str(seq2.seq)


sizeof1 = sys.getsizeof(seq_str1)
sizeof2 = sys.getsizeof(seq_str2)
print(type(seq_str1), type(seq_str2))
print(sizeof1, sizeof2)

print(f"size 1: {sys.getsizeof(seq_str1)}\n len: "
      f"{len(seq_str1)}")
print(f"size 2: {sys.getsizeof(seq_str2)}\n len: "
      f"{len(seq_str2)}")

match = 2
mismatch = -1
gap_open = -10
gap_extend = -0.5

'''
aligner = Align.PairwiseAligner()
aligner.open_gap_score = gap_open
aligner.extend_gap_score = gap_extend
aligner.mode= 'local'
for alignment in aligner.align(seq_str1, seq_str2):
    print("Score = %.1f:" % alignment.score)
    print(alignment[:200])
'''



####################################################3

#pairwise2 below for protein using sub_matrix ie. blosum62

alignments = pairwise2.align.globalds(seq_str1, seq_str2, matrix,
                        gap_open, gap_extend)

match_score = 5
mismatch_score = -4
#alignments = pairwise2.align.globalds(seq_str1, seq_str2, match_score, mismatch_score,
 #                       open=gap_open, extend=gap_extend)

#alignments = pairwise2.align.localxs(seq_str1, seq_str2,
 #                       open=gap_open, extend=gap_extend)
# check types of named tuples in list alignments
alignA = alignments[0]
seqA = alignA.seqA
seqB = alignA.seqB
score = alignA.score
start = alignA.start
end = alignA.endb

print(type(seqA), type(seqB), type(score), type(start), type(end))
print(seqA, seqB)
print(score, start,end)

print(type(alignments))
print(len(alignments))

x=pairwise2.format_alignment(*alignments[0])
print(type(x))
print(f'x:{x[:200]}')

###################################
'''
#local alignment

print("local\n")
gap_extend = -1

alignments = pairwise2.align.localds(seq1.seq, seq2.seq, matrix, gap_open, gap_extend)

print(f'local len: {len(alignments)}')

y = pairwise2.format_alignment(*alignments[0], full_sequences=True)

print(f'y:\n{y}')

'''