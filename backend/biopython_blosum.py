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
matrix=substitution_matrices.load("BLOSUM62")
seq1 = SeqIO.read("example_fasta_files/homo_sapiens_lactate.fasta", "fasta")
seq2 = SeqIO.read("example_fasta_files/mus_musculus_lactate.fasta", "fasta")

#params: match, mismatch, gap open, gap extend
match = 2
mismatch = -1
gap_open = -10
gap_extend = -0.5
alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, matrix, 
                        gap_open, gap_extend)
print(type(alignments))
print(len(alignments))

x=pairwise2.format_alignment(*alignments[0])
print(type(x))
print(f'x:{x}')

###################################
#local alignment
print("local\n")
gap_extend = -1

alignments = pairwise2.align.localds(seq1, seq2, matrix, gap_open, gap_extend )

print(f'local len: {len(alignments)}')

y = pairwise2.format_alignment(*alignments[0], full_sequences=True)

print(f'y:\n{y}')