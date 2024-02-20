"""
Alternative Emboss suite in Biopython:
6.5.5  EMBOSS needle and water

The EMBOSS suite includes the water and needle tools for
Smith-Waterman algorithm local alignment, and Needleman-Wunsch
global alignment. The tools share the same style interface,
so switching between the two is trivial – we’ll just use needle here.

############################################################################
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
import os

class Needleman():
    def __init__(self, seq1, seq2, matrix_name="BLOSUM62", gap_open=-10, gap_extend=-0.5):  #rel_path='sub_matrix/'):
        self.matrix_name = matrix_name
        self.matrix = None
        self.seq1 = seq1
        self.seq2 = seq2
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.alignments = []
        #self.rel_path = rel_path

    def load_matrix(self):
        name_string = self.matrix_name
        matrix=substitution_matrices.load(name_string)
        self.matrix = matrix

    def align_global(self):
        #todo: cast to string
        seq1 = str(self.seq1)
        seq2 = str(self.seq2)
        # alignments is a list of alignments named tuples

        alignments = pairwise2.align.globalds(seq1, seq2, self.matrix,
                                              self.gap_open, self.gap_extend,
                                              penalize_end_gaps=False,
                                              one_alignment_only=True)

        self.alignments = alignments

    def align_local(self):
        seq1 = self.seq1
        seq2 = self.seq2
        # alignments is a list of alignments
        alignments = pairwise2.align.localds(seq1, seq2, self.matrix,
                                              self.gap_open, self.gap_extend,
                                             one_alignment_only=True)

        self.alignments = alignments

#todo: Heroku requires workers so take functions out of class
def matrix_load(mat_name="BLOSUM62"):
    mat = substitution_matrices.load(mat_name)
    return mat

def top_level_global(pairwise_align_ojb, seq1, seq2, matrix):
    gap_open = -10
    gap_extend = -0.5
    alignments = pairwise_align_ojb.globalds(seq1, seq2, matrix,
                                             gap_open, gap_extend,
                                             penalize_end_gaps=False,
                                             one_alignment_only=True)
    return alignments




def global_align_biop(seq1, seq2, matrix):
    gap_open = -10
    gap_extend = -0.5
    #seq1, seq2, matrix = args_tup
    alignments = pairwise2.align.globalds(seq1, seq2, matrix,
                                          gap_open, gap_extend,
                                          penalize_end_gaps=False,
                                          one_alignment_only=True)
    print(len(alignments))
    return alignments

def local_align_biop(seq1, seq2, matrix):
    gap_open = -10
    gap_extend = -0.5
    #, seq2, matrix = args_tup
    alignments = pairwise2.align.localds(seq1, seq2, matrix,
                                         gap_open, gap_extend,
                                         one_alignment_only=True)
    return alignments


homosapiens_lactate = ("MATLKDQLIVNLLKEEQAPQNKITVVGVGAVGMACAISILMKDLADELALVDVMEDKLKGEMMDLQHGSL"
+"FLKTPKIVSSKDYCVTANSKLVIITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKYSPHCKLLIVSNPV"
+"DILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHALSCHGWVLGEHGDSSVPVWSGVNVAGVS"
+"LKSLNPELGTDADKEQWKEVHKQVVDSAYEVIKLKGYTSWAIGLSVADLAESIMKNLRRVHPISTMIKGL"
+"YGINEDVFLSVPCILGQNGISDVVKVTLTPEEEARLKKSADTLWGIQKELQF")

mus_lactate = ("MATLKDQLIYNLLKEEQTPQNKITVVGVGAVGMACAISILMKDLADELALVDVIEDKLKGEMMDLQHGSL"
+"FLRTPKIVSGKVDILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHPLSCHGWVLGEHGDSSV"
+"PVWSGMNVAGVSLKTLHPDLGTDKDKEQWKEVHKQVVESAYEVIKLKGYTSWAIGLSVADLAESIMKNLR"
+"RVHPVSTMIKGLYGIKDDVFLSVPCILGQNGISDLVKVTLTSEEEARLKKSADTLWGIQKELQF")

