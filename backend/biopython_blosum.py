from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices
blosum62=substitution_matrices.load("BLOSUM62")
seq1 = SeqIO.read("example_fasta_files/homo_sapiens_lactate.fasta", "fasta")
seq2 = SeqIO.read("example_fasta_files/mus_musculus_lactate.fasta", "fasta")
alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
print(len(alignments))

print(pairwise2.format_alignment(*alignments[0]))

