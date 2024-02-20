"""
Rice University resource
"""
import numpy

# The alphabet defines the order of amino acids used for the score matrix
# so A is the first character, R is the second etc.
alphabet = "ARNDCQEGHILKMFPSTWYV"

# The BLOSUM50 matrix in third-bit units, or 3log2(odds ratio)
# We store it as a 2D array using the numpy package
blosum50 = numpy.array([
    [5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0],
    [-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3],
    [-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3],
    [-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4],
    [-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1],
    [-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3],
    [-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3],
    [0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4],
    [-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4],
    [-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4],
    [-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1],
    [-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3],
    [-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1],
    [-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1],
    [-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3],
    [1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0],
    [-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3],
    [-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1],
    [0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5],
])

# Specify a linear GAP penalty of -8
gap_penalty = -8

# Our X and Y sequences to align
x = "MAMRLLKTHL"
y = "MKNITCYL"

#Homo sapiens lactate vs Mus musculus lacate
x = ("MATLKDQLIVNLLKEEQAPQNKITVVGVGAVGMACAISILMKDLADELALVDVMEDKLKGEMMDLQHGSL"
+"FLKTPKIVSSKDYCVTANSKLVIITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKYSPHCKLLIVSNPV"
+"DILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHALSCHGWVLGEHGDSSVPVWSGVNVAGVS"
+"LKSLNPELGTDADKEQWKEVHKQVVDSAYEVIKLKGYTSWAIGLSVADLAESIMKNLRRVHPISTMIKGL"
+"YGINEDVFLSVPCILGQNGISDVVKVTLTPEEEARLKKSADTLWGIQKELQF")

y = ("MATLKDQLIYNLLKEEQTPQNKITVVGVGAVGMACAISILMKDLADELALVDVIEDKLKGEMMDLQHGSL"
+"FLRTPKIVSGKVDILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHPLSCHGWVLGEHGDSSV"
+"PVWSGMNVAGVSLKTLHPDLGTDKDKEQWKEVHKQVVESAYEVIKLKGYTSWAIGLSVADLAESIMKNLR"
+"RVHPVSTMIKGLYGIKDDVFLSVPCILGQNGISDLVKVTLTSEEEARLKKSADTLWGIQKELQF")

titin = ""
with open("titinSequenceOnly.fasta", "r") as f:
    titin=f.read()

print(type(titin), titin)

y=titin

print(y)

# Set up our F-matrix using a two dimensional score_matrix to record the
# intermediate values and a three dimensional pointer_matrix to record the
# arrows. The first and second axes of the score_matrix and pointer_matrix
# correspond to the X and Y sequences respectively. The third axis of the
# pointer_matrix records the presence (1) or absence (0) of arrows in the
# diagonal (align both), left (align X with gap), and up (align Y with gap) at
# indices 0, 1 and 2 respectively. Because with have the first row and first
# column filled in with gap penalties, the row and column indices
# corresponding to the X and Y residues will begin at 1.
score_matrix = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
pointer_matrix = numpy.zeros((len(x) + 1, len(y) + 1, 3), dtype=int)

print(f'score matrix\n{score_matrix}')
print(f'pointer matrix\n{pointer_matrix}')
# Fill in first row and column with the linear gap penalties
for xi, xaa in enumerate(x):
    score_matrix[xi + 1, 0] = gap_penalty * (xi + 1)
    pointer_matrix[xi + 1, 0, 1] = 1

for yi, yaa in enumerate(y):
    score_matrix[0, yi + 1] = gap_penalty * (yi + 1)
    pointer_matrix[0, yi + 1, 2] = 1

# Fill in middle values starting, from left to right in the top row, then from
# left to right in the second row and so on. The yi and xi indices begin at 0.
for yi, yaa in enumerate(y):
    for xi, xaa in enumerate(x):
        # Get the index (position) of the X and Y amino acids in our alphabet
        xaai = alphabet.index(xaa)
        yaai = alphabet.index(yaa)

        # These are the indices to the X and Y positions in the F-matrix,
        # which begin at 1 because of the the first row and first column of
        # gap penalties.
        fxi = xi + 1
        fyi = yi + 1

        # Look up the score for the X and Y pair of residues in our BLOSUM50 matrix
        residue_score = blosum50[xaai, yaai]

        # Calculate the possible scores for aligning both residues...
        no_gap_score = score_matrix[
                           fxi - 1, fyi - 1] + residue_score  # (fxi - 1, fyi - 1) is the diagonally up and left cell
        # ...and for aligning the X residue with a gap...
        x_gap_score = score_matrix[
                          fxi - 1, fyi] + gap_penalty  # (fxi - 1, fyi) is the cell to the left of the current cell
        # ...and for aligning the Y residue with a gap.
        y_gap_score = score_matrix[fxi, fyi - 1] + gap_penalty  # (fxi, fyi - 1) is the cell above the current cell

        # Identify the maximum score value from among the three options...
        max_score = max(no_gap_score, x_gap_score, y_gap_score)
        # ...and set the score of this F-matrix cell to that score.
        score_matrix[fxi, fyi] = max_score

        # Add arrows for any direction where choosing that direction
        # will result in the maximum score for this F-matrix cell
        if no_gap_score == max_score:
            pointer_matrix[fxi, fyi, 0] = 1
        if x_gap_score == max_score:
            pointer_matrix[fxi, fyi, 1] = 1
        if y_gap_score == max_score:
            pointer_matrix[fxi, fyi, 2] = 1


# A recursive function to regenerate the optimal alignment from the F-matrix,
# randomly sampling if there is more than one optimal alignment
def recurse_score_matrix(xi, yi):
    # Once the top left corner, return two empty sequences (an empty pairwise alignment)
    if xi == yi == 0:
        return "", ""

    # Pick one of the outgoing arrows at random. Firt, generate a random
    # permutation of the sequence [0, 1, 2] corresponding to diagonal, left
    # and up arrows respectively.
    random_pointer_order = numpy.random.permutation(3)
    # Next, find the first number in that permutation where the arrow at that
    # index of the pointer array is present.
    for random_pointer_choice in random_pointer_order:
        if pointer_matrix[xi, yi, random_pointer_choice] == 1:
            break

    # Align two residues, or align x residue with a gap, or y residue with a gap.
    # xj, yj will refer to the indices of the
    if random_pointer_choice == 0:  # along both residues (diagonal arrow)
        xj = xi - 1  # xj is the index of the next cell (one index left)
        yj = yi - 1  # yj is the index of the next cell (one index up)
        # Get the X and Y sequence characters corresponding to the current
        # cell. Remember that the indices for sequence characters begin at 0
        # for the sequence strings but begin at 1 in the F-matrix, so we
        # subtract 1 to get the correct indices in the sequence strings while
        # we are traversing the F-matrix using xi, yi.
        xc = x[xi - 1]
        yc = y[yi - 1]
    elif random_pointer_choice == 1:  # align X with gap
        xj = xi - 1
        yj = yi
        xc = x[xi - 1]
        yc = "-"  # use a dash for the gap
    else:  # align Y with gap
        xj = xi
        yj = yi - 1
        xc = "-"
        yc = y[yi - 1]

    # Add the two characters (the two residues, or a residue and a gap) to the
    # partial alignment generated by recursing through the next cell (pointed
    # to by the randomly chosen arrow), and return the new partial alignment.
    partial_x, partial_y = recurse_score_matrix(xj, yj)

    partial_x += xc
    partial_y += yc

    return partial_x, partial_y


# Call the recursive function, beginning at the bottom right hand corner. The
# coordinate of the bottom right cell is equal to the length of the X sequence
# and the length of the Y sequence
optimal_alignments = recurse_score_matrix(len(x), len(y))

# Numpy treats the first axis as rows, and the second axis as columns. But we
# are treating the axes as X and Y respectively, so if we print the score matrix,
# X0 will be the first row, X1 will be the second etc. To make the X sequence
# run side-to-side we therefore transpose the matrix before printing.
print(score_matrix.transpose())
# Print the optimal alignments in separate lines (so they are aligned in the
# output).
print("\n" + optimal_alignments[0] + "\n" + optimal_alignments[1])
