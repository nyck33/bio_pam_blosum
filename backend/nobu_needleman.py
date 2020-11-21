"""
Each NCBI Blosum matrix is the same format
gap extend is set 5 to 10 times lower than gap open
author: Nobu Kim

>NP_034829.1 L-lactate dehydrogenase A chain isoform 1 [Mus musculus]
MATLKDQLIVNLLKEEQAPQNKITVVGVGAVGMACAISILMKDLADELALVDVMEDKLKGEMMDLQHGSL
FLKTPKIVSSKDYCVTANSKLVIITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKYSPHCKLLIVSNPV
DILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHALSCHGWVLGEHGDSSVPVWSGVNVAGVS
LKSLNPELGTDADKEQWKEVHKQVVDSAYEVIKLKGYTSWAIGLSVADLAESIMKNLRRVHPISTMIKGL
YGINEDVFLSVPCILGQNGISDVVKVTLTPEEEARLKKSADTLWGIQKELQF

>NP_001128711.1 L-lactate dehydrogenase A chain isoform 2 [Homo sapiens]
MATLKDQLIYNLLKEEQTPQNKITVVGVGAVGMACAISILMKDLADELALVDVIEDKLKGEMMDLQHGSL
FLRTPKIVSGKVDILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHPLSCHGWVLGEHGDSSV
PVWSGMNVAGVSLKTLHPDLGTDKDKEQWKEVHKQVVESAYEVIKLKGYTSWAIGLSVADLAESIMKNLR
RVHPVSTMIKGLYGIKDDVFLSVPCILGQNGISDLVKVTLTSEEEARLKKSADTLWGIQKELQF
"""

import numpy as np
import os
import copy

ex_seq1 = "MATLKDQLIVNLLKEEQAPQNKITVVGVGAVGMACAISILMKDLADELALVDVMEDKLKGEMMDLQHGSL"
+"FLKTPKIVSSKDYCVTANSKLVIITAGARQQEGESRLNLVQRNVNIFKFIIPNIVKYSPHCKLLIVSNPV"
+"DILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHALSCHGWVLGEHGDSSVPVWSGVNVAGVS"
+"LKSLNPELGTDADKEQWKEVHKQVVDSAYEVIKLKGYTSWAIGLSVADLAESIMKNLRRVHPISTMIKGL"
+"YGINEDVFLSVPCILGQNGISDVVKVTLTPEEEARLKKSADTLWGIQKELQF"

ex_seq2 = "MATLKDQLIYNLLKEEQTPQNKITVVGVGAVGMACAISILMKDLADELALVDVIEDKLKGEMMDLQHGSL"
+"FLRTPKIVSGKVDILTYVAWKISGFPKNRVIGSGCNLDSARFRYLMGERLGVHPLSCHGWVLGEHGDSSV"
+"PVWSGMNVAGVSLKTLHPDLGTDKDKEQWKEVHKQVVESAYEVIKLKGYTSWAIGLSVADLAESIMKNLR"
+"RVHPVSTMIKGLYGIKDDVFLSVPCILGQNGISDLVKVTLTSEEEARLKKSADTLWGIQKELQF"

class Needleman():
    def __init__(self, seq_a="", seq_b="FX", mat_name="BLOSUM62", rel_path='sub_matrix/',
                 gap=-5, gap_open=-10, gap_extend=-2):
        # {amino alpha : amino code}, can give each amino in sequences a code
        self.rel_path = rel_path
        self.seqA = seq_a
        self.seqB = seq_b
        self.mat_name = mat_name
        self.gap = gap
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.amino_dict = {}
        # the downloaded matrix
        self.matrix = np.zeros((1, 1), dtype=int)
        # the matrix with longer sequence in the top row
        self.trav_mat = np.zeros((1, 1), dtype=int)
        self.alignment_pos_arr = []



    def parse_name(self, mode):
        """
        return absolute file path of matrix file, ex. folder/blosum62
        :param mode:
        :return:
        """
        script_dir = os.path.dirname(__file__)
        rel_path = self.rel_path
        matrix_file_path = rel_path + mode
        abs_file_path = os.path.join(script_dir, matrix_file_path)

        return abs_file_path

    def load(self, file_name):
        """
        :param file_name: full path filename of matrix
        :return:
        the score lines as arr or arrs
        """
        # get all lines from matrix file
        f = open(file_name, 'r')
        lines = f.readlines()
        f.close()
        score_lines_arr = []
        count = 0

        # build aminos dict
        for line in lines:
            if "#" not in line:
                # split on whitespace
                line_split = line.split()
                # build amino dict with amino vals as 100 to 100 something
                if count == 0:
                    for i in range(len(line_split)):
                        amino_code = 100 + i
                        self.amino_dict[str(line_split[i])] = amino_code
                else:
                    # ignore alpha to make scores arr of arrs
                    line_no_alpha = line_split[1:]
                    score_lines_arr.append(line_no_alpha)

                count += 1

        return score_lines_arr

    def build_score_matrix(self, score_lines_arr):
        """
        use dict to fill top row and left col with amino codes
        use score lines arr to fill the rest

        :param score_lines_arr: list of lists of just the nums in blosum matrix
        :return:
        """
        # get copy of amino dict
        amino_dict = self.amino_dict
        # make dim+1 for amino indices
        dim = len(amino_dict.keys()) + 1
        # intialize matrix and get rows and cols
        matrix = np.zeros((dim, dim))
        rows, cols = np.shape(matrix)
        count = 0
        for i in range(rows):
            for j in range(cols):
                # top row amino indices
                if i == 0 and j > 0:
                    for _, amino_code in amino_dict.items():
                        matrix[i, j] = amino_code
                # left col
                elif i > 0 and j == 0:
                    for _, amino_code in amino_dict.items():
                        matrix[i, j] = amino_code
                # row>0 and col>0, so fill in num scores from downloaded matrix
                else:
                    curr_line = score_lines_arr[count]
                    for k in range(curr_line):
                        matrix[i, j] = curr_line[k]

        self.matrix = np.copy(matrix)

    def initialize_row_col(self):
        """
        Put values in first row and col as multiple of gap penalty
        :return:
        """
        trav_mat = np.copy(self.trav_mat)
        gap = self.gap
        rows, cols = np.shape(trav_mat)
        trav_mat[1, 1] = 0
        for i in range(rows):
            for j in range(cols):
                # top row of nums, starts at col 2
                if i == 1 and j > 1:
                    trav_mat[i, j] = gap * (j - 1)
                # left col of nums, start at row 2
                if i > 1 and j == 1:
                    trav_mat[i, j] = gap * (i - 1)

        self.trav_mat = np.copy(trav_mat)

    def build_seq_matrix(self):
        """
        build the matrix from 2 sequences to traverse
        :return:
        """
        seq_a = self.seqA
        seq_b = self.seqB
        num_rows = len(seq_a) + 1
        num_cols = len(seq_b) + 1
        trav_mat = np.zeros((num_rows, num_cols), dtype=int)
        amino_dict = self.amino_dict
        for i in range(num_rows):
            for j in range(num_cols):
                # top row is seqb
                if i == 0 and j > 0:
                    for amino in seq_b:
                        amino_code = amino_dict[amino]
                        trav_mat[i, j] = amino_code
                # left col
                if i > 0 and j == 0:
                    for amino in seq_a:
                        amino_code = amino_dict[amino]
                        trav_mat[i, j] = amino_code

        self.trav_mat = np.copy(trav_mat)

    def find_score_at_pos(self, trav_mat, row, col):
        """

        :param trav_mat:
        :param row:
        :param col:
        :return: score at pos [row, col]
        """
        amino_left = trav_mat[row, 0]
        amino_top = trav_mat[0, col]
        # index into blosum mat
        score_mat = self.matrix
        row_idx = np.argwhere(score_mat[:, 0] == amino_left)
        col_idx = np.argwhere(score_mat[0, :] == amino_top)
        score = score_mat[row_idx, col_idx]

        return score

    def find_max(self, trav_mat, row, col):
        """

        :param trav_mat:
        :param row:
        :param col:
        :return: return the max from left, up or nw
        """
        score = self.find_score_at_pos(trav_mat, row, col)
        # compare the three options
        gap = self.gap
        score_west = trav_mat[row, col - 1] + gap
        score_north = trav_mat[row - 1, col] + gap
        score_nw = trav_mat[row - 1, col - 1] + score

        scores_arr = [score_west, score_north, score_nw]
        high = max(scores_arr)

        return high

    def traverse(self):
        """
        dynamic programming function, calls find_max
        start at (2,2) with 4 options:
        gap in s1, gap in s2, match or mismatch
        look at Blosum62, find score, compare with gap penalty
        :param: score_lines_arr is list of lists,
        :return:
        """
        score_mat = np.copy(self.matrix)
        trav_mat = np.copy(self.trav_mat)
        rows, cols = np.shape(trav_mat)

        for i in range(2, rows, 1):
            for j in range(2, cols, 1):
                trav_mat[i, j] = self.find_max(trav_mat, i, j)

        self.trav_mat = np.copy(trav_mat)

    def check_OB(self, row, col):
        if row < 1 or col < 1:
            return True  # is out of bounds
        return False

    def traceback(self):
        """
        start at bottom right of trav_mat, check which paths north, west
        and nw are correct, make new lists, copy old path over and the lists
        reamining at the end are the full paths from bottom right to top left

        :return: list of lists of paths
        """
        score_mat = np.copy(self.matrix)
        trav_mat = np.copy(self.trav_mat)
        rows, cols = np.shape(trav_mat)
        start_row = rows - 1
        start_col = cols - 1
        gap = self.gap
        # visited matrix
        visited = np.zeros((rows, cols), dtype=int)

        # init arrays
        paths_arr = []
        curr_paths = [[(start_row, start_col)]]
        done_paths = []
        done = False
        while not done:
            if len(curr_paths) > 0:
                new_paths = []
                for path in curr_paths:
                    curr_cell = path[-1]
                    row, col = curr_cell
                    # if at (1,1) this path is done
                    if row == 1 and col == 1:
                        done_path = copy.deepcopy(path)
                        done_paths.append(done_path)
                        continue
                    # current cell value
                    curr_val = self.find_score_at_pos(trav_mat, row, col)
                    # current cell's match/mismatch score value
                    score = self.find_score_at_pos(score_mat, row, col)
                    # make new arrays for any valid paths
                    # new_path_nw, new_path_north, new_path_west = [], [], []
                    # check out-of-bounds
                    if not self.check_OB(row - 1, col - 1):
                        # get the values from the "previous" cell
                        nw_val = trav_mat[row - 1, col - 1]
                        if (nw_val + score) == curr_val:
                            # copy the old path and append the new cell
                            new_path_nw = copy.deepcopy(path)
                            new_path_nw.append((row - 1, col - 1))
                            new_paths.append(new_path_nw)
                    if not self.check_OB(row - 1, col):
                        north_val = trav_mat[row - 1, col]
                        if (north_val + gap) == curr_val:
                            new_path_north = copy.deepcopy(path)
                            new_path_north.append((row - 1, col))
                            new_paths.append(new_path_north)
                    if not self.check_OB(row, col - 1):
                        west_val = trav_mat[row, col - 1]
                        if (west_val + gap) == curr_val:
                            new_path_west = copy.deepcopy(path)
                            new_path_west.append((row, col - 1))
                            new_paths.append(new_path_west)

                # curr_paths has max 2 and min 1 len of valid paths
                curr_paths.clear()
                curr_paths = copy.deepcopy(new_paths)
                new_paths.clear()
            # no more paths
            else:
                done = True
        self.alignment_pos_arr = copy.deepcopy(done_paths)


    def alignments_to_strings(self):
        """
        Align shorter sequence against longer one, step into (1,1)
        starting at bottom right
        Output both aligned areas only and full alignments
        :return:
        """
        align_pos_arr = copy.deepcopy(self.alignment_pos_arr)


    def main(self):
        # get file path of matrix like folder/blosum62
        file_name = self.parse_name(self.mat_name)
        # get arr of arrs of score num lines in blosum62
        # and build self.amino_dict
        score_lines_arr = self.load(file_name)
        # build the score matrix
        self.build_score_matrix(score_lines_arr)
        # build seq_matrix with amino_codes in top row, left col
        self.build_seq_matrix()
        # initialize first row and col
        self.initialize_row_col()
        # traverse
        self.traverse()
        self.traceback()


