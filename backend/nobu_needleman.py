"""
Each NCBI Blosum matrix is the same format
gap extend is set 5 to 10 times lower than gap open
"""

import numpy as np
import os

class Needleman():
    def __init__(self, mat_name="BLOSUM62", rel_path='sub_matrix/',
                    gap_open=-10, gap_extend=-2):
        # {amino alpha : amino code}, can give each amino in sequences a code
        self.mat_name = mat_name
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.amino_dict = {}
        self.matrix = np.zeros((1,1), dtype=int)
        self.rel_path = rel_path

    def parse_name(self,mode):
        script_dir = os.path.dirname(__file__)
        rel_path = self.rel_path
        matrix_file_path = rel_path + mode
        abs_file_path = os.path.join(script_dir, matrix_file_path)

        return abs_file_path

    def load(self,file_name):
        """
        :param file_name: full path filename of matrix
        :return: matrix initialized to zero with indices filled and
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
                #build amino dict with amino vals as 100 to 100 something
                if count == 0:
                    for i in range(len(line_split)):
                        amino_code = 100 + i
                        self.amino_dict[str(line_split[i])] = amino_code
                else:
                    # ignore alpha to make scores arr of arrs
                    line_no_alpha = line_split[1:]
                    score_lines_arr.append(line_no_alpha)

                count+=1

        return score_lines_arr

    def build_matrix(self, score_lines_arr):
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
                if i==0 and j>0:
                    for _, amino_code in amino_dict.items():
                        matrix[i,j] = amino_code
                # left col
                elif i >0 and j==0:
                    for _, amino_code in amino_dict.items():
                        matrix[i,j] = amino_code
                else: #row>0 and col>0
                    curr_line = score_lines_arr[count]
                    for k in range(curr_line):
                        matrix[i,j] = curr_line[k]

        self.matrix = np.copy(matrix)

    def initialize_row_col(self):
        """
        Put values in first row and col as multiple of gap penalty
        :return:
        """
        matrix = np.copy(self.matrix)
        gap = self.gap_open
        rows, cols = np.shape(matrix)
        matrix[1,1] = 0
        for i in range(rows):
            for j in range(cols):
                # top row of nums
                if i==1 and j>1:
                    matrix[i,j] = gap * (j-1)
                # left col of nums
                if i>1 and j==1:
                    matrix[i,j] = gap * (i-1)

        self.matrix = np.copy(matrix)

    def main(self,seq1, seq2, mode):
        file_name = self.parse_name(mode)
        score_lines_arr = self.load(file_name)
        self.build_matrix(score_lines_arr)

