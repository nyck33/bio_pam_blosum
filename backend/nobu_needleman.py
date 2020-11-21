"""
Each NCBI Blosum matrix is the same format
"""

import numpy as np
import os

class Needleman():
    def __init__(self):
        # {amino alpha : amino code}, can give each amino in sequences a code
        self.protein_dict = {}


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

        # get num alphas in line 0 and initialize matrix
        amino_line = lines[0]
        #need extra row and col for aminos indexes
        dim = len(amino_line) + 1
        matrix = np.zeros((dim,dim), dtype=int)
        rows, cols = np.shape(matrix)

        # build aminos dict
        for line in lines:
            if "#" not in line:
                # split on whitespace
                line_split = line.split()
                #build amino dict with amino vals as 100 to 100 something
                if count == 0:
                    for i in range(len(line_split)):
                        amino_code = 100 + i
                        self.protein_dict[str(line_split[i])] = amino_code
                        # leave zero for rows index
                        matrix[count, i+1] = amino_code
                else:
                    # from row 1 onwards
                    amino_code = 100 + (count-1)
                    matrix[count,0] = amino_code
                    # ignore alpha to make scores arr of arrs
                    line_no_alpha = line_split[1:]
                    assert len(line_no_alpha) == dim-1, "score line len and dim-1 mismatch"
                    score_lines_arr.append(line_no_alpha)

                count+=1

        return matrix, score_lines_arr

    def make_matrix(self, matrix, score_lines):
        """

        :param matrix:
        :param score_lines:
        :return:
        """
        # fill left index col
        rows, cols = dim, dim
        for i in range(rows):
            if i > 0:
                matrix[i, 0] =


    def parse_name(self,mode):
        script_dir = os.path.dirname(__file__)
        rel_path ='sub_matrix/'
        matrix_file_path = rel_path + mode
        abs_file_path = os.path.join(script_dir, matrix_file_path)

        return abs_file_path


    def main(self,seq1, seq2, mode, gap=-2):
        file_name = self.parse_name(mode)
        matrix = self.load(file_name)