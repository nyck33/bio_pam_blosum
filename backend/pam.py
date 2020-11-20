#!/usr/bin/env python3

import numpy as np


def compare(a, b, matrix, gap):
	if " " == a or b == " ":
		return gap
	else:
		index_a = matrix[0].index(a)
		index_b = matrix[0].index(b)
		return int(matrix[index_a+1][index_b+1])


def parse_name(mode):
	file_prefix ='./sub_matrix/'
	#if mode =='PAM30':
	#	file_name = file_prefix+'PAM30'

	return file_prefix + mode


def load(file_name):
	matrix = []

	f = open(file_name, 'r')
	lines = f.readlines()
	f.close()
	for line in lines:
		if "#" not in line:
			line=line.strip("\n")
			matrix.append(line.split())

	return matrix



#currently the matrix is looking for minimum number, thus the scoreig
def build_matrics(seq1, seq2, matrix, gap=-1):
    a = len(seq1)
    b = len(seq2)

    #initialize the matrics
    S = np.zeros((a + 1, b + 1))
    
    #set 0 to the first row and the first column, or we can initialize the number from 1 to len(sequence)
    #as mentioned in https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    
    S[0, :] = np.fromfunction(lambda x, y: gap * (x + y), (1, b + 1), dtype=int)
    S[:, 0] = np.fromfunction(lambda x, y: gap * (x + y), (1,a + 1), dtype=int)
    
    #add empty space before the sequece.
    seq1 = " " + seq1[:]
    seq2 = " " + seq2[:]
    
    
    #longer length sequence set as row
    for i in range(1, b+1 if a>b else a+1):
        #iterate loop to find the maximum score.
        for j in range(i, b+1):
            last_score = [S[i - 1, j], S[i - 1, j - 1], S[i, j - 1]]
            change_score = [gap,compare(seq1[i],seq2[j], matrix, gap),gap]

            new_score = np.add(last_score, change_score)

            S[i,j] = max(new_score)
            # print(S[i,c])
        for j in range(i+1, a+1):

            last_score = [S[j-1,i], S[j-1, i-1], S[j, i - 1]]
            change_score = [gap,compare(seq1[j],seq2[i], matrix, gap),gap]
            new_score = np.add(last_score, change_score)
            S[j,i] = max(new_score)
    return S


def trace_back(seq1, seq2, matrix, gap, S, current_x, current_y, S1, S2, t1, t2):
    
    if current_x == 1 and current_y == 1:
        S1.append(seq1[1] + t1[:])
        S2.append(seq2[1] + t2[:])
        return
    if S[current_x,current_y]-S[current_x-1,current_y] == gap:
        trace_back(seq1, seq2, matrix, gap, S, current_x-1, current_y, S1, S2, seq1[current_x]+t1[:], "-"+t2[:])
    if S[current_x,current_y]-S[current_x,current_y-1] == gap:
        trace_back(seq1, seq2, matrix, gap, S, current_x, current_y-1, S1, S2, "-" + t1[:], seq2[current_y] + t2[:])
    if S[current_x,current_y]-S[current_x-1,current_y-1] == compare(seq1[current_x], seq2[current_y], matrix, gap):
        trace_back(seq1, seq2, matrix, gap, S, current_x-1, current_y-1, S1, S2, seq1[current_x] + t1[:], seq2[current_y] + t2[:])


def main(seq1, seq2, match_score, mismatch_score, gap_penalty):
	'''
	mode='PAM30'
	file_name = parse_name(mode)
	matrix = load(file_name)
	F = build_matrics("ATGGC","ACTG", matrix, -2)
	S1 = []
	S2 = []
	t1 = ""
	t2 = ""
	trace_back(" ATGGC", " ACTG", matrix, -2, F, F.shape[0]-1, F.shape[1]-1, S1, S2, t1, t2)
	print(S1)
	print(S2)
	print(F)
	'''
    #return S1, S2, F
if __name__ == '__main__':
	seq1 = "ATGGC"
	seq2 = "ACTG"

	match_score = 5
	mismatch_score = -5
	gap_penalty = -5
	matrix = "BLOSUM 62"
	S1, S2, F = main(seq1, seq2, match_score, mismatch_score, gap_penalty, matrix)