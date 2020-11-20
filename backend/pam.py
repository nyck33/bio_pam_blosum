#!/usr/bin/env python3

import numpy as np
import os


#return gap penalty or matrix score
def compare(a, b, matrix, gap):
    if " " == a or b == " ":
        return gap
    else:
        index_a = matrix[0].index(a)
        index_b = matrix[0].index(b)
    return int(matrix[index_a+1][index_b+1])

#return file path. if the input is PAM30, it will return the matrices/PAM30
def parse_name(mode):
    script_dir = os.path.dirname(__file__)
    rel_path ='matrices/'
    matrix_file_path = rel_path + mode
    abs_file_path = os.path.join(script_dir, matrix_file_path)

    return abs_file_path

#return the matrix based on the input file path.
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


#currently the matrix is looking for minimum number
#seq1, seq2 are two input sequences.
#matrix: input matrics
def build_matrics(seq1, seq2, matrix, gap):
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


#trace back still needs improve.
#seq1, seq2: input sequences
#matrix, read matrix from previous function
#gap: gap penalty
#S scoring matrix from previous function
#current_x, current_y: current matrix row and column
#S1, S2 alignment sequence
#t1, t2 temporary alignment sequence
#score: store score for alignment    ts: temporary score
def trace_back(seq1, seq2, matrix, gap, S, current_x, current_y, S1, S2, t1, t2, score, ts):
    #if edge reached, return
    if current_x == 0 and current_y == 0:   
        #print("reached here")
        S1.append(seq1[1] + t1[:])
        S2.append(seq2[1] + t2[:])
        score.append(ts)
        return
    #test if x and y reached the edge
    if current_x>0 and current_y>0:
        #test if current value is based on x-1 and y-1
        if S[current_x,current_y]-S[current_x-1,current_y-1] == compare(seq1[current_x], seq2[current_y], matrix, gap):
            #print(current_x, current_y, t1)
            ts = ts + S[current_x,current_y]
            trace_back(seq1, seq2, matrix, gap, S, current_x-1, current_y-1, S1, S2, seq1[current_x] + t1[:], seq2[current_y] + t2[:], score, ts)
    if current_y>0:
        #test if current value is based on x and y-1, if yes, there is gap on first sequence.
        if S[current_x,current_y]-S[current_x,current_y-1] == gap:
            #print(current_x, current_y, t1)
            ts = ts + S[current_x,current_y]
            trace_back(seq1, seq2, matrix, gap, S, current_x, current_y-1, S1, S2, "-" + t1[:], seq2[current_y] + t2[:], score, ts)
    if current_x>0:
        #test if current value is based on x-1 and y, if yes, there is a gap on second sequence
        if S[current_x,current_y]-S[current_x-1,current_y] == gap:
            #print(current_x, current_y, t1)
            ts = ts + S[current_x,current_y]
            trace_back(seq1, seq2, matrix, gap, S, current_x-1, current_y, S1, S2, seq1[current_x]+t1[:], "-"+t2[:], score, ts)
        
#seq1 seq2: input sequences
#mode is what of substitution matrix is going to use. ex "PAM30"
#gap: gap penalty scores.
def main():
    seq1 = "MATLKDQLIYNLLKEEQTPQ"
    seq2 = "NKITVVGVGAVGMACAISIL"
    mode="PAM30"
    file_name = parse_name(mode)
    matrix = load(file_name)
    S = build_matrics(seq1, seq2, matrix, -2)
    score = []
    S1 = []
    S2 = []
    t1 = ""
    t2 = ""
    ts = 0
    #seq1 = " " + seq1
    #seq2 = " " + seq2
    print(S)
    trace_back(" "+seq1[:], " "+seq2[:], matrix, -2, S, S.shape[0]-1, S.shape[1]-1, S1, S2, t1, t2, score, ts)
    print(max(score))
    print(S1[score.index(max(score))])
    print(S2[score.index(max(score))])
    #print(S)
    #return S1, S2, F


if __name__ == '__main__':
    main()


