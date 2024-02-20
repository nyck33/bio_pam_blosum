"""
Pseudcode for Traceback to get one possible alignment
"""
import numpy as np
from copy import deepcopy 

# get these from previous function
scoreMat=np.zeros((1,1)) 
pointerMat=np.zeros((1,1))

#coordinates_arr_of_arrs = []

rows, cols = np.shape(scoreMat)

# start at bottom right
row = rows-1
col = cols-1

#first indices in arr of arrs of tuples
#coord_arr = [[(row,col)]]

#initialize seq's and gap
seqA = ""
seqB = ""
#seqs_arr = []
gap = "-"

# dict from code to alpha residue since numpy array doesn't hold alphas
code_to_seqs_dict={"code": "residue"}

seqs_dict = {"diag": [], "up": "", "left": ""}
while not done:
    """
    at bottom right look at idx 0,1,2 of PM[i,j]
    """
    #next_coord_arr = []
    #for z in range(len(coord_arr)):
    codeA = scoreMat[row,0]
    codeB = scoreMat[0, col]
    #translate with dict
    residueA = code_to_seqs_dict[str(codeA)]
    residueB = code_to_seqs_dict[str(codeB)]
    # check diag
    if pointerMat[row, col, 0] == 1:
        seqA += residueA
        seqB += residueB
        #next_coord_arr = deepcopy(curr_coord_arr)
        #next_coord_arr.append((row-1, col-1))
        row-=1
        col-=1
        break
    #check left
    elif pointerMat[row, col-1, 1] == 1:
        seqA += gap
        seqB += residueB
        #next_coord_arr.append((row, col-1))
        col-=1
        break
    #check up
    elif pointerMat[row-1, col, 2] == 1:
        seqA += residueA
        seqB += gap
        #next_coord_arr.append((row-1, col))
        row-=1
        break
    
    if row==0 and col==0:
        done=True

print(f'seqA:\n{seqA}')
print(f'seqB:\n{seqB}')

