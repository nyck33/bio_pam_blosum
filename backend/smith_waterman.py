#!/usr/bin/env python3

import numpy as np

def compare(m, n, match, mismatch):
    if m == n:
        return match
    else:
        return mismatch

def biuld_matrix(seq1, seq2, match, mismatch, v, u):
	a = len(seq1)
	b = len(seq2)

	S = np.zeros((a + 1, b + 1))
	H = np.zeros((a + 1, b + 1))
	E = np.zeros((a + 1, b + 1))
	F = np.zeros((a + 1, b + 1))

	seq1 = " "+seq1[:]
	seq2 = " "+seq2[:]


	for i in range(1, b+1 if a>b else a+1):
		for j in range(i, b+1):
			H[i,j] = S[i-1,j-1]+compare(seq1[i],seq2[j],match,mismatch)
			E[i,j] = max(np.add(S[0:i,j],-(np.arange(i,0,-1)*u+v)))
			F[i,j] = max(np.add(S[i,0:j],-(np.arange(j,0,-1)*u+v)))
			S[i,j] = max([0,H[i,j],E[i,j],F[i,j]])
		for j in range(i+1, a+1):
			H[j, i] = S[j-1, i-1]+compare(seq1[j], seq2[i],match,mismatch)
			E[j, i] = max(np.add(S[0:j, i], -(np.arange(j,0,-1) * u + v)))
			F[j, i] = max(np.add(S[j, 0:i], -(np.arange(i,0,-1) * u + v)))
			S[j, i] = max([0, H[j, i], E[j, i], F[j, i]])

	return S, H, E, F

def traceback(seq1, seq2, v, u, S, i, j, S1, S2, t1, t2, H, E, F):
	print(t1, " ", t2, " ", i, " ", j)
	if i == 1 and j == 1:
		S1.append(seq1[1] + t1[:])
		S2.append(seq2[1] + t2[:])
		return
	if S[i,j] == H[i,j]:
		print("first")
		traceback(seq1, seq2, v, u, S, i - 1, j - 1, S1, S2, seq1[i] + t1[:], seq2[j] + t2[:], H, E, F)
	if S[i,j] == E[i,j]:
		step = np.add(S[0:i,j],-(np.arange(i,0,-1)*u+v))
		step = list(step)
		step = i-step.index(E[i, j])
		print("second", step, E[i,j])
		traceback(seq1, seq2, v, u, S, i - step, j, S1, S2, seq1[i]+t1[:], "-"*step+t2[:], H, E, F)
	if S[i,j] == F[i,j]:
		step = np.add(S[i, 0:j], -(np.arange(j, 0, -1) * u + v))
		step = list(step)
		step = j - step.index(F[i, j])
		print("third", step)
		traceback(seq1, seq2, v, u, S, i, j-step, S1, S2, "-"*step+t1[:], seq2[j]+t2[:], H, E, F)



S, H, E, F = biuld_matrix("CTATAATCCC", "CTGTATC", 1, -1, 1, 1)

print("S = ",np.around(S,2))
print("H = ",np.around(H,2))
print("E = ",np.around(E,2))
print("F = ",np.around(F,2))
 
S1 = []
S2 = []
t1 = ""
t2 = ""
R = np.where(S == np.max(S))
for x, y in zip(R[0],R[1]):
    traceback(" CTATAATCCC", " CTGTATC", 1, 1, S, x, y, S1, S2, t1, t2, H, E, F)
print(S1)
print(S2)

