#!/usr/bin/env python3

import numpy as np
from Bio import SeqIO
import os


class smith_waterman():
	def __init__(self, seq1, seq2, matrix_name="", match=2, mismatch=-3, gap_open=-10, gap_extend=-0.5):
		self.matrix_name = matrix_name
		self.matrix = self.load_matrix()
		#add space before each sequence
		self.seq1=" "+seq1[:]
		self.seq2=" "+seq2[:]
		self.gap_open=-gap_open
		self.gap_extend=-gap_extend
		self.match = match
		self.mismatch = mismatch

		#self.H=[]
		#self.S=[]
		#self.E=[]
		#self.F=[]

	#Load the subsititution matrix to the self.matrix
	def load_matrix(self):
		matrix = []
		if self.matrix_name:
			script_dir = os.path.dirname(__file__)
			rel_path ='matrices/'
			matrix_file_path = rel_path + self.matrix_name
			abs_file_path = os.path.join(script_dir, matrix_file_path)

			f = open(abs_file_path, 'r')
			lines = f.readlines()
			f.close()
			for line in lines:
				if "#" not in line:
					line=line.strip("\n")
					matrix.append(line.split())

			return matrix

	#if two inputs are matched, return match score, else return mismatch score
	def compare(self, m, n):
		if self.matrix_name == "":
			if m == n:
				return self.match
			else:
				return self.mismatch
		else:
			index_a = self.matrix[0].index(m)
			index_b = self.matrix[0].index(n)
			return int(self.matrix[index_a+1][index_b+1])



	def biuld_matrix(self):
		a = len(self.seq1)-1
		b = len(self.seq2)-1
		#initiate main score matrix
		self.S = np.zeros((a + 1, b + 1))
		#initiate scoring matrix with no gap
		self.H = np.zeros((a + 1, b + 1))
		#initialize matrix that second sequence has gap
		self.E = np.zeros((a + 1, b + 1))
		#initialize matrix that first sequence has gap
		self.F = np.zeros((a + 1, b + 1))
		#print(self.S)

		#adding space before each sequence.



		for i in range(1, b+1 if a>b else a+1):
			for j in range(i, b+1):
				#H stores score that two input is a match
				self.H[i,j] = self.S[i-1,j-1]+self.compare(self.seq1[i],self.seq2[j])
				#E stores score that one sequence has gap
				self.E[i,j] = max(np.add(self.S[0:i,j],-(np.arange(i,0,-1)*self.gap_extend+self.gap_open)))
				#E stores score that another sequence has gap
				self.F[i,j] = max(np.add(self.S[i,0:j],-(np.arange(j,0,-1)*self.gap_extend+self.gap_open)))
				#set the highest score for current position, if lowest score is negative, set to 0
				self.S[i,j] = max([0,self.H[i,j],self.E[i,j],self.F[i,j]])
			for j in range(i+1, a+1):
				self.H[j, i] = self.S[j-1, i-1]+self.compare(self.seq1[j], self.seq2[i])
				self.E[j, i] = max(np.add(self.S[0:j, i], -(np.arange(j,0,-1) * self.gap_extend + self.gap_open)))
				self.F[j, i] = max(np.add(self.S[j, 0:i], -(np.arange(i,0,-1) * self.gap_extend + self.gap_open)))
				self.S[j, i] = max([0, self.H[j, i], self.E[j, i], self.F[j, i]])

	

	def traceback(self, i, j, S1, S2, t1, t2):
		#print(t1, " ", t2, " ", i, " ", j)
		#if it go back to start, append the sequence
		if i == 1 and j == 1:
			S1.append(self.seq1[1] + t1[:])
			S2.append(self.seq2[1] + t2[:])
			return
		#if the current score is from no gap scoring matrix
		#go back diagonoly
		if self.S[i,j] == self.H[i,j]:
			self.traceback(i - 1, j - 1, S1, S2, self.seq1[i]+t1[:], self.seq2[j]+t2[:])
		#else go up or left
		if self.S[i,j] == self.E[i,j]:
			step = np.add(self.S[0:i,j],-(np.arange(i,0,-1)*self.gap_extend+self.gap_open))
			step = list(step)
			step = i-step.index(self.E[i, j])
			self.traceback(i - step, j, S1, S2, self.seq1[i]+t1[:], "-"*step+t2[:])
		if self.S[i,j] == self.F[i,j]:
			step = np.add(self.S[i, 0:j], -(np.arange(j, 0, -1)*self.gap_extend + self.gap_open))
			step = list(step)
			step = j - step.index(self.F[i, j])
			self.traceback( i, j-step, S1, S2, "-"*step+t1[:], self.seq2[j]+t2[:])




seq1 = SeqIO.read("example_fasta_files/homo_sapiens_lactate.fasta", "fasta")
seq2 = SeqIO.read("example_fasta_files/mus_musculus_lactate.fasta", "fasta")

a_str = str(seq1.seq)
b_str = str(seq2.seq)
SW = smith_waterman(seq1=seq1, seq2=seq2, matrix_name="BLOSUM62")
#print(SW.matrix)
SW.biuld_matrix()
R = np.where(SW.S == np.max(SW.S))
t1 = ""
t2 = ""
S1 = []
S2 = []
for i, j in zip(R[0],R[1]):
	
	SW.traceback(int(i), int(j), S1, S2, t1, t2)


print(S1)
print(S2)
#R = np.where(SW.S == np.max(SW.S))

# cast to string


#print("S = ",np.around(S,2))
#print("H = ",np.around(H,2))
#print("E = ",np.around(E,2))
#print("F = ",np.around(F,2))
 
'''S1 = []
S2 = []
t1 = ""
t2 = ""
a_str = " "+a_str[:]
b_str = " "+b_str[:]
R = np.where(S == np.max(S))

#demo strings
#a= " CTATAATCCC"
#b= " CTGTATC"

for x, y in zip(R[0],R[1]):
<<<<<<< HEAD
    traceback(a_str, b_str, 1, 1/3, S, x, y, S1, S2, t1, t2, H, E, F)
=======
    traceback(a_str, b_str, 1, 1, S, x, y, S1, S2, t1, t2, H, E, F)
>>>>>>> 5a4f9c91ddf568b299b719f41229632011944dca
print(S1)
print(S2)'''

