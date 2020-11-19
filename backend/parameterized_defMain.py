'''
All the other functions are okay but I want to call function main with parameters
for (sequence 1, sequence 2, match score, mismatch score, gap penalty, matrix)
Pleae refer below
'''

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
if __name__ == '__main__':
	seq1 = "ATGGC"
	seq2 = "ACTG"

	match_score = 5
	mismatch_score = -5
	gap_penalty = -5
	matrix = "BLOSUM 62"
	main(seq1, seq2, match_score, mismatch_score, gap_penalty, matrix)