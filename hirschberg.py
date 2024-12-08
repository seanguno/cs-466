def needleman_score_linear(seq_1, seq_2):
	# init dp table
	score = [[0] * 2 for i in range(len(seq_1) + 1)]
	# init first col
	for i in range(1, len(seq_1) + 1):
		# in future could include diff penalty dictionaries 
		score[i][0] = score[i - 1][0] - 1
	print(f'initial column: {[row[0] for row in score]}')
	# fill matrix col by col
	for j in range(1, len(seq_2) + 1):
		# init first row of curr col
		score[0][1] = score[0][0] - 1
		for i in range(1, len(seq_1) + 1):
			ins_score = score[i][0] - 1
			del_score = score[i - 1][1] - 1
			match_score = score[i - 1][0] + (1 if seq_1[i - 1] == seq_2[j - 1] else -1)
			# update based on max scores
			score[i][1] = max(ins_score, del_score, match_score)
		# update two column array
		for i in range(len(seq_1) + 1):
			score[i][0] = score[i][1]
		print(f'col {j}: {[row[1] for row in score]}')
	# get last row of score
	last_col = [row[1] for row in score]
	print(f"Last col for '{seq_1}' vs '{seq_2}': {last_col}")
	return last_col

def hirschberg_helper(seq_1, seq_2, i, j, i_p, j_p, result):
	if j_p - j == 1:
		result.append((i, j))
		result.append((i_p, j_p))
		print(f"Base case 1 result: {result}")
		return
	if j_p < 1:
		print(f"Base case 2 (j_p < 1): No operation.")
		return
	mid_j = int((j + j_p) / 2)
	prefix = needleman_score_linear(seq_1[i:i_p], seq_2[j:mid_j])
	suffix = needleman_score_linear(seq_1[i:i_p][::-1], seq_2[mid_j:j_p][::-1])
	print(f"Prefix scores: {prefix}")
	print(f"Suffix scores (reversed): {suffix}")
	max_idx = 0
	max_wt = float('-inf')
	for idx in range(len(prefix)):
		wt = prefix[idx] + suffix[-idx - 1]
		if wt > max_wt:
			max_idx = idx
			max_wt = wt
	print(f"Max index: {max_idx}, Max weight: {max_wt}")
	result.append((max_idx + i, mid_j))
	print(f"Intermediate result: {result}")
	hirschberg_helper(seq_1, seq_2, i, j, max_idx + i, mid_j, result)
	hirschberg_helper(seq_1, seq_2, max_idx + i, mid_j, i_p, j_p, result)

def hirschberg(seq_1, seq_2):
	result = []
	hirschberg_helper(seq_1, seq_2, 0, 0, len(seq_1), len(seq_2), result)
	result_no_dupes = set(result)
	print(f"Result before sorting: {result_no_dupes}")
	result_sorted = sorted(list(result_no_dupes), key = lambda x: (x[0], x[1]))
	print(f"Sorted result: {result_sorted}")
	print(result_sorted)
	return result_sorted

print("Final result:", hirschberg('ATGTTAT', 'ATCGTAC'))
