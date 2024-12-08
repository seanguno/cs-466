def needleman_score_linear(seq_1, seq_2, debug=False):
	# init dp table
	score = [[0] * 2 for i in range(len(seq_1) + 1)]
	# init first col
	for i in range(1, len(seq_1) + 1):
		# in future could include diff penalty dictionaries 
		score[i][0] = score[i - 1][0] - 1
	if debug:
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
		if debug:
			print(f'col {j}: {[row[1] for row in score]}')
	# get last row of score
	last_col = [row[1] for row in score]
	if debug:
		print(f"Last col for '{seq_1}' vs '{seq_2}': {last_col}")
	return last_col

def hirschberg_helper(seq_1, seq_2, i, j, i_p, j_p, result, debug=False):
	if j_p - j == 1:
		result.append((i, j))
		result.append((i_p, j_p))
		if debug:
			print(f"Base case 1 result: {result}")
		return
	if j_p < 1:
		if debug:
			print(f"Base case 2 (j_p < 1): No operation.")
		return
	mid_j = int((j + j_p) / 2)
	prefix = needleman_score_linear(seq_1[i:i_p], seq_2[j:mid_j])
	suffix = needleman_score_linear(seq_1[i:i_p][::-1], seq_2[mid_j:j_p][::-1])
	if debug:
		print(f"Prefix scores: {prefix}")
	if debug:	
		print(f"Suffix scores (reversed): {suffix}")
	max_idx = 0
	max_wt = float('-inf')
	for idx in range(len(prefix)):
		wt = prefix[idx] + suffix[-idx - 1]
		if wt > max_wt:
			max_idx = idx
			max_wt = wt
	if debug:
		print(f"Max index: {max_idx}, Max weight: {max_wt}")
	result.append((max_idx + i, mid_j))
	if debug:
		print(f"Intermediate result: {result}")
	hirschberg_helper(seq_1, seq_2, i, j, max_idx + i, mid_j, result)
	hirschberg_helper(seq_1, seq_2, max_idx + i, mid_j, i_p, j_p, result)

def hirschberg(seq_1, seq_2, debug=False):
	result = []
	hirschberg_helper(seq_1, seq_2, 0, 0, len(seq_1), len(seq_2), result)
	result_no_dupes = set(result)
	if debug:
		print(f"Result before sorting: {result_no_dupes}")
	result_sorted = sorted(list(result_no_dupes), key = lambda x: (x[0], x[1]))
	if debug:
		print(f"Sorted result: {result_sorted}")
		print(result_sorted)
	return result_sorted

def fill_in_backtrace(seq_1, seq_2, debug=False):
	result = hirschberg(seq_1, seq_2)
	filled_in = []
	for cell_idx in range(len(result) - 1):
		filled_in.append(result[cell_idx])
		match_found = False
		if result[cell_idx + 1][0] - result[cell_idx][0] > 1:
			if debug:
				print(f'second {result[cell_idx + 1][0]} - first {result[cell_idx][0]} is {result[cell_idx + 1][0] - result[cell_idx][0]}')
			# need to account for gap
			for i in range(result[cell_idx][0], result[cell_idx + 1][0] - 1):
				if seq_1[i] == seq_2[result[cell_idx][1]]:
					if debug:
						print(f'i: {i}')
						print(f'result[cell_idx][1]: {result[cell_idx][1]}')
						print(f'curr cell: {result[cell_idx]}')
					filled_in.append((i + 1, result[cell_idx][1] + 1))
					match_found = True
				else:
					if debug:
						print(f'match found is: {match_found}')
						print(f'appending: {(i + 1, result[cell_idx][1] + 1 if match_found else result[cell_idx][1])}')
					filled_in.append((i + 1, result[cell_idx][1] + 1 if match_found else result[cell_idx][1]))
	filled_in.append(result[-1])
	return filled_in

def print_alignment(seq_1, seq_2):
	alignment_path = fill_in_backtrace(seq_1, seq_2)
	align_1 = []
	align_2 = []
	i, j = 0, 0

	for coord in alignment_path[1:]:  # start from second coord
		while i < coord[0] and j < coord[1]:
			align_1.append(seq_1[i])
			align_2.append(seq_2[j])
			i += 1
			j += 1
		while i < coord[0]:
			align_1.append(seq_1[i])
			align_2.append('-')
			i += 1
		while j < coord[1]:
			align_1.append('-')
			align_2.append(seq_2[j])
			j += 1

	# add any left over characters
	align_1.extend(seq_1[i:])
	align_2.extend(seq_2[j:])

	# pad the shorter sequence with gaps if necessary
	max_len = max(len(align_1), len(align_2))
	align_1 += ['-'] * (max_len - len(align_1))
	align_2 += ['-'] * (max_len - len(align_2))

	# create the match line
	match_line = ['|' if a == b and a != '-' else ' ' for a, b in zip(align_1, align_2)]

	# print the alignment
	print('Sequence 1:', ''.join(align_1))
	print('           ', ''.join(match_line))
	print('Sequence 2:', ''.join(align_2))

	# Calculate and print alignment statistics
	matches = sum(1 for a, b in zip(align_1, align_2) if a == b and a != '-')
	mismatches = sum(1 for a, b in zip(align_1, align_2) if a != b and a != '-' and b != '-')
	gaps = align_1.count('-') + align_2.count('-')
	print(f'Score: {matches - gaps - mismatches}')



# print("filled in", fill_in_backtrace("ATGTTAT", "ATCGTAC"))
			
# print("filled in", fill_in_backtrace("ATCGTACTTTTTT", "ATGTTAT"))
# print("filled in", fill_in_backtrace("ATGTTAT", "ATCGTACTTTTTT"))
print_alignment("ATGTAGTACAAAAAA", "ATCGTAC")



					
		

#print("Final result:", hirschberg('ATGTTAT', 'ATCGTAC'))
