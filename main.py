from hirschberg import print_alignment
from hirschberg import fill_in_backtrace
from hirschberg import hirschberg
from needleman_wunsch import global_align

import timeit

# put sequences from files into an array
num_seqs = 16   # number of files (genes)
seqs = ['' for i in range(num_seqs+1)]

# reading files
# files are named by numbers starting at 1, so index 0 of seqs remains empty
for i in range(1, num_seqs+1):
    file_name = str(i) + '.fasta'
    with open(file_name, 'r') as file:
        for line in file:
            if line[0] != '>':
                seqs[i] += line.strip()

# dictionary where key is tuple of sequences aligned (v, w)
# and value is tuple of scores (needleman score, hirschberg score)
# scores = {}
# alignment_keys = [(1, 2)]   # alignments to perform
# #(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14), (15, 16)
# for i, j in alignment_keys:
#     v = seqs[i]
#     w = seqs[j]
#     scores[(i, j)] = (global_align(v, w)[0], print_alignment(v, w, fill_in_backtrace(v, w)))

# align just sequences v & w
# print(scores)
v = seqs[15][:4000]
w = seqs[16][:4000]
#NW = global_align(v, w)
#HB = print_alignment(v, w, fill_in_backtrace(v, w))
#print("NW: " + str(NW[0]))
#print("HB: " + str(HB))


# time analysis
# array of tuples (sequence length, alignment time) 

# nw_time = []
# hb_time = []
# for i, j in alignment_keys:
#     v = seqs[i]
#     w = seqs[j]
#     # measure time to perform each alignment once, repeat 5 times to get an average
#     nw = timeit.timeit(lambda: global_align(v, w), number=1)
#     hb = timeit.timeit(lambda: print_alignment(v, w, fill_in_backtrace(v, w)), number=1)

#     len_avg = (len(v) + len(w)) / 2
#     nw_time.append((len_avg, nw))
#     hb_time.append((len_avg, hb))

# for i in range(1, 16):
#     length = i * 1000
#     v = seqs[15][:length]
#     w = seqs[16][:length]
#     nw = timeit.timeit(lambda: global_align(v, w), number=1)
#     hb = timeit.timeit(lambda: print_alignment(v, w, fill_in_backtrace(v, w)), number=1)

#     nw_time.append((length, nw))
#     hb_time.append((length, hb))
#     print(f"Length {length} done")

# print(f"NW time: {nw_time} seconds")
# print(f"HB time: {hb_time} seconds")

