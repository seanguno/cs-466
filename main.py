from hirschberg import print_alignment
from hirschberg import fill_in_backtrace
from hirschberg import hirschberg
from needleman_wunsch import global_align

import timeit
import matplotlib.pyplot as plt
import numpy as np

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
scores = {}
alignment_keys = [(15, 16)]   # alignments to perform
#(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14), (15, 16)
for i, j in alignment_keys:
    v = seqs[i]
    w = seqs[j]
    scores[(i, j)] = (global_align(v, w)[0], print_alignment(v, w, fill_in_backtrace(v, w)))

# print(scores)
#v = seqs[4]
#w = seqs[5]
#NW = global_align(v, w)
#HB = print_alignment(v, w, fill_in_backtrace(v, w))
#scores.append((global_align(seqs[8], seqs[9])[0], print_alignment(seqs[8], seqs[9], fill_in_backtrace(seqs[8], seqs[9]))))
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

# plot time analysis
# x-axis = average length of the two aligned sequences
# y-axis = average time of alignment (alignment performed 5 times, average time shown)
nw_time = [(19.0, 0.0004990999586880207), (479.0, 0.27225330006331205), (1001.5, 1.2385339997708797), (1763.5, 3.7759147002361715), (2229.0, 4.180357899982482), (6403.5, 34.75250760000199), (9307.0, 74.01478530000895), (16569.0, 236.11566699994728), (1000, 0.8331424002535641), (2000, 3.3964247000403702), (3000, 7.649441899731755), (4000, 13.700268200132996), (5000, 21.75459419982508), (6000, 32.62650189967826), (7000, 41.84669720008969), (8000, 55.991296799853444), (9000, 72.49378219991922), (10000, 86.82874789973721), (11000, 102.13932760013267), (12000, 122.52714310027659), (13000, 143.80858709989116), (14000, 219.46323610004038), (15000, 275.8117311000824)]
hb_time = [(19.0, 0.0006675003096461296), (479.0, 0.21720070019364357), (1001.5, 0.9112697001546621), (1763.5, 2.3741548000834882), (2229.0, 3.608909399714321), (6403.5, 36.08737460011616), (9307.0, 66.23498589964584), (16569.0, 223.51840800000355), (1000, 0.7073352998122573), (2000, 2.8966955998912454), (3000, 6.6584540996700525), (4000, 12.199166500009596), (5000, 19.02023210003972), (6000, 27.1168828997761), (7000, 37.55801279982552), (8000, 50.481557599734515), (9000, 65.34488129988313), (10000, 80.00038050021976), (11000, 95.23744589975104), (12000, 113.36187560018152), (13000, 135.7678286000155), (14000, 187.17405309993774), (15000, 186.1493144002743)]
x, y = zip(*nw_time)
x2, y2 = zip(*hb_time)

# fit quadratic polynomial to both
coeffs_nw = np.polyfit(x, y, 2)
poly_nw = np.poly1d(coeffs_nw)
print(coeffs_nw)

coeffs_hb = np.polyfit(x2, y2, 2)
poly_hb = np.poly1d(coeffs_hb)
print(coeffs_hb)
# x-vals for polynomial plotting
x_fit = np.linspace(min(min(x), min(x2)), max(max(x), max(x2)), 100)

# plotting
plt.scatter(x, y, color='blue', label='Needleman-Wunsch')
plt.scatter(x2, y2, color='red', label='Hirschberg')

plt.plot(x_fit, poly_nw(x_fit), color='blue', linestyle='--', label = 'Quadratic Fit for Needleman-Wunsch')
plt.plot(x_fit, poly_hb(x_fit), color='red', linestyle='--', label = 'Quadratic Fit for Hirschberg')

plt.xlabel('Average Sequence Length')
plt.ylabel('Averag Alignment Time (sec)')
plt.legend()
plt.show()
