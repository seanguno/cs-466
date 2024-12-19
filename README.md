# cs-466: Needleman-Wunsch and Hirshberg Comparison 

Our project compares the Needleman-Wunsch and Hirshberg algorithms to determine the space-time complexities. Through exploring the two algorithm’s trade-offs, we will determine if Hirschberg truly optimizes space complexity to be linear (compared to quadratic) when aligning nucleotide sequences.  
<br>

# Installation
## Prerequisites
- **Python** - Our program requires Python 3.9.12 or newer. 
## Install using PIP
Our program requires some python libraries, which can be installed using `pip`, the package manager for Python. Open a terminal or command prompt and run the following command:

```bash
$ pip install matplotlib memory_profiler
```

# Run Code
To run the code, open a terminal or command prompt and run the following command:
```bash
$ python main.py
```

# Usage Instructions
Sequences to align should be formatted as FASTA files and named using integers. For example, 1.fasta and 2.fasta. In main.py, the `num_seqs` variable should be set to the number of sequences (number of FASTA files) that will be analyzed. `alignment_keys` should be filled in with the alingments to compute, for example, `alignment_keys = [(1, 2)]` indicates 1.fasta and 2.fasta will be aligned. Note that the code will take a very long time to run if long sequences are aligned, expecially if `@profile` is not commented out beacuse memory_profiler increases the overall runtime.


To plot time and memory analysis results, run `plot.py` after filling in `nw_time` and `hb_time` for Needleman-Wunsch and Hirschberg runtimes, respectively. The same should be done for `nw_mem` and `hb_mem` to plot memory useage.
