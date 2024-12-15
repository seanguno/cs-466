import numpy as np

class PairHMM:
    def __init__(self, transition_probs, emission_probs):
        self.transition_probs = transition_probs
        self.emission_probs = emission_probs
        # M = match/mismatch ; I = insertion ; D = deletion
        self.states = ['M', 'I', 'D']

    def generate_sequences_and_score(self, length):
        seq1_aligned = []
        seq2_aligned = []
        seq1_raw = []
        seq2_raw = []
        current_state = 'M' # start in the match state
        score = 0 # lower bound for NW and HB scores 
        last_pair = None # keep track of last pair to encourage alternation

        while len(seq1_raw) < length or len(seq2_raw) < length:
            if current_state == 'M':
                pairs, probs = zip(*self.emission_probs['M'].items())
                if last_pair:
                    # case for ATAT generation
                    if ('A', 'A') in pairs:
                        # encourage alternation for matches
                        probs = [p * 3 if pair == ('T', 'T') and last_pair == ('A', 'A') \
                                 or pair == ('A', 'A') and last_pair == ('T', 'T') \
                                    else p for p, pair in zip(probs, pairs)]
                    # case for CpG island generation
                    elif ('C', 'C') in pairs:
                        # encourage alternation for matches
                        probs = [p * 3 if pair == ('C', 'C') and last_pair == ('G', 'G') \
                                 or pair == ('G', 'G') and last_pair == ('C', 'C') \
                                    else p for p, pair in zip(probs, pairs)]
                    # double probability if same letter for first sequence as last
                    probs = [p * 2 if pair != last_pair else p for p, pair in zip(probs, pairs)]
                    # renormalize probabilities
                    probs = [p / sum(probs) for p in probs]
                pair = np.random.choice(len(pairs), p=probs)
                last_pair = pairs[pair]
                seq1_aligned.append(pairs[pair][0])
                seq2_aligned.append(pairs[pair][1])
                seq1_raw.append(pairs[pair][0])
                seq2_raw.append(pairs[pair][1])
                # +1 if match, -1 else
                score += 1 if pairs[pair][0] == pairs[pair][1] else -1
            elif current_state == 'I':
                symbols, probs = zip(*self.emission_probs['I'].items())
                symbol = np.random.choice(len(symbols), p=probs)
                seq1_aligned.append(symbols[symbol])
                seq2_aligned.append('-')
                seq1_raw.append(symbols[symbol])
                # gap penalty
                score -= 1  
            elif current_state == 'D':
                symbols, probs = zip(*self.emission_probs['D'].items())
                symbol = np.random.choice(len(symbols), p=probs)
                seq1_aligned.append('-')
                seq2_aligned.append(symbols[symbol])
                seq2_raw.append(symbols[symbol])
                # gap penalty
                score -= 1  
 
            current_state = np.random.choice(self.states, p=self.transition_probs[current_state])

        return (''.join(seq1_aligned[:length]), ''.join(seq2_aligned[:length])), \
               (''.join(seq1_raw[:length]), ''.join(seq2_raw[:length])), score

# Parameters for ATAT repeats
transition_probs_atat = {
    'M': [0.9, 0.05, 0.05],
    'I': [0.5, 0.5, 0.0],
    'D': [0.5, 0.0, 0.5]
}

emission_probs_atat = {
    'M': {('A', 'T'): 0.1, ('T', 'A'): 0.1, ('A', 'A'): 0.3, ('T', 'T'): 0.3, ('C', 'G'): 0.1, ('G', 'C'): 0.1},
    'I': {'A': 0.4, 'T': 0.4, 'C': 0.1, 'G': 0.1},
    'D': {'A': 0.4, 'T': 0.4, 'C': 0.1, 'G': 0.1}
}

# Parameters for CpG islands
transition_probs_cpg = {
    'M': [0.8, 0.1, 0.1],
    'I': [0.5, 0.5, 0.0],
    'D': [0.5, 0.0, 0.5]
}

emission_probs_cpg = {
    'M': {('C', 'G'): 0.1, ('G', 'C'): 0.1, ('C', 'C'): 0.3, ('G', 'G'): 0.3, ('A', 'T'): 0.1, ('T', 'A'): 0.1},
    'I': {'A': 0.1, 'T': 0.1, 'C': 0.4, 'G': 0.4},
    'D': {'A': 0.1, 'T': 0.1, 'C': 0.4, 'G': 0.4}
}

# Function to save sequences and scores to files
def save_sequence_and_score(sequence, score, filename):
    with open(filename, "w") as file:
        file.write(f"Sequence: {sequence}\n")
        file.write(f"Score: {score}\n")

def save_sequence_fasta(sequence, score, identifier, filename):
    with open(filename, "w") as file:
        file.write(f">{identifier}_score: {score}\n")
        # Write sequence in lines of 70 characters
        for i in range(0, len(sequence), 70):
            file.write(sequence[i:i+70] + "\n")

# Generate ATAT repeat sequences
pair_hmm_atat = PairHMM(transition_probs_atat, emission_probs_atat)
(seq1_atat_aligned, seq2_atat_aligned), (seq1_atat_raw, seq2_atat_raw), atat_score = pair_hmm_atat.generate_sequences_and_score(1000)

# Generate CpG island sequences
pair_hmm_cpg = PairHMM(transition_probs_cpg, emission_probs_cpg)
(seq1_cpg_aligned, seq2_cpg_aligned), (seq1_cpg_raw, seq2_cpg_raw), cpg_score = pair_hmm_cpg.generate_sequences_and_score(1000)


# Save ATAT sequences and scores
save_sequence_fasta(seq1_atat_raw, atat_score, "seq1_ATAT_raw", "seq1_ATAT_raw.fasta")
save_sequence_fasta(seq2_atat_raw, atat_score, "seq2_ATAT_raw", "seq2_ATAT_raw.fasta")
save_sequence_fasta(seq1_atat_aligned, atat_score, "seq1_ATAT_aligned", "seq1_ATAT_aligned.fasta")
save_sequence_fasta(seq2_atat_aligned, atat_score, "seq2_ATAT_aligned", "seq2_ATAT_aligned.fasta")

# Save CpG sequences and scores
save_sequence_fasta(seq1_cpg_aligned, cpg_score, "seq1_CpG_aligned", "seq1_CpG_aligned.fasta")
save_sequence_fasta(seq2_cpg_aligned, cpg_score, "seq2_CpG_aligned", "seq2_CpG_aligned.fasta")
save_sequence_fasta(seq1_cpg_raw, cpg_score, "seq1_CpG_raw", "seq1_CpG_raw.fasta")
save_sequence_fasta(seq2_cpg_raw, cpg_score, "seq2_CpG_raw", "seq2_CpG_raw.fasta")
