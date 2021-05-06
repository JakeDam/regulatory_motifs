import math

# Calculates hamming distance between two nucleotide sequences 
def hamming_dist(seqA, seqB):
    h_dist = 0 
    for i in range(0, len(seqA)):
        if seqA[i] != seqB[i]:
            h_dist += 1
    return h_dist


# Generates the d-neighborhood of pattern
def neighbors(pattern, d):
    if d == 0:
        return set([pattern])
    if len(pattern) == 1:
        return set(('A', 'T', 'C', 'G'))
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if hamming_dist(pattern[1:], text) < d:
            for base in ['A', 'T', 'C', 'G']:
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood

# Identifies all (k,d)-motifs that appear in a collection of DNA sequences
def motif_enumeration(dna, k, d):
    patterns = set()
    for i in range(0, len(dna[0]) - k + 1):
        neighborhood = neighbors(dna[0][i:i + k], d)
        for neighbor in neighborhood:
            count = 0
            for sequence in dna:
                for j in range(0, len(sequence) - k + 1):
                    if hamming_dist(neighbor, sequence[j: j + k]) <= d:
                        count += 1
                        break
            if count == len(dna):
                patterns.add(neighbor)
    return patterns


# Calculates the entropy of a motif matrix
def entropy(matrix):
    column_entropies = []
    for i in range(0, len(matrix[0])):
        column = {}
        entropy = 0
        for sequence in matrix:
            if sequence[i] not in column:
                column[sequence[i]] = 1
            else:
                column[sequence[i]] += 1
        for value in column.values():
            if value > 0:
                value = value / len(matrix)
                value = value * math.log2(value)
                entropy += value
        entropy *= -1
        column_entropies.append(entropy)
    return sum(column_entropies)
    
# Calcuates the hamming distance of a pattern to a collection of DNA sequences
def pattern_sequence_dist(pattern, dna):
    k = len(pattern)
    distance = 0
    for sequence in dna:
        hamming_distance = float('inf')
        for i in range(len(sequence) - k + 1):
            if hamming_distance > hamming_dist(pattern, sequence[i:i+k]):
                hamming_distance = hamming_dist(pattern, sequence[i:i+k])
        distance += hamming_distance
    return distance

# Returns a k-mer pattern that minimizes the hamming distance between a Pattern and a collection of DNA sequences
def median_string(k, dna):
    distance = float('inf')
    median = None
    neighborhood = neighbors('A' * k, k)
    for i in neighborhood:
        if distance > pattern_sequence_dist(i, dna):
            distance = pattern_sequence_dist(i, dna)
            median = i
    return median

# Computes the probablity of a given k-mer from a profile matrix of k-mers
def probability(text, matrix):
    prob = 1
    for i in range(0, len(text)):
        if text[i] == 'A':
            prob *= matrix['A'][i]
        elif text[i] == 'C':
            prob *= matrix['C'][i]
        elif text[i] == 'G':
            prob *= matrix['G'][i]
        else:
            prob *= matrix['T'][i]
    return prob

# Determines the profile most probable k-mer from a given DNA sequence and profile matrix 
def prof_most_prob(text, k, prof):
    most_prob = 0
    most_prob_kmer = ''
    for i in range(0, len(text) - k + 1):
        k_mer = text[i:i + k]
        prob = probability(k_mer, prof)
        if prob > most_prob:
            most_prob = prob
            most_prob_kmer = k_mer
    if most_prob_kmer == '':
        return text[0:k]
    else:
        return most_prob_kmer

# Generates a profile matrix from a group of motifs
def generate_matrix(motifs):
    profile = {}
    A, C, G, T = [], [], [], []
    for j in range(len(motifs[0])):
        count_A, count_C, count_G, count_T = 1, 1, 1, 1
        for motif in motifs:
            if motif[j] == "A":
                count_A += 1
            elif motif[j] == "C":
                count_C += 1
            elif motif[j] == "G":
                count_G += 1
            elif motif[j] == "T":
                count_T += 1
        A.append(count_A)
        C.append(count_C)
        G.append(count_G)
        T.append(count_T)
    profile["A"] = A
    profile["C"] = C
    profile["G"] = G
    profile["T"] = T
    return profile

# Determines the consensus string from a group of motifs 
def consensus_string(motifs):
    consensus = ""
    for i in range(len(motifs[0])):
        count_A, count_C, count_G, count_T = 0, 0, 0, 0
        for motif in motifs:
            if motif[i] == "A":
                count_A += 1
            elif motif[i] == "C":
                count_C += 1
            elif motif[i] == "G":
                count_G += 1
            elif motif[i] == "T":
                count_T += 1
        if count_A >= max(count_C, count_G, count_T):
            consensus += "A"
        elif count_C >= max(count_A, count_G, count_T):
            consensus += "C"
        elif count_G >= max(count_C, count_A, count_T):
            consensus += "G"
        elif count_T >= max(count_C, count_G, count_A):
            consensus += "T"
    return consensus

# Scores a motif based on closeness to consensus string of motif matrix
def score(motifs):
  consensus = consensus_string(motifs)
  score = 0
  for motif in motifs:
    score += hamming_dist(consensus, motif)
  return score

# Returns a collection of k-mers representing the best motifs (most representing of each other) from a collection of DNA sequences
def greedy_motif_search(dna, k, t):
    best_motifs = []
    for sequence in dna:
        best_motifs.append(sequence[:k])
    base_motif = dna[0]
    other_motifs = dna[1:]
    for i in range(0, len(dna[0]) - k + 1):
        motifs = []
        motifs.append(dna[0][i:i + k])
        for motif in other_motifs:
            profile_matrix = generate_matrix(motifs)
            next_motif = prof_most_prob(motif, k, profile_matrix)
            motifs.append(next_motif)
        if score(motifs) < score(best_motifs):
            best_score = score(motifs)
            best_motifs = motifs
    return best_motifs

      


        



