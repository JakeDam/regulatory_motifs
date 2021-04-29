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
    return most_prob_kmer

# Returns a collection of k-mers representing the best motifs (most representing of each other) from a collection of DNA sequences
def greedy_motif_search(dna, k, t):
    best_motifs = []
    for sequence in dna:
        best_motifs.append(sequence[0:k])
    for i in range(0, len(dna[0]) - k + 1):
        motif = []
        motif.append(dna[0][i:i + k])



        
        

        


        



