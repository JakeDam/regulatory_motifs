# Calculates hamming distance between two nucleotide sequences 
def hamming_dist(seqA, seqB):
    h_dist = 0 
    for i in range(0, len(seqA)):
        if seqA[i] != seqB[i]:
            h_dist += 1
    return h_dist

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


seqs = ['CGGGTAACCTCCCTGCTAAGGGTCC', 'ATTCCACATGCTCAGTAAGTAATAA', 'AGTCAGACGGCTCAGAGCGTTTCAA', 'CTAAGCAGAGAAATGTCGGGTGGAC', 'TAATGGTTATCTAAGTCTGAGTGGG', 'TCTTTTTACGCTCAGGTCACTTTTC']

print(*motif_enumeration(seqs, 5, 1))


