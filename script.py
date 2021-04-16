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
