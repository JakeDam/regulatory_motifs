# Regulatory Motifs

## Coding Challenges from UCSD's Bioinformatics I Week 3

Python script that contains tools for characterizing regulatory motifs:
- hamming_dist: Calculates the Hamming Distance between two DNA sequences
- neighbors: Generates the d-neighborhood (set of all k-mers whose hamming distance from the pattern does not exceed d) of a DNA sequence 
- motif_enumeration: Identifies all (k,d)-motifs that appear in a collection of DNA sequences
- entropy: Calculates the entropy of a motif matrix
- pattern_sequence_dist: Calcuates the hamming distance of a pattern to a collection of DNA sequences
- median_string: Returns a k-mer pattern that minimizes the hamming distance between a Pattern and a collection of DNA sequences
- probability: Computes the probablity of a given k-mer from a profile matrix of k-mers
- prof_most_prob: Determines the profile most probable k-mer from a given DNA sequence and profile matrix 
- generate_matrix: Generates a profile matrix from a group of motifs
- consensus_string: Determines the consensus string from a group of motifs 
- score: Scores a motif based on closeness to consensus string of motif matrix
- greedy_motif_search: Returns a collection of k-mers representing the best motifs (most representing of each other) from a collection of DNA sequences