"""
Script to compute the Sum-of-Pair score of a multiple sequence alignment.
"""


from profile import blosum
import sys

if len(sys.argv) != 2:
    print('Usage: python3 sp_score.py msa_file')
    exit()

with open(sys.argv[1]) as f:
    alignments = [line.strip() for line in f.readlines()]

sp_score = 0
for x in alignments:
    for y in alignments:
        for xi, yi in zip(x, y):
            sp_score += blosum[xi, yi]

print('Sum-of-Pair Score:', sp_score)