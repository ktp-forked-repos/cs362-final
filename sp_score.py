from profile import blosum

with open('progressive_alignment_comp.txt') as f:
    alignments = [line.strip() for line in f.readlines()]

sp_score = 0
for x in alignments:
    for y in alignments:
        for xi, yi in zip(x, y):
            sp_score += blosum[xi, yi]

print('Sum-of-Pair Score:', sp_score)