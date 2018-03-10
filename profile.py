START = -1
INSERT = 0
DELETE = 1
SUBSTITUTE = 2


with open('blosum62.txt') as f:
    blosum_file = f.readlines()

proteins = blosum_file[0].split()

blosum = {}
for line in blosum_file[1:]:
    line = line.split()
    for index, entry in enumerate(line[1:]):
        blosum[line[0], proteins[index]] = int(entry)


class Profile:
    def __init__(self, sequences, alignments):
        n = len(alignments)
        self.length = len(alignments[0])

        self.sequences = sequences
        self.alignments = alignments
        self.frequency = {(protein, i): 0 for protein in proteins
                          for i in range(len(self))}

        for alignment in alignments:
            for index, protein in enumerate(alignment):
                self.frequency[protein, index] += 1 / n

    def __len__(self):
        return self.length

    def __getitem__(self, item):
        return self.frequency[item]


def psp_empty(p, i):
    return sum(p[x, i] * blosum[x, '-'] for x in proteins)


def psp(p1, p2, i, j):
    return sum(p1[x, i]*p2[y, j]*blosum[x, y]
               for x in proteins for y in proteins)


def profile_align(p1, p2):
    """
    Finds an optimal global alignment of p1 and p2 using Needleman-Wunsch.
    :param p1: first profile to align
    :param p2: other profile to align
    :return: a tuple (score, v_aligned, w_aligned) with the optimal score
    and aligned strings
    """
    m = len(p1)
    n = len(p2)

    # Initialization; D[i][j][0] contains the max alignment score of the
    # ith prefix of p1 and the jth of p2; D[i][j][1] contains the back pointer.
    D = [[(0, START) for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = (D[i-1][0][0] + psp_empty(p1, i-1), DELETE)

    for j in range(1, n + 1):
        D[0][j] = (D[0][j-1][0] + psp_empty(p2, j-1), INSERT)

    # Recurrence
    for i in range(1, m + 1):
        print(i)
        for j in range(1, n + 1):
            insert = D[i][j - 1][0] + psp_empty(p2, j-1)
            delete = D[i - 1][j][0] + psp_empty(p1, i-1)
            substitute = D[i - 1][j - 1][0] + psp(p1, p2, i-1, j-1)

            # Set D[i][j] to the max of the recurrences
            if insert > delete and insert > substitute:
                D[i][j] = (insert, INSERT)
            elif delete > substitute:
                D[i][j] = (delete, DELETE)
            else:
                D[i][j] = (substitute, SUBSTITUTE)

    # Traceback
    p1_alignments = [''] * len(p1.alignments)
    p2_alignments = [''] * len(p2.alignments)
    i, j = m, n
    back_pointer = D[m][n][1]
    while back_pointer != START:
        if back_pointer == INSERT:
            j -= 1
            p1_alignments = ['-' + alignment for alignment in p1_alignments]
            p2_alignments = [p2.alignments[k][j] + p2_alignments[k]
                             for k in range(len(p2_alignments))]
        elif back_pointer == DELETE:
            i -= 1
            p1_alignments = [p1.alignments[k][i] + p1_alignments[k]
                             for k in range(len(p1_alignments))]
            p2_alignments = ['-' + alignment for alignment in p2_alignments]

        elif back_pointer == SUBSTITUTE:
            i -= 1
            j -= 1
            p1_alignments = [p1.alignments[k][i] + p1_alignments[k]
                             for k in range(len(p1_alignments))]

            p2_alignments = [p2.alignments[k][j] + p2_alignments[k]
                             for k in range(len(p2_alignments))]

        back_pointer = D[i][j][1]

    combined = Profile(p1.sequences + p2.sequences,
                       p1_alignments + p2_alignments)

    print('.', end='', flush=True)

    return combined
