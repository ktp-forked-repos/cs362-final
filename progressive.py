from parse import parse_fasta
from profile import Profile, profile_align, blosum
from neighbor_join import construct_tree

START = -1
INSERT = 0
DELETE = 1
SUBSTITUTE = 2


def align_score(v, w):
    """
    Finds the optimal alignment score of v and w.
    :param v: first string to align
    :param w: other string to align
    :return: the score of the alignment
    """
    m = len(v)
    n = len(w)

    # Initialization; D[i][j] contains the max alignment score of the
    # ith prefix of v and the jth of w
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = D[i-1][0] + blosum[v[i - 1], '-']

    for j in range(1, n + 1):
        D[0][j] = D[0][j-1] + blosum['-', w[j - 1]]

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = D[i][j - 1] + blosum['-', w[j - 1]]
            delete = D[i - 1][j] + blosum[v[i - 1], '-']
            substitute = D[i - 1][j - 1] + blosum[v[i - 1], w[j - 1]]

            D[i][j] = max(insert, delete, substitute)

    return D[m][n]


def multiple_align(node, sequences):
    """
    Recursively perform multiple sequence alignment along the guide tree
    :param node: the root node of the guide tree
    :param sequences: the dictionary of seqeucnes to align
    :return:
    """
    if len(node.children) == 0:
        return Profile([sequences[node.label]], [sequences[node.label]])

    return profile_align(multiple_align(node.children[0], sequences),
                         multiple_align(node.children[1], sequences))


def main():
    print('Reading fasta file...')
    sequences = parse_fasta('COMP.txt')

    print('Computing pairwise edit distances...', end='', flush=True)
    D = {(a, b): -align_score(sequences[a], sequences[b])
         for a in sequences for b in sequences}
    print()

    for a in sequences:
        print(a, end=' ')

    for a in sequences:
        print()
        print(a, end=' ')
        for b in sequences:
            print(D[a, b], end=' ')

    print()

    print('Constructing guide tree...')
    root = construct_tree(D, sequences)

    print('Performing progressive alignment...', end='', flush=True)
    profile = multiple_align(root, sequences)
    with open('progressive_alignment.txt', 'w') as f:
        f.write('\n'.join(profile.alignments))
    print()


if __name__ == '__main__':
    main()