START = -1
INSERT = 0
DELETE = 1
SUBSTITUTE = 2


def global_align(v, w, match, mismatch, indel):
    """
    Finds an optimal global alignment of v and w using Needleman-Wunsch.
    :param v: first string to align
    :param w: other string to align
    :param match: score for matches
    :param mismatch: score for mismatches
    :param indel: score for insertions and deletions
    :return: a tuple (score, v_aligned, w_aligned) with the optimal score
    and aligned strings
    """
    m = len(v)
    n = len(w)

    # Initialization; D[i][j][0] contains the max alignment score of the
    # ith prefix of v and the jth of w; D[i][j][1] contains the back pointer.
    D = [[(0, START) for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = (i * indel, DELETE)

    for j in range(1, n + 1):
        D[0][j] = (j * indel, INSERT)

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = D[i][j-1][0] + indel
            delete = D[i-1][j][0] + indel
            substitute = D[i-1][j-1][0] + (match if v[i-1] == w[j-1]
                                           else mismatch)

            # Set D[i][j] to the max of the recurrences
            if 0 > insert and 0 > delete and 0 > substitute:
                D[i][j] = (0, START)
            elif insert > delete and insert > substitute:
                D[i][j] = (insert, INSERT)
            elif delete > substitute:
                D[i][j] = (delete, DELETE)
            else:
                D[i][j] = (substitute, SUBSTITUTE)

    for row in D:
        print(row)

    # Traceback starting at max cell and ending at first 0 encountered
    i, j = m, n
    v_aligned = ''
    w_aligned = ''
    back_pointer = D[i][j][1]
    while back_pointer != START:
        if back_pointer == INSERT:
            j -= 1
            v_aligned = '-' + v_aligned
            w_aligned = w[j] + w_aligned
        elif back_pointer == DELETE:
            i -= 1
            v_aligned = v[i] + v_aligned
            w_aligned = '-' + w_aligned
        elif back_pointer == SUBSTITUTE:
            i -= 1
            j -= 1
            v_aligned = v[i] + v_aligned
            w_aligned = w[j] + w_aligned
        back_pointer = D[i][j][1]

    return D[m][n][0], v_aligned, w_aligned



def centerStar(listofSeq):
    
    match = 0
    mistmatch = 1
    indel = 1
    
    seqLen = len(listofSeq)
    pwMatrix = [["-"]*seqLen for i in range(seqLen)]
    
    for seq in listofSeq:
        for seq2 in listofSeq:
            # in1 gives row, in2 gives column 
            in1 = listofSeq.index(seq)
            in2 = listofSeq.index(seq2)
            pwMatrix[in1][in2] = global_align(seq, seq2, match, mismatch, indel)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
