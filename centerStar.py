from parse import parse_fasta

START = -1
INSERT = 0
DELETE = 1
SUBSTITUTE = 2


# Load in BLOSUM62 matrix
with open('blosum62.txt') as f:
    blosum_file = f.readlines()

proteins = blosum_file[0].split()

blosum = {}
for line in blosum_file[1:]:
    line = line.split()
    for index, entry in enumerate(line[1:]):
        blosum[line[0], proteins[index]] = int(entry)

    
def pairwise(string_v, string_w):
    """
    Finds an pairwise distance of v and w using Needleman-Wunsch and taking gaps
    into consideration.
    :param string_v: first string to align
    :param string_w: other string to align
    :return: a tuple (score, v_aligned, w_aligned) with the optimal score
    and aligned strings
    """
    m = len(string_v)
    n = len(string_w)

    # Initialization; D[i][j][0] contains the max alignment score of the
    # ith prefix of v and the jth of w; D[i][j][1] contains the back pointer.
    D = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = i

    for j in range(1, n + 1):
        D[0][j] = j

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = D[i][j-1] + 1
            delete = D[i-1][j] + 1
            substitute = D[i-1][j-1] + (0 if (string_v[i-1] == string_w[j-1]) else 1)
            # Set D[i][j] to the max of the recurrences
            if insert < delete and insert < substitute:
                D[i][j] = insert
            elif delete < substitute:
                D[i][j] = delete
            else:
                D[i][j] = substitute

    return D[m][n]

def findCenterSeq(dictofSeq):
    """
    Finds the center sequence by taking the sequence with minimum of sum of edit
    distances.
    :param dictofSeq: Sequences passed by parse.py
    :return: the Name of center sequence
    """
    seqLen = len(dictofSeq)
    pwMatrix = [["-"]*seqLen for i in range(seqLen)]
    listofSeq = []
    for key in dictofSeq:
        listofSeq.append(dictofSeq.get(key))
    
    findMin = []
    acc = 0
    for seq in listofSeq:
        for seq2 in listofSeq:
            # in1 gives row, in2 gives column 
            in1 = listofSeq.index(seq)
            in2 = listofSeq.index(seq2)
            pwMatrix[in1][in2] = pairwise(seq, seq2)
            acc += pwMatrix[in1][in2]
            #TypeError: 'int' object is not subscriptable
        findMin.append(acc)
        acc = 0
    posSeq = findMin.index(min(findMin))
    refString = listofSeq[posSeq]
    refName = ""
    
    for name, seq in dictofSeq.items():
        if seq == refString:
            refName = name
    
    print(refName)
    
    return refName

def sequence_align(string_v, string_w):
    """
    Finds an optimal global alignment of string v and string w using Needleman-Wunsch.
    :param string_v: first string to align
    :param string_w: other string to align
    :return: a tuple (v_aligned, w_aligned) of aligned strings
    """
    m = len(string_v)
    n = len(string_w)

    # Initialization; D[i][j][0] contains the max alignment score of the
    # ith prefix of v and the jth of w; D[i][j][1] contains the back pointer.
    D = [[(0, START) for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = (D[i - 1][0][0] + blosum['-', string_v[i - 1]], DELETE)

    for j in range(1, n + 1):
        D[0][j] = (D[0][j - 1][0] + blosum['-', string_w[j - 1]], INSERT)

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = D[i][j-1][0] + blosum['-', string_w[j - 1]]
            delete = D[i-1][j][0] + blosum[string_v[i - 1], '-']
            substitute = D[i-1][j-1][0] + blosum[string_v[i - 1], string_w[j - 1]]
            # Set D[i][j] to the max of the recurrences
            if insert > delete and insert > substitute:
                D[i][j] = (insert, INSERT)
            elif delete > substitute:
                D[i][j] = (delete, DELETE)
            else:
                D[i][j] = (substitute, SUBSTITUTE)

    i, j = m, n
    v_aligned = ''
    w_aligned = ''
    back_pointer = D[i][j][1]
    while back_pointer != START:
        if back_pointer == INSERT:
            j -= 1
            v_aligned = '-' + v_aligned
            w_aligned = string_w[j] + w_aligned

            
        elif back_pointer == DELETE:
            i -= 1
            v_aligned = string_v[i] + v_aligned
            w_aligned = '-' + w_aligned

        elif back_pointer == SUBSTITUTE:
            i -= 1
            j -= 1
            v_aligned = string_v[i] + v_aligned
            w_aligned = string_w[j] + w_aligned

                    
        back_pointer = D[i][j][1]
        
    return v_aligned, w_aligned

def gap_align(center, string_w):
    """
    Updates global alignment of string v and string w using Needleman-Wunsch.
    Prevents insertion during this alignment because string_v is always the center string.
    :param string_v: first string to align
    :param string_w: other string to align
    :return: a tuple (v_aligned, w_aligned) of aligned strings
    """
    m = len(center)
    n = len(string_w)

    # Initialization; D[i][j][0] contains the max alignment score of the
    # ith prefix of v and the jth of w; D[i][j][1] contains the back pointer.
    D = [[(0, START) for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        D[i][0] = (D[i - 1][0][0] + blosum['-', center[i - 1]], DELETE)

    for j in range(1, n + 1):
        D[0][j] = (D[0][j - 1][0] + blosum['-', string_w[j - 1]], INSERT)

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            delete = D[i-1][j][0] + blosum[center[i - 1], '-']
            substitute = D[i-1][j-1][0] + blosum[center[i - 1], string_w[j - 1]]
            # Set D[i][j] to the max of the recurrences
            if delete > substitute:
                D[i][j] = (delete, DELETE)
            else:
                D[i][j] = (substitute, SUBSTITUTE)

    i, j = m, n
    w_aligned = ''
    back_pointer = D[i][j][1]
    while back_pointer != START:
        if back_pointer == DELETE:
            i -= 1
            w_aligned = '-' + w_aligned

        elif back_pointer == SUBSTITUTE:
            i -= 1
            j -= 1
            w_aligned = string_w[j] + w_aligned

                    
        back_pointer = D[i][j][1]
        
    return w_aligned

def centerStar_align(refName, dictofSeq):
    """
    Aligns all the sequences with Center Star MSA using Needleman-Wunsch.
    :param refString: the center sequence
    :param listofSeq: all the sequences need to be aligned
    :return: a list of aligned sequences
    """
    dictofFinalStr = {}
    refString = dictofSeq.pop(refName)
    #remove the center sequence from the list of sequence so it won't align to itself
    centerString = refString
    #construct a pointer to center squence
    for name in dictofSeq:
        alignment = sequence_align(centerString, dictofSeq.get(name))
        centerString = alignment[0]
        #print(centerString)
        strAligned = alignment[1]
        #print(strAligned)
        dictofFinalStr[name] = strAligned
        #print(len(listofFinalStr))

    for seq in dictofFinalStr:
        #Aligns all the sequence to the final center sequence with all the gaps inserted
        finalScore = gap_align(centerString, dictofFinalStr[seq])
        finalStr = finalScore
        dictofFinalStr[seq] = finalStr

    dictofFinalStr[refName] = (centerString)
    return dictofFinalStr

def main():
    print('Reading fasta file...')
    sequences = parse_fasta('COMP.txt')
    
    print('Performing center star alignment...', end='', flush=True)
    center = findCenterSeq(sequences)
    alignment = centerStar_align(center, sequences)

    with open('center_start_alignment.txt', 'w') as f:
        f.write('\n'.join(['{}: {}'.format(key, value)
                           for key, value in alignment.items()]))

    for name in alignment:
        print(alignment[name])

if __name__ == '__main__':
    main()
