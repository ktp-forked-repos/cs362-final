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

    
def pairwise(v, w):
    """
    Finds an pairwise distance of v and w using Needleman-Wunsch and taking gaps
    into consideration.
    :param v: first string to align
    :param w: other string to align
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
            insert = D[i][j-1][0] + 1
            delete = D[i-1][j][0] + 1
            substitute = D[i-1][j-1][0] + (0 if ((v[i-1] == w[j-1]) or (v[i-i] == "-") or (w[j-1] == "-")) else 1)
            # Set D[i][j] to the max of the recurrences
            if insert < delete and insert < substitute:
                D[i][j] = (insert, INSERT)
            elif delete < substitute:
                D[i][j] = (delete, DELETE)
            else:
                D[i][j] = (substitute, SUBSTITUTE)

    return D[m][n][0]

def findCenterSeq(listofSeq):
    """
    Finds the center sequence by taking the sequence with minimum of sum of edit
    distances.
    :param listofSeq: Sequences passed by parse.py
    :return: the center sequence
    """
    seqLen = len(listofSeq)
    pwMatrix = [["-"]*seqLen for i in range(seqLen)]
    
    findMin = []
    acc = 0
    for seq in listofSeq:
        for seq2 in listofSeq:
            # in1 gives row, in2 gives column 
            in1 = listofSeq.index(seq)
            in2 = listofSeq.index(seq2)
            pwMatrix[in1][in2] = pairwise(seq, seq2)
            acc += pwMatrix[in1][in2][0]
        findMin.append(acc)
        acc = 0
    posSeq = findMin.index(min(findMin))
    
    return listofSeq[posSeq]

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
        D[i][0] = (i * indel, DELETE)

    for j in range(1, n + 1):
        D[0][j] = (j * indel, INSERT)

    # Recurrence
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            insert = D[i][j-1][0] + blosum['-', w[j - 1]]
            delete = D[i-1][j][0] + blosum[v[i - 1], '-']
            substitute = D[i-1][j-1][0] + blosum[v[i - 1], w[j - 1]]
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
        
    return v_aligned, w_aligned

def centerStar_align(refString, listofSeq):
    """
    Aligns all the sequences with Center Star MSA using Needleman-Wunsch.
    :param refString: the center sequence
    :param listofSeq: all the sequences need to be aligned
    :return: a list of aligned sequences
    """
    match = 0
    listofFinalStr = []
    listofSeq.remove(refString)
    #remove the center sequence from the list of sequence so it won't align to itself
    centerString = refString
    #construct a pointer to center squence
    l = len(listofSeq)
    for i in range(l):
        alignment = sequence_align(centerString, listofSeq[i])
        centerString = alignment[1]
        #print(centerString)
        strAligned = alignment[2]
        #print(strAligned)
        listofFinalStr.append(strAligned)
        #print(len(listofFinalStr))
    
    for j in range(len(listofFinalStr)):
        #Aligns all the sequence to the final center sequence with all the gaps inserted
        finalScore = sequence_align(centerString, listofFinalStr[j])
        finalStr = finalScore[2]
        listofFinalStr[j] = finalStr

    listofFinalStr.append(centerString)
    return listofFinalStr

def main():
    print('Reading fasta file...')
    sequences = parse_fasta('COMP.txt')
    
    print('Performing center star alignment...', end='', flush=True)
    center = findCenterSeq(sequences)
    #For parsing, it does not pass in a list, still need to work on change dictionary into a list
    alignment = centerStar_align(center, sequences)

if __name__ == '__main__':
    main()
