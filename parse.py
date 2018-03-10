from re import search


def parse_fasta(path):
    """
    Extract a list of sequences from a .fasta file.
    :param path: the path to the file
    :return: a dictionary with species name -> sequence
    """
    try:
        with open(path) as f:
            lines = [line.strip().upper() for line in f.readlines()]
    except FileNotFoundError:
        print('Error: no such file: {}'.format(path))
        exit(1)

    sequences = {}
    current = ''
    name = ''
    for line in lines:
        if line.startswith('>'):
            if len(current) > 0:
                sequences[name] = current
            current = ''
            name = search('\[(.*)\]', line).group(1)
            continue

        current += line

    if len(current) > 0:
        sequences[name] = current

    return sequences

