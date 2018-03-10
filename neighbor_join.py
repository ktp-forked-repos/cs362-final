
class Node:
    def __init__(self, label):
        self.label = label
        self.children = []

    def add(self, *children):
        for child in children:
            self.children.append(child)

    def traverse(self):
        print(self.label)
        for child in self.children:
            child.traverse()

    def __hash__(self):
        return hash(self.label)


def score(c1, c2, distances, clusters):
    u1 = sum(distances[c1, ck] for ck in clusters if ck != c1)
    u2 = sum(distances[c2, ck] for ck in clusters if ck != c2)

    return (len(clusters)-2)*distances[c1, c2] - u1 - u2


def write_dot(root):
    """
    Write a graph to a dot file.
    :param root: the root of the tree
    :return: nothing
    """

    nodes = ''
    edges = ''
    stack = [root]
    while len(stack) > 0:
        node = stack.pop()
        nodes += '"{}" [label="{}"];\n'.format(hash(node), node.label if not ',' in node.label else '')
        for child in node.children:
            edges += '"{}"->"{}";\n'.format(hash(node.label), hash(child.label))
            stack.append(child)

    out = 'digraph mygraph {\n' + nodes + edges + '\n}'

    with open('out.dot', 'w') as f:
        f.write(out)


def construct_tree(D, sequences):

    clusters = [Node(name) for name in sequences]
    distances = {(c1, c2): D[c1.label, c2.label] for c1 in clusters
                                               for c2 in clusters}

    while len(clusters) > 1:
        cx, cy = min(((c1, c2) for c1 in clusters for c2 in clusters if c1 != c2),
                     key=lambda x: score(x[0], x[1], distances, clusters))
        new_cluster = Node(cx.label + ',' + cy.label)
        new_cluster.add(cx, cy)

        clusters.remove(cx)
        clusters.remove(cy)

        for cz in clusters:
            new_distance = (distances[cx, cz] + distances[cy, cz]
                            - distances[cx, cy]) / 2
            distances[new_cluster, cz] = new_distance
            distances[cz, new_cluster] = new_distance

        clusters.append(new_cluster)

    write_dot(clusters[0])

    return clusters[0]


# sequences = ['A', 'B', 'C', 'D', 'E', 'F']
# M = [[0, 5, 4, 7, 6, 8],
#      [5, 0, 7, 10, 9, 11],
#      [4, 7, 0, 7, 6, 8],
#      [7, 10, 7, 0, 5, 9],
#      [6, 9, 6, 5, 0, 8],
#      [8, 11, 8, 9, 8, 0]]
#
# D = {(sequences[i], sequences[j]): M[i][j] for i in range(6) for j in range(6)}
#
# construct_tree(D, sequences)
