
class Node:
    def __init__(self, data):
        self.data = data
        self.children = []

    def add(self, *children):
        for child in children:
            self.children.append(child)

    def traverse(self):
        print(self.data)
        for child in self.children:
            child.traverse()


def score(c1, c2, distances, clusters):
    u1 = sum(distances[c1, ck] for ck in clusters if ck != c1)
    u2 = sum(distances[c2, ck] for ck in clusters if ck != c2)

    return (len(clusters)-2)*distances[c1, c2] - u1 - u2


def construct_tree(D, sequences):

    clusters = [Node(sequence) for sequence in sequences]
    distances = {(c1, c2): D[c1.data, c2.data] for c1 in clusters
                                               for c2 in clusters}

    while len(clusters) > 1:
        cx, cy = min(((c1, c2) for c1 in clusters for c2 in clusters if c1 != c2),
                     key=lambda x: score(x[0], x[1], distances, clusters))
        new_cluster = Node(cx.data + cy.data)
        new_cluster.add(cx, cy)

        clusters.remove(cx)
        clusters.remove(cy)

        for cz in clusters:
            new_distance = (distances[cx, cz] + distances[cy, cz]
                            - distances[cx, cy]) / 2
            distances[new_cluster, cz] = new_distance
            distances[cz, new_cluster] = new_distance

        clusters.append(new_cluster)

    clusters[0].traverse()


sequences = ['A', 'B', 'C', 'D', 'E', 'F']
M = [[0, 5, 4, 7, 6, 8],
     [5, 0, 7, 10, 9, 11],
     [4, 7, 0, 7, 6, 8],
     [7, 10, 7, 0, 5, 9],
     [6, 9, 6, 5, 0, 8],
     [8, 11, 8, 9, 8, 0]]

D = {(sequences[i], sequences[j]): M[i][j] for i in range(6) for j in range(6)}

construct_tree(D, sequences)
