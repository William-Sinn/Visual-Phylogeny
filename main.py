from Bio import SeqIO
from scipy import spatial
import hammingdist


class Species:
    def __init__(self, fasta_record=None):
        if fasta_record is not None:
            self.name = fasta_record.id
            self.description = fasta_record.description
            self.seq = fasta_record.seq

    def __str__(self):
        return self.name

    def __format__(self, format_spec):
        return format(self.name, format_spec)


class Cluster(Species):
    def __init__(self, name, fasta_record=None):
        super(Cluster, self).__init__(fasta_record)
        self.name = name


def matrix_gen(seq_file):
    mat = []
    row_index = 0
    species_dict = {}

    for record_y in SeqIO.parse(seq_file, "fasta"):
        line = []
        species_dict[row_index] = Species(record_y)
        col_index = 0

        for record_x in SeqIO.parse(seq_file, "fasta"):
            line.append(spatial.distance.hamming(record_x, record_y) * len(record_x))
            col_index += 1

        mat.append(line)
        row_index += 1

    return mat, species_dict


def print_matrix(matrix, species_dict=None):
    row_index = 0

    for row in matrix:
        (lambda: None if type(species_dict) is not dict
        else print(f"{species_dict[row_index]:<12}", end=""))()

        for item in row:
            print(f'{item}'.center(10),
                  end="")

        row_index += 1
        print("")
    print("")


def get_value_llt(matrix, x, y):
    if y <= len(matrix[x]) - 1:
        return matrix[x][y]
    else:
        return matrix[y][x]


def get_row_dist_llt(matrix, row_index):
    total = 0
    row_len = len(matrix)

    for col_index in range(0, row_len):
        total += get_value_llt(matrix, row_index, col_index)

    return total


def get_dists_matrix(matrix):
    dists_list = {}
    for row in range(len(matrix)):
        dists_list[row] = get_row_dist_llt(matrix, row)

    return dists_list


def neighbor_joining_matrix_gen(matrix):
    nodes = len(matrix)
    nj_matrix = []
    row_index = 0
    dist_dict = get_dists_matrix(matrix)

    for row in matrix:
        col_index = 0
        distance_row = []

        for index in row:

            if row_index > col_index:
                nodes_not_neighbors = nodes - 2
                neighbors_distance = nodes_not_neighbors * index

                total_dist_i = dist_dict[row_index]
                total_dist_j = dist_dict[col_index]

                nj_index = neighbors_distance - total_dist_i - total_dist_j
                distance_row.append(nj_index)

                col_index += 1

            elif col_index == row_index:
                nj_index = 0
                distance_row.append(nj_index)

                col_index += 1

            else:
                continue

        nj_matrix.append(distance_row)
        row_index += 1

    return nj_matrix


# for testing only!!!for testing only!!!for testing only!!!for testing only!!!for testing only!!!
def recon(matrix):
    print_matrix(matrix)
    new_mat = neighbor_joining_matrix_gen(matrix)

    for row in range(len(new_mat)):
        for i in range(len(matrix)):
            print(get_value_llt(matrix, row, i), end=" ")
        print("")


# for testing only!!!for testing only!!!for testing only!!!for testing only!!!for testing only!!!


def neighbor_joining_phylogeny(matrix, species_dict):
    matrix_len = len(matrix)

    if matrix_len == 2:
        keys = list(species_dict.keys())
        return f"({species_dict[keys[0]].name}, {species_dict[keys[1]].name})"

    matrix_prime = neighbor_joining_matrix_gen(matrix)
    min_distance = (0, 0)

    row_index = 0
    for row in matrix_prime:
        col_index = 0
        for col in row:
            if col < matrix_prime[min_distance[0]][min_distance[1]]:
                min_distance = (row_index, col_index)
            col_index += 1
        row_index += 1

    # delta_ij = (matrix_dists[min_distance[1]] - matrix_dists[min_distance[0]]) / (matrix_len - 2)
    #
    # limb_length = (lambda i, j, plus: (matrix[i][j] + delta_ij) / 2 if plus
    # else (matrix[i][j] - delta_ij) / 2)
    #
    # limb_len_i = limb_length(min_distance[0], min_distance[1], True)
    # limb_len_j = limb_length(min_distance[0], min_distance[1], False)

    group = f"({species_dict[min_distance[0]].name}, {species_dict[min_distance[1]].name})"

    new_matrix = [[0]]

    new_cluster = Cluster(group)
    new_dict = {0: new_cluster}
    new_row_index = 1

    for row in range(matrix_len):
        if row not in min_distance:
            new_dict[new_row_index] = species_dict[row]
            new_row = [(get_value_llt(matrix, row, min_distance[0])
                        + get_value_llt(matrix, row, min_distance[1])
                        - get_value_llt(matrix, min_distance[0], min_distance[1])) / 2]
            new_row_index += 1
            for col in range(matrix_len):
                if col <= row and col not in min_distance:
                    new_row.append(matrix[row][col])

            new_matrix.append(new_row)
    # print_matrix(matrix)
    # print_matrix(matrix_prime)
    # print(min_distance)
    # print(group)
    # print(matrix_prime[min_distance[0]][min_distance[1]])
    return neighbor_joining_phylogeny(new_matrix, new_dict)


if __name__ == '__main__':
    in_matrix, s_dict = matrix_gen("input_seqs/testin.fasta")
    simple_in, s_dict2 = matrix_gen("input_seqs/sample.fasta")

    sample_matrix = [
        [0, 13, 21, 22],
        [13, 0, 12, 13],
        [21, 12, 0, 13],
        [22, 13, 13, 0]
    ]
    sample_dict = {0: Cluster("i"), 1: Cluster("j"), 2: Cluster("k"), 3: Cluster("l")}

    # print_matrix(sample_matrix)
    # nj_mat = neighbor_joining_matrix_gen(sample_matrix)
    # print_matrix(nj_mat, sample_dict)
    # print(neighbor_joining_phylogeny(sample_matrix, sample_dict))
    # print()
    #
    # print_matrix(in_matrix)
    # nj_mat2 = neighbor_joining_matrix_gen(in_matrix)
    # print_matrix(nj_mat2, s_dict)
    # print("output")
    print(neighbor_joining_phylogeny(in_matrix, s_dict))
    # print()
    #
    # print_matrix(simple_in)
    # nj_mat3 = neighbor_joining_matrix_gen(simple_in)
    # print_matrix(nj_mat3, s_dict2)
    print(neighbor_joining_phylogeny(simple_in, s_dict2))
    # print()

    ch5_mat, ch5_dict = matrix_gen("input_seqs/Chap5_mtDNA_FASTA.fasta")
    print(neighbor_joining_phylogeny(ch5_mat, ch5_dict))

    ch5_mat, ch5_dict = matrix_gen("input_seqs/CompleteMTGenomes_MAFFT_STRPD.FASTA")
    print(neighbor_joining_phylogeny(ch5_mat, ch5_dict))

    # recon(in_matrix)

