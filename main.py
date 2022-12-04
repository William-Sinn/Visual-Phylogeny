from Bio import SeqIO
from scipy import spatial
import hammingdist


def matrix_gen(seq_file):
    mat = []
    row_index = 0

    for record_y in SeqIO.parse(seq_file, "fasta"):
        line = []
        col_index = 0
        for record_x in SeqIO.parse(seq_file, "fasta"):
            if row_index >= col_index:
                line.append(spatial.distance.hamming(record_x, record_y) * len(record_x))
                col_index += 1

            else:
                continue

        mat.append(line)
        row_index += 1

    return mat


def neighbor_joining_matrix_gen(matrix):
    nodes = len(matrix)
    nj_matrix = []
    row_index = 0

    for row in matrix:
        col_index = 0
        distance_row = []

        for index in row:

            if row_index > col_index:
                nodes_not_neighbors = nodes - 2
                neighbors_distance = nodes_not_neighbors * index

                total_dist_i = sum(matrix[row_index])
                total_dist_j = sum(matrix[col_index])

                nj_index = neighbors_distance - total_dist_i - total_dist_j
                distance_row.append(nj_index)

                col_index += 1

            else:
                continue

        nj_matrix.append(distance_row)
        row_index += 1

    return nj_matrix


def print_matrix(matrix):
    for row in matrix:
        for item in row:
            print(f'{item:<8} ', end="")
        print("")
    print("")


if __name__ == '__main__':
    in_matrix = matrix_gen("input_seqs/testin.fasta")

    sample_matrix = [
        [0, 13, 21, 22],
        [13, 0, 12, 13],
        [21, 12, 0, 13],
        [22, 13, 13, 0]
    ]

    nj_mat = neighbor_joining_matrix_gen(sample_matrix)

    print_matrix(in_matrix)
    print_matrix(sample_matrix)
    print_matrix(nj_mat)
