from Bio import SeqIO
import numpy as np
from scipy import spatial


if __name__ == '__main__':
    x = 0
    y = 0
    mat = np.empty(([33, 33]))

    # print(records1)
    for record_y in SeqIO.parse("CompleteMTGenomes_MAFFT_STRPD.FASTA", "fasta"):
        x = 0
        print(record_y.id)
        for record_x in SeqIO.parse("CompleteMTGenomes_MAFFT_STRPD.FASTA", "fasta"):
            mat[x, y] = spatial.distance.hamming(record_x, record_y) * len(record_x)
            x += 1
        y += 1

    print(mat)