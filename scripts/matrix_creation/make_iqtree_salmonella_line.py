import numpy as np
import pandas as pd
import sys
import math
from Bio.SeqIO.FastaIO import SimpleFastaParser




if __name__ == "__main__":
    """
    python scripts/matrix_creation/make_iqtree_salmonella_line.py
    sbatch -c 1 --mem 1G -o iq.out --wrap="python scripts/matrix_creation/make_iqtree_salmonella_line.py"
    """


    outdir = "salmonella_outgroup/"


    # path to the salmonella jellyfish output
    jf_path = "salmonella_outgroup/salmonella.fa"


    # path to the matrix's columns (kmers in order)
    mm_path = "results/wgs_standard/11mer/matrix/cols.npy"
    # list of kmers in order
    relevant_feats = np.load(mm_path,mmap_mode='r')
    # dict linking kmer name to index
    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    '''
    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*len(relevant_feats)

    # parse each file into a row
    with open(jf_path) as fasta_file:
        for kmer_count, kmer_sequence in SimpleFastaParser(fasta_file):
            if kmer_sequence in relevant_feats:
                # ensure the count is an int
                kmer_count = int(kmer_count)
                # if the kmer count (title) is over 255, it must be rounded down
                # so it can be stored as uint8
                if kmer_count > 255:
                    kmer_count = 255
                temp_row[cols_dict[kmer_sequence]] = kmer_count

    # save kmer count row
    row = np.array(temp_row)
    np.save("{}kmer_count.npy".format(outdir), row)

    print(row)
    '''

    row = np.load("{}kmer_count.npy".format(outdir))

    # binary presence absence
    # convert matrix from frequency to binary presence/absence
    row_binary = np.array(row, dtype=bool) # True/False
    row_binary = row_binary.astype(np.uint8) # 1/0
    #save
    np.save("{}kmer_count_binary.npy".format(outdir), row_binary)

    print(row_binary)

    # need to be in a format like [[0, 1, 0]]
    row_binary = row_binary.reshape(1, row_binary.shape[0])

    outfile = "{}kmer_count_binary.phy".format(outdir)
    # save binary as txt file
    np.savetxt(outfile,
               row_binary,
               fmt='%d',
               delimiter='',
               comments='')

    # insert genome name onto each line
    import fileinput

    for index,line in enumerate(fileinput.input([outfile], inplace=True)):
        genome = "Salmonella"
        sys.stdout.write('{g} {l}'.format(g=genome,l=line))
