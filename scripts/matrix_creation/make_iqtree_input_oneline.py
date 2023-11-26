#!/usr/bin/env python

# Basics
import pandas as pd
import numpy as np
import argparse
import os, sys, gc


'''
python scripts/matrix_creation/make_iqtree_input.py
sbatch -c 1 --mem 22G -o iq.out -J IQformat --wrap="python scripts/matrix_creation/make_iqtree_input_oneline.py"
'''

if __name__ == "__main__":
    # Get command line arguments
    '''
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--rows', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    '''

    kmer_len = "11" #"11"

    matrix_path = "results/wgs_standard/{}mer/matrix/matrix.npy".format(kmer_len)
    rows_path = "results/wgs_standard/{}mer/matrix/rows.npy".format(kmer_len)
    outdir = "results/wgs_standard/{}mer/iqtree/".format(kmer_len)
    #outfile_tmp = "{}binary_kmer.txt".format(outdir)
    outfile = "{}binary_kmer_alignment_oneline.phy".format(outdir)

    # if outdir doesnt exist, make it
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)


    # load matrix
    matrix_freq = np.load(matrix_path, allow_pickle=True,mmap_mode='r')
    # load genome names
    genomes = np.load(rows_path, allow_pickle=True, mmap_mode='r')

    # convert matrix from frequency to binary presence/absence
    matrix_binary = np.array(matrix_freq, dtype=bool) # True/False
    matrix_binary = matrix_binary.astype(np.uint8) # 1/0

    # Check matrix_dimensions
    print("Binary Matrix uint8")
    print("  num genomes: {}".format(len(genomes)))
    print("  matrix shape: {}".format(matrix_binary.shape))
    print("  matrix dtype: {}".format(matrix_binary.dtype))

    print(matrix_binary)

    # Get dimensions
    n_genomes = matrix_binary.shape[0]
    n_kmers = matrix_binary.shape[1]

    print(n_genomes,n_kmers)

    print("Double checking dimensions")
    print("  matrix n_genomes: {}".format(n_genomes))
    print("  matrix n_kmers: {}".format(n_kmers), flush=True)

    print("starting save",flush=True)

    header = "{} {}".format(n_genomes,n_kmers)
    # dont change the number of features to 1 when making single string!!!!
    # http://www.iqtree.org/doc/Tutorial
    #header = "{} 1".format(n_genomes,n_kmers)

    np.savetxt(outfile,
               matrix_binary,
               fmt='%d',
               delimiter='',
               header=header,
               comments='')
    # delimiter = ' ' space between 0 1
    # delimiter = '' no space

    print("done save",flush=True)

    print("starting add genome",flush=True)

    # insert genome name onto each line
    import fileinput

    for index,line in enumerate(fileinput.input([outfile], inplace=True)):

        # first line/index of file is the header
        if index!=0:
            genome = genomes[index-1]
            sys.stdout.write('{g} {l}'.format(g=genome,l=line))
        else:
            sys.stdout.write('{n} {m}\n'.format(n=n_genomes,m=n_kmers))

    print("done add genome",flush=True)


    # Make an output file, PHYLIP format
    '''
    PHYLIP format

    References:
    http://www.iqtree.org/doc/Tutorial#binary-morphological-and-snp-data
    http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html

    n_rows n_cols
    genome1 0 1 1 0 1 ...
    genome2 1 1 1 1 0 ...
    ...

    '''
