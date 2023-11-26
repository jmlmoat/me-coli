#!/usr/bin/env python

import numpy as np
import os, sys
#from sklearn.feature_selection import SelectKBest
import argparse
import math

# stop numpy from being a greedy bastard
import os
os.environ["OMP_NUM_THREADS"] = "1"

if __name__== "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-ms', help='total num of split', required=True)
    parser.add_argument('-uk', help='union_kmers', required=True)
    parser.add_argument('-o', help='outdir', required=True)
    parser.add_argument('-g', help='genome names', required=True)
    parser.add_argument('-i', nargs='+', required=True)
    parser.add_argument('-r', nargs='+', required=True)
    parser.add_argument('-c', nargs='+', required=True)
    args = parser.parse_args()

    infiles = args.i
    rowfiles = args.r
    colfiles = args.c

    max_splits = args.ms
    union_path = args.uk
    outdir = args.o
    genome_path = args.g

    # load genomes and union kmers
    all_genomes = np.load(genome_path,mmap_mode='r')
    features = np.load(union_path,mmap_mode='r')

    # matrix dimensions
    num_rows = len(all_genomes)
    num_cols = len(features)

    #matrix_dimensions = ( (num_rows,num_cols), dtype=uint8 )

    #full_matrix = []
    full_matrix = np.zeros( (num_rows,num_cols), dtype=np.uint8 )
    print(full_matrix)
    print()

    cols_dict = { features[i] : i for i in range(0, num_cols)}
    rows_dict = { all_genomes[i] : i for i in range(0, num_rows)}

    # current split number
    split_num = 1
    update_row_num = 0
    for infile,rowfile,colfile in zip(infiles,rowfiles,colfiles):

        import os, psutil
        print("------------------------------------------------")
        print(split_num)
        process = psutil.Process(os.getpid())
        print("rss",process.memory_info().rss)
        print("mem",process.memory_info(), flush=True)

        #split_matrix = np.load(infile)
        #split_rows_genomes = np.load(rowfile)

        split_matrix = np.load(infile,mmap_mode='r')
        split_rows_genomes = np.load(rowfile,mmap_mode='r')

        split_num_rows = len(split_rows_genomes)

        #with np.load(infile,mmap_mode='r') as split_matrix, \
        #np.load(rowfile,mmap_mode='r') as split_rows_genomes:

        num_per_split = math.ceil(len(all_genomes)/int(max_splits))
        #start = int(split_num)*num_per_split-1
        #stop = start+num_per_split

        for i, genome_id in enumerate(split_rows_genomes):
            # Assert that the rows are being added to the main matrix in the correct order
            assert genome_id == all_genomes[i+(split_num-1)*num_per_split]


        #print(start)
        #print(stop)

        start = ( (int(split_num)-1) * num_per_split )
        stop = start + split_num_rows

        print(start)
        print(stop)

        print()

        full_matrix[ start:stop, :] = split_matrix


        #for row in split_matrix:
        #    full_matrix[update_row_num] = row
        #    update_row_num += 1

        #if split_num == 1:
            #full_matrix = split_matrix
        #else:
            #full_matrix = np.concatenate((full_matrix,split_matrix), axis=0)

        split_num += 1

        #print(full_matrix)
        #sys.exit()

    #sys.exit()

    np.save("{}matrix.npy".format(outdir), full_matrix)
    np.save("{}rows.npy".format(outdir), all_genomes)
    np.save("{}cols.npy".format(outdir), features)
