#!/usr/bin/env python

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
import sys
import argparse
import math
from Bio.SeqIO.FastaIO import SimpleFastaParser

# stop numpy from being a greedy bastard
import os
os.environ["OMP_NUM_THREADS"] = "1"

def make_row(filename):
    """
    Given a genome file, create and return a row of kmer counts
    to be inserted into the kmer matrix.
    """

    #****************************************************************************
    #import os, psutil
    #process = psutil.Process(os.getpid())
    #rss = process.memory_info().rss
    #mem = process.memory_info()
    #print("----------------------\n{}\n{}\n".format(rss,mem,flush=True))
    #****************************************************************************

    # Don't reload this, waste of resources, it's loaded in __main__
    # and this fn shouldn't be called by anything else.
    #relevant_feats = np.load(mm_path)

    # Don't redeclare this, again waste of resources, it's in __main__
    #cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*len(relevant_feats)

    # lookup dataset
    #dataset = master[master['id'] == filename]['dataset'].values[0]
    #jf_path = jelly+"{}/{}.fa".format(dataset,filename)
    jf_path = jelly+"{}.fa".format(filename)

    # parse each file into a row
    with open(jf_path) as fasta_file:
        for kmer_count, kmer_sequence in SimpleFastaParser(fasta_file):
            # ensure the count is an int
            kmer_count = int(kmer_count)
            # if the kmer count (title) is over 255, it must be rounded down
            # so it can be stored as uint8
            if kmer_count > 255:
                kmer_count = 255
            temp_row[cols_dict[kmer_sequence]] = kmer_count
    return filename, temp_row

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-uk', help='union_kmers', required=True)
    parser.add_argument('-o', help='output', required=True)
    parser.add_argument('-d', help='matrix dtype', required=True)
    parser.add_argument('-s', help='split_num', required=True)
    parser.add_argument('-ms', help='max splits', required=True)
    parser.add_argument('-i', help='input_path', required=True)
    parser.add_argument('-g', help='genome names', required=True)
    parser.add_argument('-m', help='master df', required=True)

    args = parser.parse_args()

    jelly = args.i
    matrix_dtype = args.d
    mm_path = args.uk
    split_num = args.s
    max_splits = args.ms
    genome_name_path = args.g
    outdir = args.o
    master = args.m
    master = pd.read_csv(master,sep='\t')

    num_start = 0
    num_stop = 0
    total = 0

    def progress():
        sys.stdout.write('\r')
        sys.stdout.write("Loading Genomes: {} started, {} finished, {} total".format(num_start,num_stop,total))
        #sys.stdout.flush()
        if(num_stop==total):
            print("\nAll Genomes Loaded!\n")

    relevant_feats = np.load(mm_path,mmap_mode='r')


    # Load the list of ids; this will be the row order in the matrix
    genomes = np.load(genome_name_path,mmap_mode='r')
    # How many genomes will be in each split
    num_per_split = math.ceil(len(genomes)/(int(max_splits)))
    # How many genomes will be in the current split
    start = (int(split_num)-1)*num_per_split
    stop = start+num_per_split
    genomes = genomes[start:stop]


    total = len(genomes)

    runs = [i.split('.')[0] for i in genomes]

    # declaring empty kmer matrix to fill
    kmer_matrix = np.zeros((len(genomes),len(relevant_feats)),dtype = matrix_dtype)

    # making dicts for faster indexing
    # note that rows dict is in filenames not genome/run names
    rows_dict = { genomes[i] : i for i in range(0, len(genomes))}
    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Use concurrent futures to get multiple rows at the same time
    # Then place completed rows into the matrix and update the row dictionary
    #num_start += min(16,len(genomes))
    num_start += min(16,len(genomes)) # ran out of mem; trying 8
    progress()
    with ProcessPoolExecutor(max_workers=16) as ppe:
        for genome_name,temp_row in ppe.map(make_row, genomes):
            num_stop+=1
            if(num_start<total):
                num_start+=1
            progress()
            for i, val in enumerate(temp_row):
                kmer_matrix[rows_dict[genome_name]][i] = val

    # save everything
    np.save(outdir+"matrix_split_{}.npy".format(split_num), kmer_matrix)
    np.save(outdir+"matrix_split_{}_rows.npy".format(split_num), runs)
    np.save(outdir+"matrix_split_{}_cols.npy".format(split_num), relevant_feats)

    # testing df version
    #df = pd.DataFrame(data=kmer_matrix, columns=relevant_feats, index=runs)
    #df.to_csv("{}matrix_split_{}.tsv".format(outdir,split_num),index=True,sep='\t')
