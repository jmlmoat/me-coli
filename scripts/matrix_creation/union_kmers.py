#!/usr/bin/env python

import os, sys
#from Bio import Seq, SeqIO
import numpy as np
import time
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import chain
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
import math
import argparse

# stop numpy from being a greedy bastard
import os
os.environ["OMP_NUM_THREADS"] = "1"

def parse_fasta(oth_args):
    genome, kmer_length = oth_args
    current_multi_mers = []
    # parse each file
    with open(genome) as fasta_file:
        for title, sequence in SimpleFastaParser(fasta_file):
            if(len(sequence) == int(kmer_length)):
                current_multi_mers.append(sequence)
    return current_multi_mers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-k', help='kmer_length', required=True)
    parser.add_argument('-o', help='output_path', required=True)
    parser.add_argument('-s', help='set_num', required=True)
    parser.add_argument('-ms', help='max_sets', required=True)
    parser.add_argument('-i', nargs='+', required=True)

    args = parser.parse_args()

    kmer_length = args.k
    out = args.o
    set_num = int(args.s)-1
    max_sets = args.ms

    genomes = args.i

    num_per_split = math.ceil(len(genomes)/int(max_sets))
    start_of_split = int(set_num)*num_per_split
    end_of_split = start_of_split+num_per_split

    genomes = genomes[start_of_split:end_of_split]

    multi_mers = []

    print("Starting parse")
    start = time.time()
    genome_counter = 0
    times_printed = 0
    with ProcessPoolExecutor(max_workers = min(16, cpu_count())) as ppe:
        for current_multi_mer in ppe.map(parse_fasta, zip(genomes,[kmer_length for i in range(len(genomes))])):
            genome_counter+=1
            multi_mers.append(current_multi_mer)
            if(genome_counter>100):
                times_printed += 1
                print("done:", times_printed*100)
                genome_counter -= 100


    print("Parse took {}s, starting unions".format(time.time()-start))
    #start = time.time()
    #master_mers = list(set().union(*multi_mers))
    #print("Union Time:",time.time()-start)
    start = time.time()
    master_mers = np.array(list(set(chain(*multi_mers))))
    bool_mask = [len(master_mers[i]) for i in range(len(master_mers))]
    bool_mask = np.array([i == int(kmer_length) for i in bool_mask])
    master_mers = master_mers[bool_mask]
    print("Chain Time:",time.time()-start)

    print("Found {} unique multi-mers".format(len(master_mers)))

    np.save(out.split('.')[0]+'.npy', master_mers)
