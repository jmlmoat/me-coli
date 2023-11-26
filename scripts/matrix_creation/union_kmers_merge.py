#!/usr/bin/env python
"""
takes the 5 folds of seen kmers and finds a master list
"""

import numpy as np
import pandas as pd
import os, sys
from itertools import chain
import argparse

# stop numpy from being a greedy bastard
import os
os.environ["OMP_NUM_THREADS"] = "1"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-k', help='kmer_length', required=True)
    parser.add_argument('-o', help='output_path', required=True)
    parser.add_argument('-i', nargs='+', required=True)
    args = parser.parse_args()

    kmer_length = args.k
    out = args.o
    infiles = args.i

    set_nums = [str(i) for i in range(0,len(infiles))]#args.set

    out_folder = out.rsplit('/',1)[0]

    feat_sets = {}
    for set_num in set_nums:
        #feat_sets[set_num] = np.load(out_folder+"/{}mers_{}.npy".format(kmer_length,set_num))
        # I didn't try mmap_mode='r' for this, not sure if it would have an
        # effect because the whole thing is being used
        feat_sets[set_num] = np.load(infiles[int(set_num)])
        '''
        # why doesnt context manager work with np.load anymore?
        file = infiles[int(set_num)]
        print("*******************************")
        print(file)
        print(type(file))
        y = np.load(file)
        print(y)
        print(np.__version__)
        with (np.load(file)) as x:
            print(x)
            feat_sets[set_num] = x
        sys.exit()
        '''

    all_feats = [feat_sets[i] for i in set_nums]

    master_mers = np.array(list(set(chain(*all_feats))))
    bool_mask = [len(master_mers[i]) for i in range(len(master_mers))]
    bool_mask = np.array([i == int(kmer_length) for i in bool_mask])
    master_mers = master_mers[bool_mask]

    np.save(out, master_mers)
