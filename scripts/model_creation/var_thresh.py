#!/usr/bin/env python

# Basics
import pandas as pd
import numpy as np
import argparse
import pickle
import os, sys, gc
import joblib
from joblib import dump, load
from ruamel.yaml import YAML
import psutil
import time

# Helper scripts for grabbing data etc
from model_helper import *

# Filter out the constant columns (eg all 0s etc)
from sklearn.feature_selection import VarianceThreshold



if __name__ == "__main__":
    # Get command line arguments
    #parser.add_argument('--rows', required=True)
    #parser.add_argument('--labels', required=True)
    #parser.add_argument('--modelname', required=True)
    #parser.add_argument('--nsplits',required=True)

    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--cols', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    #model_name = args.modelname.lower()
    matrix = args.matrix
    #rows = args.rows
    #labels = args.labels
    cols = args.cols
    outdir = args.outdir
    #nsplits = int(args.nsplits)

    # --------------------------------------------------------------------------
    # Set up Directories

    # if outdir doesnt exist, make it
    if not os.path.exists(outdir):
        os.makedirs(out_path, exist_ok=True)

    # folder for the split pieces
    outdir_tmp = "{}thresh_chunks/".format(outdir)
    if not os.path.exists(outdir_tmp):
        os.makedirs(outdir_tmp, exist_ok=True)

    # --------------------------------------------------------------------------
    # Load Data
    print("\nStarting: load data", flush=True)

    # load data
    x = np.load(matrix, allow_pickle=True,mmap_mode='r')
    #labels = np.load(labels, allow_pickle=True)
    #rows = np.load(rows, allow_pickle=True,mmap_mode='r')
    k = np.load(cols, allow_pickle=True, mmap_mode='r')

    print("  x shape", x.shape)
    print("  k shape", k.shape)


    # --------------------------------------------------------------------------
    # Split
    print("\nStarting Split", flush=True)

    nsplits = 5

    # Define the split locations
    # Can't just tell np.hsplit to make 5 splits etc, because it requires
    # the splits to be perfectly equal

    num_kmers = x.shape[1]
    step_size = int(num_kmers / nsplits)
    split_starts = [step_size*i for i in range(nsplits+1)]
    split_ends = [step_size*i for i in range(1,nsplits+1)]
    # make the last value the max value (num kmers) to include any remainders
    split_ends[-1] = num_kmers
    print()
    print("  num columns (kmers): {}".format(num_kmers))
    print("  step size: {}".format(step_size))
    print("  split markers - start: {}".format(split_starts))
    print("  split markers - end: {}".format(split_ends))


    # Split into pieces column-wise; result is a list of numpy matrices
    #x = np.hsplit(x, nsplits)
    x = np.hsplit(x, split_ends)
    #k = np.hsplit(k, nsplits)
    k = np.hsplit(k, split_ends)

    print()
    print("  len x:", len(x))
    print("  x[0] shape:", x[0].shape)
    print("  len k:", len(k))
    print("  k[0] shape:", k[0].shape)
    #sys.exit()

    '''
    print()
    for s in x:
        print(s)
    for q in k:
        print(q)
    '''

    # --------------------------------------------------------------------------
    # Remove columns containing all zeros
    print("\nStarting Column Masking", flush=True)

    for i in range(nsplits):
        print("  Masking split # {}".format(i),flush=True)

        x_tmp = x[i]
        k_tmp = k[i]

        num_cols = x_tmp.shape[1]
        num_cols_2 = len(k_tmp)
        print("    num cols x and k:", num_cols, num_cols_2)


        print()
        print("    getting indices of all 0 cols",flush=True)
        start = time.time()
        # array of indices for columns that contain all zeroes
        zero_col_index_np = np.where(~x_tmp.any(axis=0)) #[0]
        zero_col_index = zero_col_index_np[0]
        end = time.time()
        print("      time: ", end-start)
        # convert to list
        zero_col_index = zero_col_index.tolist()
        print("      num of zero columns", len(zero_col_index))
        #print("      zero col indices", type(zero_col_index))

        # instead of mask just use delete

        # apply column mask to matrix and columns
        print("    applying mask", flush=True)

        print("      k before mask ",k_tmp.shape, flush=True)
        start = time.time()
        k_tmp = np.delete(k_tmp,zero_col_index)
        end = time.time()
        print("      k after mask  ",k_tmp.shape, flush=True)
        print("        time: ", end-start)

        col_mask = np.zeros(num_cols, dtype=bool)
        col_mask[zero_col_index_np] = 1
        reverse_col_mask = ~col_mask

        print("      x before mask ",x_tmp.shape, flush=True)
        start = time.time()
        #x_tmp = np.delete(x_tmp,zero_col_index,axis=1)
        x_tmp = x_tmp[:,reverse_col_mask]
        end = time.time()
        print("      x after mask  ",x_tmp.shape, flush=True)
        print("        time: ", end-start)

        #sys.exit()

        # Saving
        print("  Saving split # {}".format(i),flush=True)
        start = time.time()
        np.save("{}matrix{}_filtered.npy".format(outdir_tmp,i), x_tmp)
        np.save("{}cols{}_filtered.npy".format(outdir_tmp,i), k_tmp)
        end = time.time()
        print("    time: ", end-start)

        del x_tmp, k_tmp, zero_col_index, zero_col_index_np, col_mask, reverse_col_mask
        gc.collect()

    del i
    gc.collect()

    #sys.exit()


    # --------------------------------------------------------------------------
    # Merge
    print("\nStarting Merge", flush=True)

    x_new = []
    k_new = []

    start = time.time()
    # Load a piece of a time and merge it
    for i in range(nsplits):
        print("  Merging split # {}".format(i),flush=True)

        if i == 0:
            x_new = np.load("{}matrix{}_filtered.npy".format(outdir_tmp,i))
        else:
            x_train_tmp = np.load("{}matrix{}_filtered.npy".format(outdir_tmp,i))
            x_new = np.concatenate((x_new,x_train_tmp),axis=1)

            del x_train_tmp
            gc.collect()

        k_tmp = np.load("{}cols{}_filtered.npy".format(outdir_tmp,i))
        for z in k_tmp:
            k_new.append(z)
    end = time.time()

    x_new = np.array(x_new, dtype=np.uint8)
    #k_new = np.array(k_new, dtype=np.uint8)

    print("  x shape:", x_new.shape)
    print("  k length:", len(k_new))
    print("  time:", end-start)

    np.save("{}matrix_vt.npy".format(outdir,i), x_new)
    np.save("{}cols_vt.npy".format(outdir,i), k_new)





    '''
    # --------------------------------------------------------------------------
    # Try np mask on whole matrix
    print("Starting: np mask on whole matrix", flush=True)

    # array of indices for columns that contain all zeroes
    zero_col_index = np.where(~x.any(axis=0))[0]
    print(zero_col_index)

    # get number of columns present, need for mask length
    num_cols = x.shape[1]
    num_cols_2 = len(k)
    print(num_cols,num_cols_2,flush=True)

    # create a mask with the column indices
    col_mask = np.zeros(num_cols, dtype=bool)
    col_mask[zero_col_index] = 1
    print(col_mask,flush=True)

    # apply column mask to matrix and columns

    print("k ", k.shape, flush=True)
    k = k[:,col_mask.astype('bool')]
    print("k ", k.shape, flush=True)

    sys.exit()

    print("x ", x.shape)
    x = x[:,col_mask.astype('bool')]
    print("x ", x.shape)

    np.save("{}matrix_vt.npy".format(out_dir,i), x)
    np.save("{}cols_vt.npy".format(out_dir,i), k)
    '''

    '''
    # --------------------------------------------------------------------------
    # Save the splits
    for i in range(nsplits):
        print()
        print("  Saving split # {}".format(i),flush=True)

        np.save("{}matrix{}_not_filtered.npy".format(outdir_tmp,i), x[i])
        np.save("{}cols{}_not_filtered.npy".format(outdir_tmp,i), k[i])

    del x, k, i
    gc.collect()


    # --------------------------------------------------------------------------
    # Remove columns containing all zeros
    print("\nStarting Column Masking", flush=True)

    for i in range(nsplits):
        print("  Saving split # {}".format(i),flush=True)

        x = np.load("{}matrix{}_not_filtered.npy".format(outdir_tmp,i))
        k = np.load("{}cols{}_not_filtered.npy".format(outdir_tmp,i))

        # array of indices for columns that contain all zeroes
        zero_col_index = np.where(~x.any(axis=0))[0]
        print(zero_col_index)

        # get number of columns present, need for mask length
        num_cols = x.shape[1]
        num_cols_2 = len(k)
        print(num_cols,num_cols_2)

        # create a mask with the column indices
        col_mask = np.zeros(num_cols, dtype=bool)
        col_mask[zero_col_index] = 1
        print(col_mask)

        # apply column mask to matrix and columns

        print("    k before mask ",k.shape)
        k = k[:,col_mask.astype('bool')]
        print("    k after mask  ",k.shape)


        print("    x before mask ",x.shape)
        x = x[:,col_mask.astype('bool')]
        print("    x after mask  ",x.shape)

        np.save("{}matrix{}_filtered.npy".format(outdir_tmp,i), x)
        np.save("{}cols{}_filtered.npy".format(outdir_tmp,i), k)

    del x, k, i
    gc.collect()
    '''

    """
    # --------------------------------------------------------------------------
    # Split
    print("Starting Split")

    nsplits = 5

    # define split
    num_kmers = matrix.shape[1]
    step_size = int(num_kmers / nsplits)
    split_starts = [step_size*i for i in range(nsplits+1)]
    print("num columns (kmers): {}".format(num_kmers))
    print("step size: {}".format(step_size))
    print("split markers: {}".format(split_starts))

    for i in range(nsplits):
        print("  Saving split # {}".format(i),flush=True)

        x = ((matrix.T)[split_starts[i]:split_starts[i+1]]).T
        np.save("{}matrix{}_not_filtered.npy".format(outdir_tmp,i), x)

        k = cols[split_starts[i]:split_starts[i+1]]
        np.save("{}cols{}_not_filtered.npy".format(outdir_tmp,i), cols)

        del x, k
        gc.collect()

    del matrix, cols, i
    gc.collect()

    # --------------------------------------------------------------------------
    # Thresh
    print("Starting Selector (thresh)")

    for i in range(nsplits):
        print("  Selector on split # {}".format(i),flush=True)

        x = np.load("{}matrix{}_not_filtered.npy".format(outdir_tmp,i))
        k = np.load("{}cols{}_not_filtered.npy".format(outdir_tmp,i))

        # define selector
        selector = VarianceThreshold(threshold=0)

        '''
        # array of indices for columns that contain all zeroes
        zero_col_index = np.where(~x.any(axis=0))[0]
        print(zero_col_index)

        num_cols = x.shape[1]
        num_cols_2 = len(k)
        print(num_cols,num_cols_2)

        # create a mask with the column indices
        col_mask = np.zeros(num_cols, dtype=bool)
        col_mask[zero_col_index] = 1
        print(col_mask)

        # apply column mask to matrix and columns
        print("x ",x.shape)
        print("k ",k.shape)

        x = x[:,col_mask.astype('bool')]
        k = k[:,col_mask.astype('bool')]
        '''

        #x_new = np.ma.masked_array(x,col_mask)
        #k_new = np.ma.masked_array(k,col_mask)


        print("x_new ",x_new.shape)
        print("k_new ",k_new.shape)

        sys.exit()


        # perform selector
        x = selector.fit_transform(x)
        k = selector.transform(k.reshape(1,-1)).flatten()

        np.save("{}matrix{}_filtered.npy".format(outdir_tmp,i), x)
        np.save("{}cols{}_filtered.npy".format(outdir_tmp,i), k)

        del x, k
        gc.collect()


    # --------------------------------------------------------------------------
    # Merge
    print("Starting Merge")

    x_new = []
    k_new = []

    # Load a piece of a time and merge it
    for i in range(nsplits):
        print("  Merging split # {}".format(i),flush=True)

        if i == 0:
            x_new = np.load("{}matrix{}_filtered.npy".format(outdir_tmp,i))
        else:
            x_train_tmp = np.load("{}matrix{}_filtered.npy".format(outdir_tmp,i))
            x_new = np.concatenate((x_new,x_train_tmp),axis=1)

        k_tmp = np.load("{}cols{}_filtered.npy".format(outdir_tmp,i))
        for z in k_tmp:
            k_new.append(z)

    np.save("{}matrix_vt.npy".format(out_dir,i), x_new)
    np.save("{}cols_vt.npy".format(out_dir,i), k_new)
    """
