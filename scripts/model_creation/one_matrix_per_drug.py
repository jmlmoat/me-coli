#!/usr/bin/env python

# Basics
import pandas as pd
import numpy as np
import argparse
import pickle
import os, sys
import joblib
from joblib import dump, load
from ruamel.yaml import YAML

# Helper scripts for grabbing data etc
from model_helper import *

# Filter out the constant columns (eg all 0s etc)
from sklearn.feature_selection import VarianceThreshold

if __name__ == "__main__":
    # Get command line arguments
    #args = get_args()

    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        '-m',
        '--matrix',
        help='path to k-mer frequency matrix',
        required=True)
    parser.add_argument(
        '-r',
        '--rows',
        help='path to k-mer frequency matrix rows (genome names)',
        required=True)
    parser.add_argument(
        '-l',
        '--labels',
        help='path to k-mer frequency matrix labels',
        required=True)
    parser.add_argument(
        '-c',
        '--cols',
        help='path to k-mer frequency matrix columns (k-mer sequences)',
        required=True)
    parser.add_argument(
        '-mdf',
        '--metadata',
        help='path to metadata master file',
        required=True)
    parser.add_argument(
        '-p',
        '--predfor',
        help='thing to predict for ie SIRAMP',
        required=True)
    parser.add_argument(
        '-o',
        '--outdir',
        help='where to save output',
        required=True)
    parser.add_argument(
        '-fc',
        '--forcecol',
        help='name of column to filter the data for; case must match metadata',
        required=True)
    parser.add_argument(
        '-cr',
        '--criteria',
        help='list of criteria to filter for in the given column',
        nargs='+',
        required=True)

    args = parser.parse_args()

    # Required
    #model_name = args.modelname.lower()
    matrix = args.matrix
    rows = args.rows
    labels = args.labels
    cols = args.cols
    master_path = args.metadata
    #n_feats = args.feats
    pred_for = args.predfor
    outdir = args.outdir
    force_col = args.forcecol
    criteria = args.criteria

    # --------------------------------------------------------------------------
    # Load the master df
    # Load only the relevant columns of the master df
    if pred_for == 'clade':
        relevant_columns = [ 'id', 'collection', 'dataset', pred_for]
    else:
        predfor_col = pred_for.replace('SIR','SIR_').upper()
        relevant_columns = [ 'id', 'collection', 'dataset', predfor_col]
    # Load the master df
    master = pd.read_csv(master_path, sep='\t', usecols=relevant_columns, low_memory=False)
    # Strip underscores from column names (they confuse snakemake regex)
    master.columns = master.columns.str.replace('_', '')
    # Make all contents of master uppercase for easier comparison
    master = master.apply(lambda x: x.astype(str).str.upper())

    # --------------------------------------------------------------------------
    # Grab the data and handle it as needed

    # Load the full matrix data, remove info not relevant for the pred_for

    X, Y, Z = get_data(pred_for, master, matrix, labels, rows, force_col, criteria)

    cols = np.load(cols, allow_pickle=True)

    # --------------------------------------------------------------------------
    # Save
    if not os.path.exists(outdir):
        os.makedirs(out_path, exist_ok=True)
    np.save("{}matrix.npy".format(outdir), X)
    np.save("{}rows.npy".format(outdir), Z)
    np.save("{}labels.npy".format(outdir), Y)
    np.save("{}cols.npy".format(outdir), cols)
