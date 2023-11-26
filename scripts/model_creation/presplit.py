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
import psutil

# Helper scripts for grabbing data etc
from model_helper import *

# ML models
from sklearn.svm import SVC # svm
import xgboost as xgb # xgb
from keras.wrappers.scikit_learn import KerasClassifier # ann

# sklearn basics
import sklearn as skl
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import cross_val_score, cross_validate
#from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.base import clone
from sklearn.pipeline import Pipeline

# sklearn scoring metrics
from sklearn import metrics
from sklearn.metrics import precision_score, accuracy_score, f1_score, recall_score, classification_report, precision_recall_fscore_support




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
    #force_col = args.forcecol
    #criteria = args.criteria

    # Optional
    #param_configfile = args.paramconfig
    #model_type = args.modeltype.lower()
    #n_threads = args.threads
    #verbose = args.verbose

    # --------------------------------------------------------------------------
    # Suppress certain warnings
    # Warnings that are meaningless for me and just make it harder to see
    # real errors.
    # Can put them back in after done testing

    # Tensorflow's excess information
    # Every single fold spams the
    #   'Your CPU supports instructions that this TensorFlow binary was not
    #   compiled to use: ...'

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    # 0 = all messages are logged (default behavior)
    # 1 = INFO messages are not printed
    # 2 = INFO and WARNING messages are not printed
    # 3 = INFO, WARNING, and ERROR messages are not printed


    # --------------------------------------------------------------------------
    # Load the master df
    master = pd.read_csv(master_path,sep='\t',low_memory=False)
    # Strip underscores from column names (they confuse snakemake regex)
    master.columns = master.columns.str.replace('_', '')
    # Make all contents of master uppercase for easier comparison
    master = master.apply(lambda x: x.astype(str).str.upper())


    #--------------------------------------------------------------------------
    # Load the config file and base settings
    '''
    yaml = YAML(typ='safe')
    with open(param_configfile, 'r') as f:
        config = yaml.load(f)

    # Get the basic settings
    n_folds = config['n_folds']
    # scoring
    scoring_list = config['scoring_settings'][model_type]['scoring']
    refit_metric = config['scoring_settings'][model_type]['refit']
    # resources
    n_jobs = config['resource_settings'][model_type]['n_jobs']
    pre_dispatch = config['resource_settings'][model_type]['pre_dispatch']
    use_dask = config['resource_settings'][model_type]['dask']
    '''

    # --------------------------------------------------------------------------
    # Grab the data and handle it as needed

    print("Loading & filtering input matrix")

    # Grab the data
    #X, Y, Z = get_data(model_name, pred_for, master, matrix, labels, rows, force_col)
    #X, Y, Z = get_data(pred_for, master, matrix, labels, rows, force_col, criteria)
    # Load the matrix, labels, and rows (genomes)
    # 25mer ones were already filtered and saved
    X = np.load(matrix, allow_pickle=True,mmap_mode='r')
    Y = np.load(labels, allow_pickle=True)
    Z = np.load(rows, allow_pickle=True,mmap_mode='r')

    # bin if necessary
    class_count_threshold = 25 #n_folds_inner*n_folds_outer
    X, Y, Z, binned, n_classes, class_dict, label_count, bin_dict = bin(X,Y,Z,class_count_threshold,pred_for)

    # Withold 20% for validating at the end
    X, X_withheld, Y, Y_withheld, Z, Z_withheld = train_test_split(X, Y, Z, test_size=0.20, shuffle=True, stratify=Y)


    print("  n_classes: ",n_classes)
    print("  class dict: ",class_dict)


    if not os.path.exists(outdir):
        os.makedirs(out_path, exist_ok=True)

    np.save("{}X.npy".format(outdir), X)
    np.save("{}X_withheld.npy".format(outdir), X_withheld)

    np.save("{}Z.npy".format(outdir), Z)
    np.save("{}Z_withheld.npy".format(outdir), Z_withheld)

    np.save("{}Y.npy".format(outdir), Y)
    np.save("{}Y_withheld.npy".format(outdir), Y_withheld)


    binned_file = "{}binned.pkl".format(outdir)
    with open(binned_file, 'wb') as fh:
        pickle.dump(binned, fh)

    nclasses_file = "{}n_classes.pkl".format(outdir)
    with open(nclasses_file, 'wb') as fh:
        pickle.dump(n_classes, fh)

    class_dict_file = "{}class_dict.pkl".format(outdir)
    with open(class_dict_file, 'wb') as fh:
        pickle.dump(class_dict, fh)

    label_count_file = "{}label_count.pkl".format(outdir)
    with open(label_count_file, 'wb') as fh:
        pickle.dump(label_count, fh)

    bin_dict_file = "{}bin_dict.pkl".format(outdir)
    with open(bin_dict_file, 'wb') as fh:
        pickle.dump(bin_dict, fh)
