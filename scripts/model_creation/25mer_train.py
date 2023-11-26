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
import gc

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


def print_psutil(msg):
    print("---------------------------------------------------------------")
    print(msg)
    process = psutil.Process(os.getpid())
    print('  RAM memory % used:', psutil.virtual_memory()[2])
    print("  rss",process.memory_info().rss)
    print("  mem",process.memory_info(), flush=True)
    return


if __name__ == "__main__":
    # Get command line arguments
    args = get_args()

    # Required
    model_name = args.modelname.lower()
    matrix = args.matrix
    rows = args.rows
    labels = args.labels
    cols = args.cols
    master_path = args.metadata
    n_feats = args.feats
    pred_for = args.predfor
    outdir = args.outdir
    force_col = args.forcecol
    criteria = args.criteria

    # Optional
    param_configfile = args.paramconfig
    model_type = args.modeltype.lower()
    n_threads = args.threads
    verbose = args.verbose

    # --------------------------------------------------------------------------
    # SO MANY PATHS

    # Folder of all things matrix
    matrix_folder = matrix.split("matrix.npy")[0]

    # Columns that apply to both the training and withheld data
    col_path = cols

    # Training data paths
    #X_path = "{}X.npy".format(matrix_folder)
    X_path = "{}sparse_test/X_sparse.npz".format(matrix_folder)
    Y_path = "{}Y.npy".format(matrix_folder)
    Z_path = "{}Z.npy".format(matrix_folder)

    # Validation data paths
    X_withheld_path = "{}X_withheld.npy".format(matrix_folder)
    Y_withheld_path = "{}Y_withheld.npy".format(matrix_folder)
    Z_withheld_path = "{}Z_withheld.npy".format(matrix_folder)

    # Validation column info
    #X_val_prefilter_path = "{}validation_filtered/X_filtered.npy".format(matrix_folder)
    #cols_val_prefilter_path = "{}validation_filtered/cols_filtered.npy".format(matrix_folder)

    binned_file = "{}binned.pkl".format(matrix_folder)
    nclasses_file = "{}n_classes.pkl".format(matrix_folder)
    class_dict_file = "{}class_dict.pkl".format(matrix_folder)
    label_count_file = "{}label_count.pkl".format(matrix_folder)
    bin_dict_file = "{}bin_dict.pkl".format(matrix_folder)

    str_num_feats = str(n_feats)
    scores_path = "{o}{p}/nested_cv/{f}feats/{m}/".format(
        o=outdir,
        p=pred_for,
        f=str_num_feats,
        m=model_type)
    if not os.path.exists(scores_path):
        os.makedirs(scores_path, exist_ok=True)


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
    print_psutil("Starting training")
    print("  model name {}, model type {}, classifying for {}".format(model_name,model_type,pred_for),flush=True)


    # --------------------------------------------------------------------------
    # Load the master df
    master = pd.read_csv(master_path,sep='\t',low_memory=False)
    # Strip underscores from column names (they confuse snakemake regex)
    master.columns = master.columns.str.replace('_', '')
    # Make all contents of master uppercase for easier comparison
    master = master.apply(lambda x: x.astype(str).str.upper())


    #--------------------------------------------------------------------------
    # Load the config file and base settings

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

    print("Resource settings")
    print("  n_jobs: ", n_jobs)
    print("  n_jobs: ", pre_dispatch)

    # --------------------------------------------------------------------------
    # Grab the data and handle it as needed

    print("Loading & filtering input matrix",flush=True)

    # Grab the data
    #X, Y, Z = get_data(model_name, pred_for, master, matrix, labels, rows, force_col)
    #X, Y, Z = get_data(pred_for, master, matrix, labels, rows, force_col, criteria)
    # Load the matrix, labels, and rows (genomes)
    # 25mer ones were already filtered and saved
    #X = np.load(matrix, allow_pickle=True,mmap_mode='r')
    #Y = np.load(labels, allow_pickle=True)
    #Z = np.load(rows, allow_pickle=True,mmap_mode='r')

    # bin if necessary
    class_count_threshold = 25 #n_folds_inner*n_folds_outer
    #X, Y, Z, binned, n_classes, class_dict, label_count, bin_dict = bin(X,Y,Z,class_count_threshold,pred_for)

    # Withold 20% for validating at the end
    #X, X_withheld, Y, Y_withheld, Z, Z_withheld = train_test_split(X, Y, Z, test_size=0.20, shuffle=True, stratify=Y)
    import scipy.sparse as sp

    print_psutil("Before Load")

    #X = np.load(X_path, allow_pickle=True,mmap_mode='r')
    X = sp.load_npz(X_path)
    Y = np.load(Y_path, allow_pickle=True,mmap_mode='r')
    Z = np.load(Z_path, allow_pickle=True,mmap_mode='r')

    print_psutil("After Load")

    with open(binned_file, 'rb') as fh:
        binned = pickle.load(fh)
    with open(nclasses_file, 'rb') as fh:
        n_classes = pickle.load(fh)
    with open(class_dict_file, 'rb') as fh:
        class_dict = pickle.load(fh)
    with open(label_count_file, 'rb') as fh:
        label_count = pickle.load(fh)
    with open(bin_dict_file, 'rb') as fh:
        bin_dict = pickle.load(fh)

    print("  n_classes: ", n_classes)
    print("  class dict: ", class_dict)
    print("  binned: ", binned)
    print("  label_count: ", label_count)
    print("  bin_dict: ", bin_dict)

    print()

    print("  x shape: ", X.shape)

    #sys.exit()

    '''
    colnames = np.load(cols)
    X = pd.DataFrame(X,columns=colnames)
    print(X)

    import xgboost as xgb
    X = xgb.DMatrix(X, label=Y, feature_names=colnames)
    print(X)
    '''


    #--------------------------------------------------------------------------
    # Set up the param grid

    # Load the parameters
    print("Loading parameters")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # If using best params from 11mers
    '''
    # Load best 11mer model params
    best_11mer_file = "results/wgs_standard/11mer/models/munchkin/best_{}_nested_cv_results.tsv".format(model_type.upper())
    best_11mer_df = pd.read_csv(best_11mer_file,sep='\t')
    # Filter for drug
    param_grid = best_11mer_df.loc[ (best_11mer_df['predfor'] == pred_for), 'params' ].values[0]

    #import ast
    param_grid = ast.literal_eval(param_grid)

    n_feats = param_grid['selector__k']
    '''
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Using default params
    n_feats = int(config['n_feats'])

    param_grid = config[model_type]
    if 'selector__k' in param_grid:
        raise Exception(
            "Don't specify selector__k in the config file, instead define the list n_feats."
            "Each num of feats is run separately and are compared afterwards.")
    if n_feats==0:
        param_grid['selector__k']="all"
    else:
        param_grid['selector__k']=n_feats
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Only defined for ANN
    batch_size = 0
    epochs = 0
    patience = 0

    # If the model is an ANN, then need to expand the params
    if model_type == 'ann':
        epochs = param_grid['epochs']
        batch_size = param_grid['batch_size']
        patience = param_grid['patience']
        # Get the callbacks using patience
        #callbacks = get_keras_callbacks(patience)
        # Get number of hidden layers to search
        n_hidden_layers = param_grid['n_hidden_layers']
        # Make the expanded param grid
        #param_grid = build_keras_grid(n_hidden_layers, param_grid, n_classes, n_feats, patience)

        #abcd = param_grid.pop('epochs', None)
        #efgh = param_grid.pop('batch_size', None)
        #ijkl = param_grid.pop('patience', None)

        param_grid = build_keras_grid(n_hidden_layers, param_grid, n_classes, n_feats, patience)
        param_grid = param_grid[0]
        # take smallest option
        print("*************************")
        print(param_grid)
        #print(param_grid['model_neurons_1'])
        #param_grid['model__neurons_0'] = [param_grid['model__neurons_0'][1]]
        #param_grid['model__neurons_1'] = [param_grid['model__neurons_1'][0]]

        print("*************************")
        print(param_grid)
        print("*************************")


        # Print param grid
        for key,value in param_grid.items():
            print("    ",key,value)
            #param_grid[key] = param_grid[key][0]
            #print("    ",key,value)

        #sys.exit()
    else:
        # Print param grid
        for key,value in param_grid.items():
            print("    ",key,value)


    # --------------------------------------------------------------------------
    # Set up a pipeline of the selector & model
    print("Preparing model")

    # Define the basic model structure
    model = get_base_model(
        model_type,
        n_classes,
        batch_size,
        epochs,
        patience,
        verbose
    )

    # Make a pipeline of the feature selection and the model
    pipe = Pipeline([
        ('selector', SelectKBest(f_classif,k=n_feats)),
        ('model', model)
    ])

    if model_type == 'ann':
        print("ignore the printout of 'param grid'")
        print("the default params in the build_keras_model are actually used")
    else:
        # Set model params
        pipe.set_params(**param_grid)


    # --------------------------------------------------------------------------
    # Single cross validation

    # If defalt scoring and refit is requested
    if scoring_list[0].lower() == 'default':
        scoring_metrics = None
        refit = True
    # If multimetric scoring is requested
    elif scoring_list[0].lower() == 'multi':
        # Define scoring & refit metric(s)
        # Use f1 for binary classification, and f1_micro for multiclass
        if n_classes == 2:
            # Define scoring metrics
            scoring_metrics = {
                'accuracy': 'accuracy',
                'precision': 'precision',
                'recall': 'recall',
                'f1': 'f1'
            }
        else:
            # Define scoring metrics
            scoring_metrics = {
                'accuracy': 'accuracy',
                'precision': 'precision_micro',
                'recall': 'recall_micro',
                'f1': 'f1_micro'
            }
        # Define the metric that will be used to find the best params
        # Note this uses the defined name from the above dicts, so in the case
        # of multiclass, it is using the f1_micro
        #refit_metric = 'f1'
    else:
        raise Exception("Invalid scoring list {}".format(scoring_list))

    print("  scoring metrics: {}".format(scoring_metrics))
    print("  refit metric: {} (taken from above dict)".format(refit_metric))


    print_psutil("BEFORE CROSS VALIDATION")
    #print("--Before cv----------------------------------------------")
    #process = psutil.Process(os.getpid())
    #print('RAM memory % used:', psutil.virtual_memory()[2])
    #print("rss",process.memory_info().rss)
    #print("mem",process.memory_info(), flush=True)


    # Define the isplitting
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True)

    # Give the grid search object to the cross validator; use outer cv.
    # outer_scores is the estimated performance for the final model
    outer_scores = cross_validate(
        pipe,
        X, Y,
        cv=cv,
        n_jobs=n_jobs,
        pre_dispatch=pre_dispatch,
        verbose=verbose,
        scoring=scoring_metrics,
        #scoring='accuracy'
        error_score="raise",
    )


    print("Collected outer scores", flush=True)
    # save outer scores separately just in case
    np.save("{}outer_scores.npy".format(scores_path), outer_scores)
    print("Saved outer scores ************************************", flush=True)


    # **************************************************************************
    # temporary, so can skip redoing the cv step for a typo later on
    #outer_scores_npy = np.load("{}outer_scores.npy".format(scores_path), allow_pickle=True)
    #outer_scores = outer_scores_npy.item()
    #print(outer_scores)


    print_psutil("AFTER CROSS VALIDATION")

    from keras import backend as K
    del X
    K.clear_session()
    gc.collect()


    print_psutil("AFTER CROSS VALIDATION & CLEAR")
    #print("--After cv----------------------------------------------")
    #process = psutil.Process(os.getpid())
    #print('RAM memory % used:', psutil.virtual_memory()[2])
    #print("rss",process.memory_info().rss)
    #print("mem",process.memory_info(), flush=True)



    # --------------------------------------------------------------------------
    # final model making
    print("Preparing final model/pipeline", flush=True)

    #X = np.load(X_val_prefilter_path, allow_pickle=True,mmap_mode='r')
    #X = np.load(X_path, allow_pickle=True,mmap_mode='r')
    X = sp.load_npz(X_path)
    #Z = np.load(Z_path, allow_pickle=True,mmap_mode='r')

    print_psutil("AFTER LOADING THE VALIDATION X")
    #print("--Fit Final Model ----------------------------------------------")
    #process = psutil.Process(os.getpid())
    #print('RAM memory % used:', psutil.virtual_memory()[2])
    #print("rss",process.memory_info().rss)
    #print("mem",process.memory_info(), flush=True)

    # Get best params
    final_params = param_grid #grid.best_params_


    if model_type == 'ann':
        print("ignore the printout of 'param grid'")
        print("the default params in the build_keras_model are actually used")
    else:
        # Set model params
        pipe = pipe.set_params(**final_params)

    # Set up the final model
    #pipe = pipe.set_params(**final_params)
    # Train the final model on the whole dataset
    pipe.fit(X, Y)

    # Save pipe now jusst in case error later on
    save_file('final_model',pipe)

    #*********************************************************************************
    # temp so dont have to redo the fitting
    #mt_lower = model_type.lower()
    #pipe = joblib.load("{}{}_final_model.joblib".format(scores_path,mt_lower))

    print_psutil("AFTER FITTING FINAL MODEL")


    # --------------------------------------------------------------------------
    # validation
    print("Starting validation", flush=True)

    X_withheld = np.load(X_withheld_path, allow_pickle=True,mmap_mode='r')
    X_withheld = sp.csr_matrix(X_withheld)
    Y_withheld = np.load(Y_withheld_path, allow_pickle=True,mmap_mode='r')
    Z_withheld = np.load(Z_withheld_path, allow_pickle=True,mmap_mode='r')

    '''
    # need to make the withheld features match the training set

    # originally 14376033 columns
    original_cols = np.load(cols, allow_pickle=True,mmap_mode='r')
    # trimmed only 25000 columns
    new_cols = np.load(cols_val_prefilter_path, allow_pickle=True,mmap_mode='r')
    # make a column mask
    ori_col_mask = np.in1d(original_cols, new_cols)
    ori_col_mask = np.isin(original_cols, new_cols)
    # apply mask
    X_withheld = X_withheld[:,ori_col_mask]
    '''

    print("*******************************")
    print(X.shape)
    print(X_withheld.shape)
    print("*******************************",flush=True)
    '''
    # save just in case, so dont need to repeat this step later
    np.save("{}X_withheld_25000_filter.npy".format(X_withheld_path),X_withheld)
    '''

    # **************************************************************************************
    #X_withheld = np.load("{}X_withheld_25000_filter.npy".format(X_withheld_path), allow_pickle=True,mmap_mode='r')


    # Now run the withheld data through it to evaluate it
    Y_prediction = pipe.predict(X_withheld)
    # labels and names from the class dict
    #print(class_dict)
    #print()
    labels = list(class_dict.values())
    target_names = list(class_dict.keys())
    # Evaluate
    prfs_report = classification_report(Y_withheld, Y_prediction, labels=labels, target_names=target_names)
    print("\n",prfs_report,"\n",)
    prfs_report_dict = classification_report(Y_withheld, Y_prediction, labels=labels, target_names=target_names, output_dict=True)
    final_results_df = pd.DataFrame(prfs_report_dict)
    prfs_report_df = pd.DataFrame(prfs_report_dict)
    # after scikit v0.20.0, micro avg was replaced with "accuracy"
    # just in case add it as a column
    # it will be the same value for every entry
    blah = precision_recall_fscore_support(Y_withheld, Y_prediction, average='micro')
    final_results_df['micro avg'] = list(blah)
    prfs_report_df['micro avg'] = list(blah)
    print("\n",final_results_df,"\n",)

    # --------------------------------------------------------------------------
    # Look up the values in the results tables
    # the way the tables are laid out, accuracy takes up a whole column with itentical values

    print("--Math Time----------------------------------------------")
    process = psutil.Process(os.getpid())
    print('RAM memory % used:', psutil.virtual_memory()[2])
    print("rss",process.memory_info().rss)
    print("mem",process.memory_info(), flush=True)

    acc_micro  = final_results_df.loc["precision", ["accuracy"]].values[0]
    prec_micro = final_results_df.loc["precision", ["micro avg"]].values[0]
    rec_micro  = final_results_df.loc["recall", ["micro avg"]].values[0]
    f1_micro   = final_results_df.loc["f1-score", ["micro avg"]].values[0]

    acc_macro  = final_results_df.loc["precision", ["accuracy"]].values[0]
    prec_macro = final_results_df.loc["precision", ["macro avg"]].values[0]
    rec_macro  = final_results_df.loc["recall", ["macro avg"]].values[0]
    f1_macro   = final_results_df.loc["f1-score", ["macro avg"]].values[0]

    # Put the results into a dictionary & turn it into a pandas df
    final_results_micro = {
        'predfor':[pred_for],
        'model':[model_type],
        'feats': [n_feats],
        'acc':[acc_micro],
        'precision':[prec_micro],
        'recall':[rec_micro],
        'f1':[f1_micro],
        'params':[final_params],
        'binned': [binned],
    }
    final_results_micro = {**final_results_micro,**label_count}
    final_results_micro_df = pd.DataFrame.from_dict(final_results_micro)

    #print()
    #print(label_count)
    #print()
    #print(final_results_micro)
    #print()
    #print(final_results_micro_df)
    #print()

    #sys.exit()

    final_results_macro = {
        'predfor':[pred_for],
        'model':[model_type],
        'feats': [n_feats],
        'acc':[acc_macro],
        'precision':[prec_macro],
        'recall':[rec_macro],
        'f1':[f1_macro],
        'params':[final_params],
        'binned': [binned],
    }
    final_results_macro = {**final_results_macro,**label_count}
    final_results_macro_df = pd.DataFrame.from_dict(final_results_macro)

    # --------------------------------------------------------------------------
    # Get the final cv performance estimates & results of param tuning

    print("Formatting Results", flush=True)

    # Estimated performance of the final model
    # Calculate the mean, standard deviation, and variance for:
    #   accuracy, precision, f1 score, fit time
    def mean_std_var(value_list):
        """
        Calculate and return the mean, standard deviation, and variance of
        the input list
        """
        m = np.mean(value_list)
        s = np.std(value_list)
        v = np.var(value_list)
        return m,s,v

    #***************************************************************************
    print(outer_scores)
    #***************************************************************************
    #***************************************************************************
    if scoring_list[0].lower() == 'default':
        # Accuracy
        accuracy_list = outer_scores['test_score']
        accuracy_mean,accuracy_std,accuracy_var = mean_std_var(accuracy_list)
        # Precision
        precision_list = np.nan
        precision_mean = np.nan
        precision_std = np.nan
        precision_var = np.nan
        # Recall
        recall_list = np.nan
        recall_mean = np.nan
        recall_std = np.nan
        recall_var = np.nan
        # F1 Score
        f1_list = np.nan
        f1_mean = np.nan
        f1_std = np.nan
        f1_var = np.nan
        # Fit Time
        fit_time_list = outer_scores['fit_time']
        fit_time_mean,fit_time_std,fit_time_var = mean_std_var(fit_time_list)
        # Score Time
        score_time_list = outer_scores['score_time']
        score_time_mean,score_time_std,score_time_var = mean_std_var(score_time_list)
    # If multimetric scoring is requested
    elif scoring_list[0].lower() == 'multi':
        # Accuracy
        accuracy_list = outer_scores['test_accuracy']
        accuracy_mean,accuracy_std,accuracy_var = mean_std_var(accuracy_list)
        # Precision
        precision_list = outer_scores['test_precision']
        precision_mean,precision_std,precision_var = mean_std_var(precision_list)
        # Recall
        recall_list = outer_scores['test_recall']
        recall_mean,recall_std,recall_var = mean_std_var(recall_list)
        # F1 Score
        f1_list = outer_scores['test_f1']
        f1_mean,f1_std,f1_var = mean_std_var(f1_list)
        # Fit Time
        fit_time_list = outer_scores['fit_time']
        fit_time_mean,fit_time_std,fit_time_var = mean_std_var(fit_time_list)
        # Score Time
        score_time_list = outer_scores['score_time']
        score_time_mean,score_time_std,score_time_var = mean_std_var(score_time_list)
    #***************************************************************************

    # Get the final parameters
    #final_params = grid.best_params_
    #print(grid.best_estimator_)
    #print(grid.best_score_)
    #print(grid.best_params_)

    # Inner loop detailed results, will save just in case
    #inner_results_dict = grid.cv_results_
    #inner_results_df = pd.DataFrame.from_dict(inner_results_dict)

    # Put the results into a dictionary & turn it into a pandas df
    outer_results = {
        'predfor':[pred_for],
        'model':[model_type],
        'feats': [n_feats],
        'acc_mean':[accuracy_mean],
        'acc_std':[accuracy_std],
        'acc_var':[accuracy_var],
        'precision_mean':[precision_mean],
        'precision_std':[precision_std],
        'precision_var':[precision_var],
        'recall_mean':[recall_mean],
        'recall_std':[recall_std],
        'recall_var':[recall_var],
        'f1_mean':[f1_mean],
        'f1_std':[f1_std],
        'f1_var':[f1_var],
        'params':[final_params],
        'binned': [binned],
        'acc_all': [accuracy_list],
        'prec_all': [precision_list],
        'f1_all': [f1_list],
        'fit_time_mean':[fit_time_mean],
        'fit_time_std':[fit_time_std],
        'fit_time_var':[fit_time_var],
        'score_time_mean':[score_time_mean],
        'score_time_std':[score_time_std],
        'score_time_var':[score_time_var],
        'fit_time_all': [fit_time_list],
        'score_time_all': [score_time_list]
    }
    outer_results = {**outer_results,**label_count}
    outer_results_df = pd.DataFrame.from_dict(outer_results)


    # --------------------------------------------------------------------------
    # Get the top feats

    print("Extracting feats. Currently, importances only saved if model is XGB.")

    # Get the top N feats that were used by the selector
    # Get the bool list of whether or not a feature was used
    mask = pipe['selector'].get_support()
    # Load matrix column headers (kmer sequences)
    #kmer_cols = np.load(cols)
    #kmer_cols = np.load(cols_val_prefilter_path, allow_pickle=True,mmap_mode='r')
    kmer_cols = np.load(cols, allow_pickle=True)
    # Apply the mask to the list of kmers to obtain the top n_feats
    final_feats = kmer_cols[mask]

    # Get the feat importances
    if model_type.lower() == 'xgb':
        print("Extracting feature importance.")
        # currently cant actually handle gblinear
        if pipe['model'].booster == 'gblinear':
            # INCOMPLETE / NOT WORKING
            model_coef = pipe['model'].coef_
            feat_imptces = dict(zip(final_feats, model_coef))
            feat_imptces = {k: abs(v / sum(model_coef)) for k, v in feat_imptces.items()}
            feat_imptces = [ abs(v / sum(model_coef)) for v in model_coef ]
            feat_imptces = np.array(feat_imptces)
            #print('gblinear')
            #print(feat_imptces)
            raise Exception("Currently can only handle gbtree")
        else:
            #print('gbtree')
            feat_imptces = pipe['model'].feature_importances_
            #print(feat_imptces)


    # --------------------------------------------------------------------------
    # Save Everything

    print("Saving Results")

    # Save the class dict
    save_file('class_dict',class_dict)
    # Save the clade binning if needed
    if 'clade' in pred_for.lower() and len(label_count)>2:
        save_file('bin_dict',bin_dict)
    # Save the final model
    save_file('final_model',pipe)
    # Save the final selectkbest feats
    save_file('selectkbest_feats',final_feats)
    # Save the genome names of the withheld Z
    save_file('Z_withheld',Z_withheld)
    save_file('Z_train',Z)
    # Save the final results
    save_file('prfs_report',prfs_report_df)
    save_file('micro_results',final_results_micro_df)
    save_file('macro_results',final_results_macro_df)
    # Save the outer results
    save_file('outer_results',outer_results_df)
    # Save the inner results
    #save_file('inner_results',inner_results_df)

    if model_type.lower() == 'xgb':
        final_feats_vstack = np.vstack((final_feats.flatten(), feat_imptces))
        save_file('feature_importance',final_feats_vstack)



    ############################################################################
    # Check that model can be loaded
    '''
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("LOADING MODEL")
    # Get & parse args
    args = get_args()
    outdir = args.outdir
    pred_for = args.predfor
    str_num_feats = str(args.feats)
    model_type = args.modeltype.lower()
    # Want to save in labelled sub folders
    out_path = "{o}{p}/nested_cv/{f}feats/{m}/".format(
        o=outdir,
        p=pred_for,
        f=str_num_feats,
        m=model_type)
    final_model_file = "{}{}_final_model.joblib".format(out_path,model_type)
    pipe = load(final_model_file)
    print("LOADED MODEL SUCCESSFULLY")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    '''
    ############################################################################



    '''
    # Number of feats as a string so it can be used in file paths
    #str_n_feats = str(len(final_feats))
    str_n_feats = str(n_feats)

    # Want to save in labelled sub folders
    out_path = "{o}{p}/nested_cv/{f}feats/{m}/".format(
        o=outdir,
        p=pred_for,
        f=str_n_feats,
        m=model_type)
    if not os.path.exists(out_path):
        os.makedirs(out_path, exist_ok=True)

    # 1. Save the class dictionary
    class_dict_file = "{}class_dict.pkl".format(out_path)
    with open(class_dict_file, 'wb') as fh:
        pickle.dump(class_dict, fh)

    #1b. if the clades were binned, save the bin dict
    if 'clade' in pred_for.lower() and len(label_count)>2:
        bin_dict_file = "{}bin_dict.pkl".format(out_path)
        with open(bin_dict_file, 'wb') as fh:
            pickle.dump(bin_dict, fh)

    # 2. Save the inner detailed results in case we want to look at them later
    inner_results_file = "{}inner_results_detail.tsv".format(out_path)
    inner_results_df.to_csv(inner_results_file, sep='\t', index=False)

    # 3. Save the final model's performance & specs
    final_results_file = "{}final_model_results.tsv".format(out_path)
    final_results_df.to_csv(final_results_file, sep='\t', index=False)

    # 4. Save the final model; joblib is recommended by sklearn
    final_model_file = "{}{}_final_model.joblib".format(out_path,model_type)
    dump(pipe,final_model_file)

    # 5. Save the top n feats
    final_feats_file = "{}/selectkbest_{}_feats.pkl".format(out_path,str_n_feats)
    with open(final_feats_file, 'wb') as fh:
        pickle.dump(final_feats, fh)

    # 6. If the model is xgb, also save the importance ranking
    if model_type.lower() == 'xgb':
        impt_feat_file = "{}/feature_importance.npy".format(out_path,str_n_feats)
        np.save(impt_feat_file, np.vstack((final_feats.flatten(), feat_imptces)))
    '''
