
import pandas as pd
import numpy as np
import pickle
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ruamel.yaml import YAML
from joblib import dump, load
from model_helper import *




if __name__== "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--id', help='id', required=True)
    parser.add_argument('--infile', help='infile', required=True)
    parser.add_argument('--outfile', help='outfile', required=True)
    parser.add_argument('--modeltype', help='model type', required=True)
    parser.add_argument('-k', help='kmer len', required=True)
    parser.add_argument('-g', help='model name', required=True)
    args = parser.parse_args()


    id = args.id
    in_file = args.infile
    out_file = args.outfile
    model_type = args.modeltype
    k = args.k
    g = args.g
    w = 'standard'


    drugs = [ 'AMC', 'AMP', 'AMX', 'AZT',
             'CET', 'CFZ', 'CHL', 'CIP', 'CPM', 'CRO', 'CST', 'CTX', 'CTZ', 'CXA', 'CXM',
             'EFX', 'ETP', 'FFC', 'FOX', 'GEN', 'IMP', 'KAN', 'LVX', 'MER', 'NAL', 'NIT',
             'SAM', 'SOX', 'STR', 'SXT',
             'TBM', 'TET', 'TGC', 'TIO', 'TMP', 'TZP' ]

    #best = "results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=model_type)
    best = "results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=model_type.upper())
    df = pd.read_csv(best,sep='\t')


    # init the dict that will contain drug: S/I/R
    pred_dict ={'id':id}
    for drug in drugs:
        print("------------------------------")
        print(drug)
        # get the path to the best xgb model for current drug
        #g= 'munchkin'#'all'
        #w = wgs_type
        #k = str(config['kmer_size'])
        #best = "/Drives/K/janmoat/mecoli/results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=model_type)
        #best = "results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=model_type)

        #df = pd.read_csv(best,sep='\t')
        # find the path to the best model for current pred_for
        path = df.loc[df['predfor'] == 'SIR'+drug, 'path'].values[0]
        mpath = path + model_type + "/"
        # get model and class dict for this model
        model_file = mpath + model_type + "_final_model.joblib"
        class_dict_file = mpath + "class_dict.pkl"
        # Load the class dict to decode the prediction
        with open(class_dict_file, 'rb') as fh:
            class_dict = pickle.load(fh)

        print(class_dict)

        # Load the input sample
        data = np.load(in_file)
        # "Reshape your data either using array.reshape(-1, 1) if your data has
        #a single feature or array.reshape(1, -1) if it contains a single sample."
        data = data.reshape(1,-1)
        # Load the model pipeline
        #from model_creation.model_helper import *
        pipe = load(model_file)
        # Apply the model pipeline to the sample
        predn = pipe.predict(data)
        if isinstance(predn[0],(list,pd.core.series.Series,np.ndarray)):
            predn = predn[0]
            print(predn)
        print(predn)
        # Decode the prediction
        predn = decode_labels(predn, class_dict)

        # Add the prediction to the dict
        pred_dict[drug+'_'+model_type] = predn
    # convert to dataframe
    df = pd.DataFrame.from_dict(pred_dict)
    # save
    df.to_csv(out_file, index=False, sep='\t')
