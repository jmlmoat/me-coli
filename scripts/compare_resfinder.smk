

import pandas as pd
import numpy as np
import pickle
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ruamel.yaml import YAML
from joblib import dump, load
from model_creation.model_helper import *

configfile: "config/config.yaml"

# Settings
wgs_type = config['wgs_type']
#ksnp_min_frac = config['ksnp_min_frac'] # smaller = faster
kmer_len = config['kmer_size']
matrix_dtype = config['kmer_matrix_dtype']

# Drugs we have models for
drugs = config['drugs_we_care_about'] # drugs to predict on
#drugs = [ AMC, AMP, AMX, AZT,
#         CET, CFZ, CHL, CIP, CPM, CRO, CST, CTX, CTZ, CXA, CXM,
#         EFX, ETP, FFC, FOX, GEN, IMP, KAN, LVX, MER, NAL, NIT,
#         SAM, SOX, STR, SXT,
#         TBM, TET, TGC, TIO, TMP, TZP ]

# ******************************************************************************
# ResFinder 4.0 paper's 390 genome val set (DE)
'''
predgroup = 'resf_val_data_DE'
resf_prj = "metadata/public/PRJNA6164_metadata.tsv"
resf_df = pd.read_csv(resf_prj, sep='\t')
ids = resf_df['id'].tolist()
#resf_val_drugs = ['AMP','CHL','CIP','CTX','CTZ','GEN','NAL','TET']
# Full list of drugs available for the resf DE dataset = ['AMP','CHL','CIP','CPM','CTX','CTZ','ETP','GEN','IMP','MER','NAL','TET','TGC']
resf_val_drugs = ['AMP','CHL','CIP','CPM','CTX','CTZ','ETP','GEN','IMP','MER','NAL','TET','TGC']
'''
#******************************************************************************
# ResFinder 4.0 paper's DK genome val set

predgroup = 'resf_val_data_DK'
resf_prj = "metadata/public/PRJEB22091_metadata.tsv"
resf_df = pd.read_csv(resf_prj, sep='\t')
ids = resf_df['id'].tolist()
#resf_val_drugs = ['AMP','CHL','CIP','CTX','CTZ','FOX','GEN','NAL','TET','TMP']
# Full list of drugs available for the resf DE dataset = ['AMP','CHL','CIP','CPM','CST','CTX','CTZ','ETP','FOX','GEN','IMP','MER','NAL','SMZ','TET','TMP']
resf_val_drugs = ['AMP','CHL','CIP','CPM','CST','CTX','CTZ','ETP','FOX','GEN','IMP','MER','NAL','TET','TMP']


#******************************************************************************


# Paths
JELLY = "results/wgs_{}/{}mer/jellyfish/".format(wgs_type,kmer_len)
PRED = "results/wgs_{}/{}/".format(wgs_type,predgroup)
EXT = "results/wgs_{}/external_amr_tools/".format(wgs_type)
RESF_OUTPUT = "results/wgs_{}/external_amr_tools/resfinder/".format(wgs_type)
RAW   = "data/genomes/raw/" # raw genomes
#raw_test = "results/wgs_{}/{}/raw_test/".format(wgs_type,predgroup)

# Path to resfinder location; will be downloaded if not present
resfinder_path = config['resfinder_path']
if '~' in resfinder_path:
    resfinder_path = os.path.expanduser(resfinder_path)

# Download genomes
subworkflow download:
    workdir: "../"
    snakefile: "download.smk"

# Assemble WGS data
subworkflow assemble:
    workdir: "../"
    snakefile: "assemble.smk"

# k-mer counting
subworkflow jellyfish:
    workdir: "../"
    snakefile: "jellyfish.smk"

# Comparison with other amr databases
subworkflow external_amr_tools:
    workdir: "../"
    snakefile: "external_amr_tools.smk"

rule all:
    input:
    ## Download & Process Data
    # Download the resfinder paper genomes (390 from ENA)
        #download(expand("data/genomes/raw/{id}/{id}_{pe}.fastq",
        #    id=resf_ids,
        #    pe=[1,2])),
    # Assemble (will cause download if not already there)
        #assemble(expand("data/genomes/asmbl/{id}/scaffolds.fasta",
        #    id=ids))
    # Count k-mers for all wgs data
        #jellyfish(expand(JELLY+"{id}.fa", id=ids)),
    # Format input for the models
        #expand("{f}model_input/{id}.npy",f=PRED,id=ids),
    # --------------------------------------------------------------------------
    ## Predict with ML Models
        #expand("{f}model_output/{id}/{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm']),
    # put every sample into one tsv
        #expand("{f}prediction_{mt}.tsv",f=PRED,mt=['xgb','svm'])
    # Model accuracy for each model type
        #expand("{f}ecoff_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm']),
        #expand("{f}clsi_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm']),
        #expand("{f}eucast_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm'])
        #************************************************************************************
        expand("{f}summary_model/{standard}_{mt}.tsv",
            f=PRED,
            mt=['xgb','svm', 'ann'],
            standard=['clsi','eucast','ecoff']),
    # --------------------------------------------------------------------------
    ## Predict with Resfinder 4.0 (not abricate)
    # extracts a minimalist df from the summary_by_sample output
        #expand("{f}resf_output/{id}.tsv",f=PRED,id=ids),
    # Get accuracy for each standard and each of raw/assembled
        #************************************************************************************
        expand(PRED+"summary_resfinder/{db}_acc_{standard}.tsv",
            db=['rawcgeresfinder','asmblcgeresfinder'],
            standard=['clsi','eucast','ecoff']),
    # --------------------------------------------------------------------------
    ## Comparison
        #expand("{f}model_output/{id}/{mt}.tsv",f=PRED,id=ids,mt=['xgb']),
        #expand("{f}resf_output/{id}.tsv",f=PRED,id=ids),
        #expand("{f}ecoff_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm','resfinder']),
        #expand("{f}clsi_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm','resfinder']),
        #expand("{f}eucast_{mt}.tsv",f=PRED,id=ids,mt=['xgb','svm','resfinder'])

################################################################################
# ResFinder Prediction & Accuracy
################################################################################

def get_all_cge_resistance_lists(wcs):
    if wcs.db == 'asmblcgeresfinder':
        return external_amr_tools(expand(RESF_OUTPUT+"asmbl/{id}/cgeresfinder_resistance.pkl",id=ids))
    elif wcs.db == 'rawcgeresfinder':
        return external_amr_tools(expand(RESF_OUTPUT+"raw/{id}/cgeresfinder_resistance.pkl",id=ids))
    else:
        raise Exception("Wrong db")

# NOT the same as the one in external_amr_tools, this is adjusted for the
# different bioproject metadata! different col names too, and drug list
rule comp_res_fast_table:
    group: 'comp_res_fast_table'
    input:
        get_all_cge_resistance_lists
    output:
        out = PRED+"summary_resfinder/{db}_predictions_{standard}.tsv"
    run:
        # Database name
        db = wildcards.db
        #standard comparing to
        standard = wildcards.standard

        # Load the metadata master
        master_df = pd.read_csv(resf_prj,sep='\t')
        print(master_df)

        # Get list of only drugs we made models for
        drugs = resf_val_drugs
        # Load drug information; only need certain columns
        drug_df = pd.read_csv(config['antibiotic_master'],sep='\t',
            usecols=['code', 'name', 'name_alt', 'group_broad', 'group',
                     'class', 'subclass', 'card_class', 'combination'])
        # Filter for the drugs in the list
        drug_df = drug_df[drug_df['code'].isin(drugs)]
        drug_df = drug_df.reset_index(drop=True)
        # Make names lowercase
        drug_df['name'] = drug_df['name'].astype(str).apply(lambda x: x.lower())
        drug_df['name_alt'] = drug_df['name_alt'].astype(str).apply(lambda x: x.lower())
        # Make the drugs list match the order of the df
        #  - dont remember why i did this
        drugs = drug_df['code'].tolist()

        # Initialize a list to store dicts; will be used to make final dataframe
        row_list = []

        # For each file (resistance list) in the input
        for file in input:
            # Gt the id from the filename
            id = file.split('/')[-2]

            # Load the list of resistances
            with open(file, 'rb') as r:
                res_list = pkl.load(r)

            # Initialize a dict for the current isolate
            new_row = { 'id':id }

            # Add the SIR prediction and the SIR lab value to the dict
            for drug in drugs:
                # Look up the lab SIR value
                metadata_row = master_df.loc[master_df['id'] == id]
                amr_lab = metadata_row[drug+"_"+standard].values[0]
                # Append the lab value to the dict; eg AMP_lab
                new_row[drug+'_lab'] = amr_lab

                # Look up the SIR prediction in the resistance list

                # Get the drug info row of the df for current drug
                drug_row = drug_df.loc[drug_df['code'] == drug]

                # Create a list of the drug's names, ensure no nans
                ori_drug_name = drug_row['name'].values.tolist()[0]
                drug_names = drug_row['name'].values.tolist()
                # the alts cell has items delimited by commas but is not
                # actually in list form yet
                alts = drug_row['name_alt'].values.tolist()[0].split(',')
                # strip whitepaces
                alts = [a.strip(' ') for a in alts]
                # make one big list of drug names
                drug_names.extend(alts)
                # no nans!
                drug_names = list(filter(lambda a: a != 'nan', drug_names))
                if np.nan in drug_names or 'nan' in drug_names or 'NaN' in drug_names:
                    raise Exception("nans!")

                # Get other info for the drug
                group_broad = drug_row['group_broad'].values[0]
                group = drug_row['group'].values[0]
                drug_class = drug_row['class'].values[0]
                subclass = drug_row['subclass'].values[0]
                card_class = drug_row['card_class'].values[0]

                # Drug class may be a list of items ie 'sulfonamide, sulphonamide'
                # Convert it into a list
                if ',' in drug_class:
                    drug_class = drug_class.split(',')
                    drug_class = [i.strip(' ') for i in drug_class]
                else:
                    drug_class = [drug_class]

                # If no resistance is found, fill with nan
                amr_pred = np.nan
                db_class = np.nan

                # TGC is a special case for resfinder
                if drug == 'TGC':
                    print(drug_class, db_class)
                    if db in ['asmblcgeresfinder','rawcgeresfinder','resfinder']:
                        if 'tigecycline' in res_list:
                            amr_pred = 'R'
                            db_class = 'tigecycline'

                # SXT is a special case
                elif drug == 'SXT':
                    # The dbs have no single entry for SXT (trimethoprim with
                    # sulfamethoxazole) but detect both components individually.
                    # cge
                    #   cgeresfinder: sulphonamide, trimethoprim
                    #   resfinder: sulfamethoxazole, trimethoprim
                    # ncbi
                    #   ncbi: sulfonamide, trimethoprim
                    #   amf: sulfonamide, trimethoprim
                    # card
                    #   card: sulfonamide, diaminopyrimidine
                    #   rgi: sulfonamide, diaminopyrimidine
                    if db in ['asmblcgeresfinder','rawcgeresfinder']:
                        if 'sulfamethoxazole' in res_list \
                        and 'trimethoprim' in res_list:
                            amr_pred = 'R'
                            db_class = 'sulfamethoxazole, trimethoprim'
                    elif db == 'resfinder':
                        if 'sulfamethoxazole' in res_list \
                        and 'trimethoprim' in res_list:
                            amr_pred = 'R'
                            db_class = 'sulfamethoxazole, trimethoprim'
                    elif db == 'ncbi' or db == 'amf':
                        if 'sulfonamide' in res_list \
                        and 'trimethoprim' in res_list:
                            amr_pred = 'R'
                            db_class = 'sulfonamide, trimethoprim'
                    elif db == 'card' or db == 'rgi':
                        if 'sulfonamide' in res_list \
                        and 'diaminopyrimidine' in res_list:
                            amr_pred = 'R'
                            db_class = 'sulfonamide, diaminopyrimidine'

                # Special cases where we only want to look at cephalosporin
                # but not beta-lactam (large accuracy boost)
                elif db in ['amf','ncbi'] \
                and drug in ['CET','CRO','CTX','CTZ','CXM','TIO']:
                    # check names first
                    if any(d in drug_names for d in res_list):
                        amr_pred = 'R'
                        # same drug name as it appears in res_list
                        #db_class += [i for i in drug_names if i in res_list]
                        # what names appeared in both lists
                        overlap = list(set(drug_names).intersection(res_list))
                        db_class = ','.join(overlap)
                    # if the drug class is in the res list
                    # class is a list for [sulfonamide, sulphonamide] because
                    # multiple spellings
                    #elif drug_class in res_list:
                    elif any(d in drug_class for d in res_list):
                        amr_pred = 'R'
                        #db_class += [drug_class]
                        overlap = list(set(drug_class).intersection(res_list))
                        db_class = ','.join(overlap)

                # rest of the drugs proceed the same
                else:
                    # if one of the drug names is in the resistant list
                    if any(d in drug_names for d in res_list):
                        amr_pred = 'R'
                        # same drug name as it appears in res_list
                        #db_class += [i for i in drug_names if i in res_list]
                        # what names appeared in both lists
                        overlap = list(set(drug_names).intersection(res_list))
                        db_class = ','.join(overlap)
                    # if the drug class is in the res list
                    # class is a list for [sulfonamide, sulphonamide] because
                    # multiple spellings
                    #elif drug_class in res_list:
                    elif any(d in drug_class for d in res_list):
                        amr_pred = 'R'
                        #db_class += [drug_class]
                        overlap = list(set(drug_class).intersection(res_list))
                        db_class = ','.join(overlap)
                    # if the card class is in the res list
                    # CARD considers nalidixic acid fluoroquinolone while others dont
                    elif (db in ['card','rgi']) and card_class in res_list:
                        amr_pred = 'R'
                        db_class = card_class
                    # group eg cephem
                    elif group in res_list:
                        amr_pred = 'R'
                        db_class = group
                    # group_broad group eg beta-lactam
                    elif group_broad in res_list:
                        amr_pred = 'R'
                        db_class = group_broad

                # Append the prediction to the dict; eg AMP_pred
                new_row['res_list'] = ','.join(res_list)
                new_row[drug+'_pred'] = amr_pred
                new_row[drug+'_db_class'] = db_class

            # append the new row to the row list
            row_list.append(new_row)

        # Turn the list of rows into a dataframe
        out_df = pd.DataFrame(row_list)

        #pd.set_option('display.max_columns', 500)
        #pd.set_option('display.width', 1000)
        #print(out_df)
        #sys.exit()

        # move the column containing the resistance list to the second column,
        # next to the ids, and before any SIR columns
        col_list = list(out_df)
        col_list.insert(1, col_list.pop(col_list.index('res_list')))
        out_df = out_df.reindex(columns=col_list)
        print(out_df)

        # Save dataframe
        out_df.to_csv(output.out, index=False, sep='\t')

rule comp_res_fast_table_compare:
    group: 'comp_res_fast_table_compare'
    input:
        pred = PRED+"summary_resfinder/{db}_predictions_{standard}.tsv"
    output:
        out = PRED+"summary_resfinder/{db}_acc_{standard}.tsv"
    run:
        # Load the file
        all_df = pd.read_csv(input.pred,sep='\t',low_memory=False)
        #print(all_df)

        # Load the antibiotic drug master file to add names to the rows
        #drug_df = pd.read_csv(config['antibiotic_master'],sep='\t',low_memory=False)
        drug_df = pd.read_csv(config['antibiotic_master'],sep='\t',
            usecols=['code', 'name', 'name_alt', 'group_broad', 'group',
                     'class', 'subclass', 'card_class', 'combination'])
        # Ensure all drug names and classes are lowercase
        drug_df['name'] = drug_df['name'].astype(str).apply(lambda x: x.lower())
        drug_df['name_alt'] = drug_df['name_alt'].astype(str).apply(lambda x: x.lower())

        # Initialize a list to hold one dict for each drug; to be made into a df
        row_list = []

        # Get list of only drugs we made models for
        drugs = resf_val_drugs
        # For each drug, get the accuracy
        for drug in drugs:
            drug_name = drug_df.loc[drug_df['code'] == drug,'name'].item()
            drug_name_alt = drug_df.loc[drug_df['code'] == drug,'name_alt'].item()

            # Names of the lab amr column and the db prediction column in the df
            lab_col = drug+'_lab'
            pred_col = drug+'_pred'
            db_class_col = drug+'_db_class'

            # Filter dataframe for current drug
            #df = all_df[all_df['code']==drug]
            #df = all_df.set_index('id') # make id the index
            #df = df.filter(regex=drug) # filter for the drug's columns
            #print(df)
            #df = df.reset_index(drop=False) # add id index back as a column
            df = all_df.filter(items=['id',db_class_col,lab_col,pred_col])

            # Filter to keep only the samples which have
            # experimental/lab values for SIR
            df = df[df[lab_col].notnull()]

            # When the db doesn't find R / resistance, then mark it as S
            df[[pred_col]] = df[[pred_col]].fillna(value='S')

            # Number of samples with lab values for this drug
            n_samples = len(df['id'].tolist())
            # Number of samples that have a lab value of Resistant
            #n_lab_R = df[lab_col].value_counts()['R']
            val_counts = df[lab_col].value_counts()
            if 'R' in val_counts:
                n_lab_R = val_counts['R']
            else:
                n_lab_R=0

            # Make a column to indicate if the result matched the lab amr
            df[drug+'_correct'] = np.where(df[lab_col]==df[pred_col], 1, 0)
            # Count the correct predictions by summing the column
            correct_count = df[drug+'_correct'].sum()

            # Make a col to indicate if the result was R
            df[drug+'_predicted_R'] = np.where(df[pred_col]=='R', 1, 0)
            # Count the number of R predictions
            # In the case of the penicillins, they're all R
            predicted_R_count = df[drug+'_predicted_R'].sum()

            # What did the db classify the drug as? eg is db using penicillin,
            # penam, or beta-lactam to make the choice of R?
            db_class = df[db_class_col].tolist()
            # Remove duplicates
            db_class = list(set(db_class))
            # Remove nan
            #db_class = list(filter(lambda a: a != 'nan', db_class))
            db_class = [i for i in db_class if str(i) != 'nan']
            # Make into a string for the df
            db_class = ', '.join(db_class)

            temp = {
                'drug': drug,
                'name': drug_name,
                'alt_name': drug_name_alt,
                'db_class':db_class,
                'n_samples': n_samples,
                'n_correct': correct_count,
                'n_predicted_R': predicted_R_count,
                'n_lab_R': n_lab_R,
                'perc_correct': (correct_count/n_samples)*100,
                'perc_predictions_R': (predicted_R_count/n_samples)*100
                }

            row_list.append(temp)

        new_df = pd.DataFrame(row_list)

        pd.set_option('display.max_columns', 500)
        pd.set_option('display.width', 1000)

        new_df.to_csv(output.out, index=False, sep='\t')


################################################################################
# Model prediction & Accuracy
################################################################################

# get appropriate list of kmers
def get_union_kmer_file(wcs):
    return "results/wgs_{}/{}mer/union_kmers/union_merged.npy".format(wgs_type,kmer_len)

'''
def get_best_xgb_svm_model_all(wcs):
    # load the tsv of the best xgb models
    #best = "results/wgs_{w}/{k}mer/models/{g}/best_xgb_nested_cv_results.tsv".format(w=wcs.wgs_type,k=wcs.kmer_len,g=wcs.group)
    w = wgs_type
    g = 'all'
    k = str(config['kmer_size'])
    m = wcs.model_type
    best = "results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=m)
    df = pd.read_csv(best,sep='\t')
    # find the path to the best model for current pred_for
    path = df.loc[df['predfor'] == 'SIR'+wcs.pred_for, 'path'].values[0]
    mpath = path + model_type + "/"
    return mpath
'''

rule format_model_input:
    group: 'format_model_input'
    input:
        jelly = jellyfish(JELLY+"{id}.fa"),
        kmers = get_union_kmer_file,
    output:
        outfile = "{f}model_input/{id}.npy",
    run:
        # Load the kmers that need to be counted
        feats = np.load(input.kmers)
        # Make a column dict
        col_dict  = { feats[i] : i for i in range(0, len(feats))}
        # Initialize the row to hold the counts
        row = [0]*len(feats)
        # Open the input jellyfish file and count the kmers for std and krak
        with open(input.jelly) as fasta_file:
            for title, sequence in SimpleFastaParser(fasta_file):
                if sequence in col_dict:
                    row[col_dict[sequence]] = title
                else:
                    print(sequence, " not in standard")
        # Convert to numpy
        row = np.asarray(row, dtype=matrix_dtype)
        # Save
        np.save(output.outfile, row)



rule model_predict:
    group: 'model_predict'
    input:
        #i = "{f}model_input/{id}.npy"
        i = "results/wgs_standard/{predgroup}/model_input/{id}.npy"
    output:
      # o = "{f}model_output/{id}/{model_type}.tsv"
       o = "results/wgs_standard/{predgroup}/model_output/{id}/{model_type}.tsv/"
    benchmark:
        "benchmarks/compare_resfinder/{predgroup}/{model_type}/{id}.benchmark.txt"
    shell:
        """
        python scripts/model_creation/model_prediciton_help.py \
        --id {wildcards.id} \
        --infile {input.i} \
        --outfile {output.o} \
        --modeltype {wildcards.model_type} \
        -k {kmer_len}\
        -g munchkin
        """

'''
rule model_predict:
    group: 'model_predict'
    input:
        i = "{f}model_input/{id}.npy"
    output:
       o = "{f}model_output/{id}/{model_type}.tsv"
    run:
        id = wildcards.id
        in_file = input.i
        out_file = output.o
        model_type = wildcards.model_type
        # init the dict that will contain drug: S/I/R
        pred_dict ={'id':id}
        for drug in drugs:
            # get the path to the best xgb model for current drug
            g= 'munchkin'#'all'
            w = wgs_type
            k = str(config['kmer_size'])
            best = "results/wgs_{w}/{k}mer/models/{g}/best_{m}_nested_cv_results.tsv".format(w=w,k=k,g=g,m=model_type)
            df = pd.read_csv(best,sep='\t')
            # find the path to the best model for current pred_for
            path = df.loc[df['predfor'] == 'SIR'+drug, 'path'].values[0]
            mpath = path + model_type + "/"
            # get model and class dict for this model
            model_file = mpath + model_type + "_final_model.joblib"
            class_dict_file = mpath + "class_dict.pkl"
            # Load the class dict to decode the prediction
            with open(class_dict_file, 'rb') as fh:
                class_dict = pickle.load(fh)
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
            # Decode the prediction
            predn = decode_labels(predn, class_dict)

            # Add the prediction to the dict
            pred_dict[drug+'_'+model_type] = predn
        # convert to dataframe
        df = pd.DataFrame.from_dict(pred_dict)
        # save
        df.to_csv(out_file, index=False, sep='\t')
        """
        One tsv per input file
        id   | AMC_xgb | ...
        ERR1 |    R    | ...
        """
'''



rule cat_model_predictions:
    group: 'cat_model_predictions'
    input:
        expand("{{f}}model_output/{id}/{{mt}}.tsv",id=ids)
    output:
        "{f}summary_model/{mt}_predictions.tsv"
    run:
        cat = pd.concat([pd.read_csv(f, sep='\t') for f in input])
        #print(cat)
        cat.to_csv(output[0], index=False, sep='\t')

rule model_compare:
    group: 'model_compare'
    input:
        p = "{f}summary_model/{mt}_predictions.tsv"
    output:
        #out = "{f}summary_model/{mt}_acc_{standard}.tsv"
        "{f}summary_model/clsi_{mt}.tsv",
        "{f}summary_model/eucast_{mt}.tsv",
        "{f}summary_model/ecoff_{mt}.tsv",
    run:
        mt = wildcards.mt

        # Load the prjna metadata
        prjdf = pd.read_csv(resf_prj,sep='\t')

        # Load the prediction df
        #tempfile = "prediction_"+mt+".tsv"
        #tempfile = "results/wgs_standard/resf_val_data_prediction/raw_test/prediction_resfinder.tsv"
        #pred_df = pd.read_csv(tempfile,sep='\t')
        pred_df = pd.read_csv(input.p,sep='\t')

        # join the dataframes together
        prjdf = prjdf.set_index('id')
        pred_df = pred_df.set_index('id')
        pred_df = pred_df.join(prjdf)

        # 3 standards to look at
        standards = ['clsi','eucast','ecoff']

        for standard in standards:

            # Initialize a list to hold one dict for each drug;
            #to be made into a df
            row_list = []

            for drug in resf_val_drugs:

                # eucast has no data for NAL and TET so just make empty row
                if drug in ['NAL','TET','FOX'] and standard == 'eucast':
                    # row to add to df
                    temp = {
                        'drug': drug,
                        'n_samples': '',
                        'n_correct': '',
                        'n_predicted_R': '',
                        'n_lab_R': '',
                        'perc_correct': '',
                        'perc_predictions_R': ''
                        }
                else:
                    # i named things dumb
                    if mt == 'resfinder':
                        x = 'resf'
                    else:
                        x = mt

                    # Filter dataframe for current drug
                    df = pred_df[[drug+'_'+x,drug+'_'+standard]]

                    # Filter to keep only the samples which have
                    # experimental/lab values for SIR
                    df = df[df[drug+'_'+standard].notnull()]

                    #*******************************************************************
                    # Models for CIP, CTZ, GEN should be S/I/R. Remaining drugs
                    # have binary S/R classifiers. So need to bin all lab value
                    # Is into Rs. Only need to bin for CLSI and EUCAST because
                    # ECOFF doesn't create and Is.
                    if standard in ['clsi','eucast'] and drug in ['AMP','CHL','CTX','GEN','NAL','TET']:
                        #if drug=='CHL':
                        #    l = df[drug+'_'+standard].tolist()
                        #    print(l)
                        # bin all I standard classifications and I predictions into R
                        df.loc[df[drug+'_'+standard] == 'I', drug+'_'+standard] = 'R'
                        #if drug=='CHL':
                        #    l = df[drug+'_'+standard].tolist()
                        #    print(l)
                        #    sys.exit()

                    # Also ECOFF doesn't actually allow for I, so bin any I predictions
                    # into R predictions
                    # Need to bin I predictions into R for comparison if ecoff
                    if standard=='ecoff' and (mt == 'xgb' or mt == 'svm' or mt == 'ann'):
                        #print(df[drug+'_'+mt].tolist())
                        df.loc[df[drug+'_'+mt] == 'I', drug+'_'+mt] = 'R'
                        #print(df[drug+'_'+mt].tolist())
                        #sys.exit()
                    #*******************************************************************

                    # Number of samples with lab values for this drug
                    n_samples = df.shape[0]

                    #b = df[df[drug+'_'+standard].str.contains('R')]
                    if 'R' in df[drug+'_'+standard].tolist():
                        n_lab_R = df[drug+'_'+standard].value_counts()['R']
                    else:
                        n_lab_R = 0

                    # Make a column to indicate if the result matched the lab amr
                    df['correct_'+standard] = np.where(df[drug+'_'+standard]==df[drug+'_'+x], 1, 0)
                    # Count the correct predictions by summing the column
                    correct_count = df['correct_'+standard].sum()

                    '''
                    #print(df[drug+'_'+standard].tolist())
                    #print(df[drug+'_'+mt].tolist())
                    #print(df['correct_'+standard].tolist())
                    if drug == 'CHL' and standard == 'clsi':
                        print('chl')
                        df = df.reset_index()
                        a = df[drug+'_'+standard].tolist()
                        b = df[drug+'_'+mt].tolist()
                        c = df['correct_'+standard].tolist()
                        d = df['id'].tolist()
                        for i in range(len(a)):
                            print(d[i],a[i],b[i],c[i])
                        #sys.exit()
                    '''

                    # Make a col to indicate if the result was R
                    df['predicted_R_'+standard] = np.where(df[drug+'_'+x]=='R', 1, 0)
                    # Count the number of R predictions
                    # In the case of the penicillins, they're all R
                    predicted_R_count = df['predicted_R_'+standard].sum()

                    # row to add to df
                    temp = {
                        'drug': drug,
                        'n_samples': n_samples,
                        'n_correct': correct_count,
                        'n_predicted_R': predicted_R_count,
                        'n_lab_R': n_lab_R,
                        'perc_correct': (correct_count/n_samples)*100,
                        'perc_predictions_R': (predicted_R_count/n_samples)*100
                        }
                # append dict to list of rows
                row_list.append(temp)


            # turn list of dicts into df
            new_df = pd.DataFrame(row_list)

            #print(new_df)

            #sys.exit()
            interm_file = "{}summary_model/{}_{}.tsv".format(wildcards.f,standard,mt)
            print(interm_file)
            new_df.to_csv(interm_file,sep='\t',index=False)
            #new_df.to_csv(output.out,sep='\t',index=False)

################################################################################
# Combine the results into one table for each model type
################################################################################
'''
def get_per_sample_predictions(wcs):
    if wcs.mt == 'resfinder':
        return expand("{{f}}resf_output/{id}.tsv",id=ids)
    elif wcs.mt == 'xgb' or wcs.mt == 'svm':
        return expand("{{f}}model_output/{id}/{{mt}}.tsv",id=ids)
    else:
        raise Exception("check model type")

rule cat_predictions:
    group: 'cat_prediction_output'
    input:
        get_per_sample_predictions
    output:
        "{f}prediction_{mt}.tsv"
    run:
        cat = pd.concat([pd.read_csv(f, sep='\t') for f in input])
        #print(cat)
        cat.to_csv(output[0], index=False, sep='\t')
'''

################################################################################
# Compare results
################################################################################
"""
rule compare:
    input:
        p = "{f}prediction_{mt}.tsv"
    output:
        "{f}clsi_{mt}.tsv",
        "{f}eucast_{mt}.tsv",
        "{f}ecoff_{mt}.tsv",
        #PRED+"{mt}_allstandards.tsv"
    run:
        mt = wildcards.mt
        print(mt)


        # Load the prjna metadata
        prjdf = pd.read_csv(resf_prj,sep='\t')

        # Load the prediction df
        #tempfile = "prediction_"+mt+".tsv"
        #tempfile = "results/wgs_standard/resf_val_data_prediction/raw_test/prediction_resfinder.tsv"
        #pred_df = pd.read_csv(tempfile,sep='\t')
        pred_df = pd.read_csv(input.p,sep='\t')

        # join the dataframes together
        prjdf = prjdf.set_index('id')
        pred_df = pred_df.set_index('id')
        pred_df = pred_df.join(prjdf)

        # 3 standards to look at
        standards = ['clsi','eucast','ecoff']

        for standard in standards:

            # Initialize a list to hold one dict for each drug;
            #to be made into a df
            row_list = []

            for drug in resf_val_drugs:

                # eucast has no data for NAL and TET so just make empty row
                if drug in ['NAL','TET','FOX'] and standard == 'eucast':
                    # row to add to df
                    temp = {
                        'drug': drug,
                        'n_samples': '',
                        'n_correct': '',
                        'n_predicted_R': '',
                        'n_lab_R': '',
                        'perc_correct': '',
                        'perc_predictions_R': ''
                        }
                else:
                    # i named things dumb
                    if mt == 'resfinder':
                        x = 'resf'
                    else:
                        x = mt

                    # Filter dataframe for current drug
                    df = pred_df[[drug+'_'+x,drug+'_'+standard]]

                    # Filter to keep only the samples which have
                    # experimental/lab values for SIR
                    df = df[df[drug+'_'+standard].notnull()]

                    #*******************************************************************
                    # Models for CIP, CTZ, GEN should be S/I/R. Remaining drugs
                    # have binary S/R classifiers. So need to bin all lab value
                    # Is into Rs. Only need to bin for CLSI and EUCAST because
                    # ECOFF doesn't create and Is.
                    if standard in ['clsi','eucast'] and drug in ['AMP','CHL','CTX','GEN','NAL','TET']:
                        #if drug=='CHL':
                        #    l = df[drug+'_'+standard].tolist()
                        #    print(l)
                        # bin all I standard classifications and I predictions into R
                        df.loc[df[drug+'_'+standard] == 'I', drug+'_'+standard] = 'R'
                        #if drug=='CHL':
                        #    l = df[drug+'_'+standard].tolist()
                        #    print(l)
                        #    sys.exit()

                    # Also ECOFF doesn't actually allow for I, so bin any I predictions
                    # into R predictions
                    # Need to bin I predictions into R for comparison if ecoff
                    if standard=='ecoff' and (mt == 'xgb' or mt == 'svm' or mt == 'ann'):
                        #print(df[drug+'_'+mt].tolist())
                        df.loc[df[drug+'_'+mt] == 'I', drug+'_'+mt] = 'R'
                        #print(df[drug+'_'+mt].tolist())
                        #sys.exit()
                    #*******************************************************************

                    # Number of samples with lab values for this drug
                    n_samples = df.shape[0]

                    #b = df[df[drug+'_'+standard].str.contains('R')]
                    if 'R' in df[drug+'_'+standard].tolist():
                        n_lab_R = df[drug+'_'+standard].value_counts()['R']
                    else:
                        n_lab_R = 0

                    # Make a column to indicate if the result matched the lab amr
                    df['correct_'+standard] = np.where(df[drug+'_'+standard]==df[drug+'_'+x], 1, 0)
                    # Count the correct predictions by summing the column
                    correct_count = df['correct_'+standard].sum()

                    '''
                    #print(df[drug+'_'+standard].tolist())
                    #print(df[drug+'_'+mt].tolist())
                    #print(df['correct_'+standard].tolist())
                    if drug == 'CHL' and standard == 'clsi':
                        print('chl')
                        df = df.reset_index()
                        a = df[drug+'_'+standard].tolist()
                        b = df[drug+'_'+mt].tolist()
                        c = df['correct_'+standard].tolist()
                        d = df['id'].tolist()
                        for i in range(len(a)):
                            print(d[i],a[i],b[i],c[i])
                        #sys.exit()
                    '''

                    # Make a col to indicate if the result was R
                    df['predicted_R_'+standard] = np.where(df[drug+'_'+x]=='R', 1, 0)
                    # Count the number of R predictions
                    # In the case of the penicillins, they're all R
                    predicted_R_count = df['predicted_R_'+standard].sum()

                    # row to add to df
                    temp = {
                        'drug': drug,
                        'n_samples': n_samples,
                        'n_correct': correct_count,
                        'n_predicted_R': predicted_R_count,
                        'n_lab_R': n_lab_R,
                        'perc_correct': (correct_count/n_samples)*100,
                        'perc_predictions_R': (predicted_R_count/n_samples)*100
                        }
                # append dict to list of rows
                row_list.append(temp)


            # turn list of dicts into df
            new_df = pd.DataFrame(row_list)

            print(new_df)

            #sys.exit()
            interm_file = "{}{}_{}.tsv".format(wildcards.f,standard,mt)
            print(interm_file)
            new_df.to_csv(interm_file,sep='\t',index=False)

        '''
        clsi = pd.read_csv(PRED+mt+"_clsi.tsv",sep='\t', usecols=['drug','perc_correct'])
        clsi.rename(columns={'perc_correct':'clsi'})
        clsi.set_index('drug')

        eucast = pd.read_csv(PRED+mt+"_eucast.tsv",sep='\t', usecols=['drug','perc_correct'])
        eucast.rename(columns={'perc_correct':'eucast'})
        eucast.set_index('drug')

        ecoff = pd.read_csv(PRED+mt+"_ecoff.tsv",sep='\t', usecols=['drug','perc_correct'])
        ecoff.rename(columns={'perc_correct':'ecoff'})
        ecoff.set_index('drug')

        out_df = clsi.join(eucast)
        out_df = out_df.join(ecoff)

        out_df.to_csv(PRED+mt+"_allstandards.tsv",sep='\t',index=False)
        '''
"""
