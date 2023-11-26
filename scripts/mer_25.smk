import pandas as pd
import numpy as np
import math
#import pickle
import yaml
import scipy.sparse as sp

################################################################################
# Load the metadata and config files
################################################################################

# The main config file
configfile: "config/config.yaml"

# Load the metadata master and get the list of ids from it
master_path = config['metadata_master']
master = pd.read_csv(master_path,sep='\t')
# list of ids
ids = master['id'].tolist()

# Load the model config file
model_configfile = config['model_composition_configfile']
# Load the model compositon settings
with open(model_configfile, 'r') as fh:
    model_config = yaml.safe_load(fh)


################################################################################
# Settings
################################################################################

# ------------------------------------------------------------------------------
# 25mer
# ------------------------------------------------------------------------------

matrix_dtype = config['kmer_matrix_dtype']
model_names = ['raisin']
model_name = 'raisin'
model_configfile = config['model_composition_configfile']
# Load the model compositon settings
yaml = YAML(typ='safe')
with open(model_configfile, 'r') as f:
    model_config = yaml.load(f)
drugs_we_care_about = model_config[model_name]['drugs_we_care_about']
sir_drugs_we_care_about = ['SIR'+drug for drug in drugs_we_care_about]
data_column = model_config[model_name]['data_column']
data_criteria = model_config[model_name]['data_criteria']
kmer_len = model_config[model_name]['kmer_len']
nfeats = model_config[model_name]['num_feats']
model_types = model_config[model_name]['model_types']
n_top_feats = model_config[model_name]['num_top_feats']

param_configfile = config['param_fixed_configfile']
with open(param_configfile, 'r') as fh:
    param_config = yaml.safe_load(fh)


################################################################################

# These rules don't benefit from being submitted separately on the cluster
# Better to run them as part of the submission job
localrules: genome_names, label_rows # Matrix
localrules: model_cat_results #model_cat_results_macro # Model results
localrules: impt_feats, dense, query_top_feats # Top features


rule all:
    input:
    # --------------------------------------------------------------------------
    # k-mer counting
        # Jellyfish - kmer counting
            #expand("results/wgs_standard/{k}mer/jellyfish/{id}.fa", k=kmer_len, id=ids),
        # Kat - k-mer spectra
            #expand("results/wgs_standard/kat_kmer_spectra/{k}mer.png",
            #    k=kmer_len),
    # --------------------------------------------------------------------------
    # Matrix
        # Union kmers
            #expand("results/wgs_standard/{k}mer/union_kmers/union_merged.npy",
            #    k=kmer_len),
        # All in one matrix
            #expand("results/wgs_standard/{k}mer/matrix/matrix.npy",
            #    k=kmer_len),
        # Matrix Labels
            #expand("results/wgs_standard/{k}mer/matrix/row_labels/{p}.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # Matrix - one per drug
            #expand("results/wgs_standard/{k}mer/matrix/{p}/matrix.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # Matrix - no zeroes
            #expand("results/wgs_standard/{k}mer/matrix/{p}/matrix_vt.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # Withhold
            #expand("results/wgs_standard/{k}mer/matrix/{p}/X_withheld.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # Make the folds for training
            #expand("results/wgs_standard/{k}mer/matrix/{p}/folds/f1/filtered/x_train.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # Try sparse matrix
            #expand("results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse.npy",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
        # try loading
            #expand("results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse_blah.npz",
            #    k=kmer_len,
            #    p=sir_drugs_we_care_about),
    # ---------------- ----------------------------------------------------------
    # Models
        # Train models
            #expand("results/wgs_standard/{k}mer/models/{mn}/{p}/nested_cv/{f}feats/{mt}/final_model_results_macro.tsv",
             #   k = kmer_len,
             #   mn = model_names,
             #   p = sir_drugs_we_care_about,
             #   f = nfeats,
             #   mt = model_types),
        # Cat results
            #expand("results/wgs_standard/{k}mer/models/{mn}/all_nested_cv_results_{result_type}.tsv",
             #   k = kmer_len,
             #   mn = model_names,
             #   result_type = ['micro','macro']),
            #expand("results/wgs_standard/{k}mer/models/{mn}/all_nested_cv_results_macro.tsv",
            #    k = kmer_len,
            #    mn = model_names),
            #expand("results/wgs_standard/{k}mer/models/{mn}/all_nested_cv_results_micro.tsv",
            #    k = kmer_len,
            #    mn = model_names),
    # --------------------------------------------------------------------------
    # Feature Extraction
        # Extract top features
            #expand("results/wgs_standard/{k}mer/models/{mn}/{p}/nested_cv/{f}feats/{mt}/top_{t}_feats.tsv",
             #   k = kmer_len,
             #   mn = model_names,
             #   p = sir_drugs_we_care_about,
             #   f = nfeats,
             #   mt = ['xgb'],
             #   t = n_top_feats),
        # Top feature summary
            #expand("results/wgs_standard/{k}mer/top_feats/{mn}/top_feats_dense.tsv",
            #    k = kmer_len,
            #    mn = model_names),
        # Top feature query
            #expand("results/wgs_standard/{k}mer/top_feats/{mn}/jellyfish_query/{id}.txt",
            #    k = kmer_len,
            #    mn = model_names,
            #    id = ids),
        # Top feature frequency table
            expand("results/wgs_standard/{k}mer/top_feats/{mn}/feat_count_table.tsv",
                k = kmer_len,
                mn = model_names),
    # --------------------------------------------------------------------------
    # Feature Annotation
        # Annotate genomes with prokka
            #expand("results/wgs_standard/annotated_assemblies/{id}_df.pkl",
            #    id = ids),
        # Find hits (find the top kmers in the annotated assembly files)
            #expand("results/wgs_standard/{k}mer/models/{mn}/{p}/nested_cv/{f}feats/annotation/top_{t}_feats_hits.pkl",
            #    k = kmer_len,
            #    mn = model_names,
            #    p = sir_drugs_we_care_about,
            #    f = nfeats,
            #    t = n_top_feats),
        # Summarize the top hits
            #expand("results/wgs_standard/{k}mer/models/{mn}/{p}/nested_cv/{f}feats/annotation/brief_top_{t}_feats_hits_{p}.tsv",
            #    k = kmer_len,
            #    mn = model_names,
            #    p = sir_drugs_we_care_about,
            #    f = nfeats,
            #    t = n_top_feats),


################################################################################
# k-mer counting
################################################################################

# ------------------------------------------------------------------------------
# Jellyfish - k-mer counting
# ------------------------------------------------------------------------------

# count k-mers with Jellyfish
rule jelly_count:
  group: 'jelly_count'
  input:
    "data/genomes/asmbl/{id}/scaffolds.fasta"
  output:
    "results/wgs_standard/{k}mer/jellyfish/{id}.jf"
  threads:
    2
  shell:
    "jellyfish count -C -m {wildcards.k} -s 100M -t {threads} {input} -o {output}"

# dump the Jellyfish results as fasta files
rule jelly_dump:
    group: 'jelly_dump'
    input:
        "results/wgs_standard/{k}mer/jellyfish/{id}.jf"
    output:
        "results/wgs_standard/{k}mer/jellyfish/{id}.fa"
    shell:
        "jellyfish dump {input} > {output}"


# ******************************************************************************
# Kat - k-mer frequency - in progress
# ******************************************************************************
'''
jelly_input = ["'(data/genomes/asmbl/"+id+"/scaffolds.fasta)'" for id in ids]


# Get the kmer counts across the whole dataset
# For kmer distribution
rule kat:
    group: 'kat'
    input:
    #use the raw fastq
        ?.fastq.gz
    # jellyfish jf
        #"results/wgs_standard/{k}mer/jellyfish/{id}.jf"
    # assembly fasta
        #expand("data/genomes/asmbl/{id}/scaffolds.fasta",id=ids)
    output:
        "results/wgs_standard/{k}mer/kat_kmer_spectra/{id}.png"
    params:
        prefix = "results/wgs_standard/kat_kmer_spectra/{k}mer",
        input_quoted = lambda wildcards, input: [ f"'{a}'" for a in input ],
        input_all = "results/wgs_standard/{k}mer/jellyfish/?.jf"
    threads:
        2
    shell:
        "kat hist --threads {threads} --mer_len {wildcards.k} --verbose --output_type png --output_prefix {params.prefix} {input}"
'''
'''
"kat hist --threads {threads} --mer_len {wildcards.k} --verbose --output_type png --output_prefix {params.prefix} ({params.input_quoted})"


"kat hist --threads 1 --mer_len 11 --verbose --output_type png --output_prefix kattest results/wgs_standard/11mer/jellyfish/?.jf"
'''
# --canonical

'''
sbatch -c 2 --mem 4G -o kat_test.out -J kat_test --wrap="kat hist --output_prefix kat_test results/wgs_standard/25mer/jellyfish/*.jf"

sbatch -c 2 --mem 4G -o kat_test.out -J kat_test --wrap="kat hist --verbose --mer_len 25 --output_prefix kat_test data/genomes/asmbl/*/scaffolds.fasta"

"kat hist --output_prefix kat_test 'results/wgs_standard/25mer/jellyfish/*.jf'"

"kat hist --verbose --mer_len 25 --output_prefix kat_test data/genomes/asmbl/*/scaffolds.fasta"

'''

################################################################################
# Matrix building
################################################################################

# ------------------------------------------------------------------------------
# Union k-mers: Determine the set (union) of kmers present in dataset
# ------------------------------------------------------------------------------

"""
For kmers longer than 11, it is more memory efficient to precomute the present
kmers instead of computing all possible kmers from ATCG. Also more efficient to
parse the files in parallel then combine them.
"""

## Math
# For resource efficiency, we need to split the matrix creation into chunks.
# The ideal number of genomes per chunk is 128 because of waffles capacity
# and how the chunk scripts are set up
num_per_split = 128
num_genomes = len(ids)
num_splits = math.ceil(num_genomes/num_per_split)
splits = [str(i) for i in range(1,num_splits+1)]


# Find the union of kmers in batches (splits)
rule union_kmers:
    group: 'union_kmers'
    input:
        #jellyfish(expand(JELLY+"{id}.fa", id=ids))
        expand("results/wgs_standard/{{k}}mer/jellyfish/{id}.fa", id=ids)
    output:
        "results/wgs_standard/{k}mer/union_kmers/union_set_{split}.npy"
    params:
        #real_in = JELLY+"*.fa"
        real_in = "results/wgs_standard/{k}mer/jellyfish/*.fa"
    shell:
        "python scripts/matrix_creation/union_kmers.py \
        -k {wildcards.k} \
        -o {output} \
        -s {wildcards.split} \
        -ms {num_splits} \
        -i {params.real_in}"

# Merge the batches (splits) into one
rule merge_union_kmers:
    group: 'merge_union_kmers'
    input:
        expand("results/wgs_standard/{{k}}mer/union_kmers/union_set_{split}.npy", split=splits)
    output:
        "results/wgs_standard/{k}mer/union_kmers/union_merged.npy"
    shell:
        "python scripts/matrix_creation/union_kmers_merge.py \
        -k {wildcards.k} \
        -o {output} \
        -i {input}"


# ------------------------------------------------------------------------------
# Build k-mer frequency matrix
# ------------------------------------------------------------------------------

# Will be the order of genomes in the matrix (order is the list of ids)
rule genome_names:
    group: 'genome_names'
    output:
        "results/wgs_standard/{k}mer/matrix/genome_names.npy"
    run:
        np.save(output[0], ids)

# Make matrix chunks
rule make_matrix:
    group: 'matrix'
    input:
        union_merged = "results/wgs_standard/{k}mer/union_kmers/union_merged.npy",
        genome_names = "results/wgs_standard/{k}mer/matrix/genome_names.npy",
        #infiles = jellyfish(expand(JELLY+"{id}.fa", id=ids)),
        infiles = expand("results/wgs_standard/{{k}}mer/jellyfish/{id}.fa", id=ids)
    output:
        matrix = "results/wgs_standard/{k}mer/matrix/splits/matrix_split_{split}.npy",
        rows = "results/wgs_standard/{k}mer/matrix/splits/matrix_split_{split}_rows.npy",
        cols ="results/wgs_standard/{k}mer/matrix/splits/matrix_split_{split}_cols.npy",
        #matrix2 = MATRIX+"splits/matrix_split_{split}.df"
    params:
        jelly_path = "results/wgs_standard/{k}mer/jellyfish/",
        outdir = "results/wgs_standard/{k}mer/matrix/splits/"
    shell:
        "python scripts/matrix_creation/matrix_chunks.py \
        -i {params.jelly_path} \
        -o {params.outdir} \
        -d {matrix_dtype} \
        -g {input.genome_names} \
        -uk {input.union_merged} \
        -s {wildcards.split} \
        -ms {num_splits} \
        -m {master_path}"

rule merge_matrices:
    group: 'merge_matrices'
    input:
        union_merged = "results/wgs_standard/{k}mer/union_kmers/union_merged.npy",
        genome_names = "results/wgs_standard/{k}mer/matrix/genome_names.npy",
        matrices = expand(
                "results/wgs_standard/{{k}}mer/matrix/splits/matrix_split_{split}.npy",
                split=splits),
        rows = expand(
                "results/wgs_standard/{{k}}mer/matrix/splits/matrix_split_{split}_rows.npy",
                split=splits),
        cols = expand(
                "results/wgs_standard/{{k}}mer/matrix/splits/matrix_split_{split}_cols.npy",
                split=splits),
    params:
        outdir = "results/wgs_standard/{k}mer/matrix/"
    output:
        "results/wgs_standard/{k}mer/matrix/matrix.npy",
        "results/wgs_standard/{k}mer/matrix/rows.npy",
        "results/wgs_standard/{k}mer/matrix/cols.npy"
    shell:
        "python scripts/matrix_creation/matrix_chunks_merge.py \
        -i {input.matrices} \
        -r {input.rows} \
        -c {input.cols} \
        -o {params.outdir} \
        -g {input.genome_names} \
        -uk {input.union_merged} \
        -ms {num_splits}"


# ------------------------------------------------------------------------------
# Label Rows
# ------------------------------------------------------------------------------

# turns every column, except 'id', in the metadata master into labels for
# the matrix
rule label_rows:
    group: 'label_rows'
    input:
        rows = "results/wgs_standard/{k}mer/matrix/rows.npy"
    output:
        labels = "results/wgs_standard/{k}mer/matrix/row_labels/{pred_for}.npy"
    params:
        outdir = "results/wgs_standard/{k}mer/matrix/row_labels/"
    run:
        master_path = config['metadata_master']
        tmp_master = pd.read_csv(master_path,sep='\t',low_memory=False)
        # thing to predict for
        pred_for = wildcards.pred_for
        # underscores mess with the regex in snakemake later, so remove them
        tmp_master.columns = tmp_master.columns.str.replace('_', '')
        # genome names in order they appear in the matrix
        genome_names = np.load(input.rows)
        # make the index the IDs column
        tmp_master = tmp_master.set_index('id')
        # reorder the metadata sheet to match the matrix
        tmp_master = tmp_master.reindex(genome_names)
        # If column is empty raise error
        if len(tmp_master[pred_for].value_counts()) <= 0:
            import warnings
            raise Exception("column {} has no data".format(pred_for))
        # convert column to numpy
        labels = tmp_master[pred_for].to_numpy()
        # Save without encoding
        np.save(output.labels, labels)


################################################################################
# Matrix Manipulation
################################################################################

# ------------------------------------------------------------------------------
# Matrix Reduction/Filtering
# ------------------------------------------------------------------------------
"""
For kmer_len > 11, it is better to have 1 matrix per antibiotic, to reduce its
size (and therefore RAM usage)

Also kmer len >11 means less likely to occur by chance, so kmers may be absent
in many genomes. So after making 1 matrix per antibiotic, also remove any
all-zero columns to reduce size further.
"""


# Make one matrix per drug
rule one_matrix_per_drug:
    group: 'one_m_per'
    input:
        matrix = "results/wgs_standard/{k}mer/matrix/matrix.npy",
        rows = "results/wgs_standard/{k}mer/matrix/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/row_labels/{p}.npy",
        cols = "results/wgs_standard/{k}mer/matrix/cols.npy"
    output:
        temp("results/wgs_standard/{k}mer/matrix/{p}/matrix.npy"),
        rows = "results/wgs_standard/{k}mer/matrix/{p}/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/{p}/labels.npy",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols.npy"
    params:
        outdir = "results/wgs_standard/{k}mer/matrix/{p}/",
        data_column = data_column, #model_config[model_name]['data_column'],
        #num_feats_for_chunk = model_config[model_name]['num_feats'],
        data_criteria = data_criteria #model_config[model_name]['data_criteria']
    shell:
        "python scripts/model_creation/one_matrix_per_drug.py \
        --matrix {input.matrix} \
        --rows {input.rows} \
        --labels {input.labels} \
        --cols {input.cols} \
        --metadata {master_path} \
        --predfor {wildcards.p} \
        --outdir {params.outdir} \
        --forcecol {params.data_column} \
        --criteria {params.data_criteria}"

#--modelname {model_name} \
#--feats {params.num_feats_for_chunk} \

# Remove columns for which every entry is 0
rule var_thresh:
    group: 'rm0_thresh'
    input:
        matrix = "results/wgs_standard/{k}mer/matrix/{p}/matrix.npy",
        rows = "results/wgs_standard/{k}mer/matrix/{p}/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/{p}/labels.npy",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols.npy"
    output:
        matrix = "results/wgs_standard/{k}mer/matrix/{p}/matrix_vt.npy",
        #rows = "results/wgs_standard/{k}mer/matrix/{p}/rows_vt.npy",
        #labels = "results/wgs_standard/{k}mer/matrix/{p}/labels_vt.npy",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols_vt.npy"
    params:
        outdir = "results/wgs_standard/{k}mer/matrix/{p}/"
    shell:
        "python scripts/model_creation/var_thresh.py \
        --matrix {input.matrix} \
        --cols {input.cols} \
        --outdir {params.outdir} "

#--modelname {model_name} \


# ------------------------------------------------------------------------------
# Withhold a piece for validation after training
# ------------------------------------------------------------------------------

rule withhold:
    group: 'withheld11'
    input:
        matrix = "results/wgs_standard/{k}mer/matrix/{p}/matrix_vt.npy",
        rows = "results/wgs_standard/{k}mer/matrix/{p}/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/{p}/labels.npy",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols_vt.npy",
    output:
        x = "results/wgs_standard/{k}mer/matrix/{p}/X.npy",
        x_withheld = "results/wgs_standard/{k}mer/matrix/{p}/X_withheld.npy",
        y = "results/wgs_standard/{k}mer/matrix/{p}/Y.npy",
        y_withheld = "results/wgs_standard/{k}mer/matrix/{p}/Y_withheld.npy",
        z = "results/wgs_standard/{k}mer/matrix/{p}/Z.npy",
        z_withheld = "results/wgs_standard/{k}mer/matrix/{p}/Z_withheld.npy",
    params:
        outdir = "results/wgs_standard/{k}mer/matrix/{p}/",
        #data_column = data_column #model_config[model_name]['data_column'],
        #num_feats_for_chunk = model_config[model_name]['num_feats'],
        #data_criteria = data_criteria #model_config[model_name]['data_criteria']
    shell:
        "python scripts/model_creation/presplit.py \
        --matrix {input.matrix} \
        --rows {input.rows} \
        --labels {input.labels} \
        --cols {input.cols} \
        --metadata {master_path} \
        --predfor {wildcards.p} \
        --outdir {params.outdir} "

#--modelname {model_name} \
#--feats {params.num_feats_for_chunk} \
#--forcecol {params.data_column} \
#--criteria {params.data_criteria}

# ------------------------------------------------------------------------------
# Convert to a sparse matrix
# ------------------------------------------------------------------------------
rule make_sparse:
    group:"make_sparse"
    input:
        x = "results/wgs_standard/{k}mer/matrix/{p}/X.npy",
        #y = "results/wgs_standard/{k}mer/matrix/{p}/Y.npy",
        #z = "results/wgs_standard/{k}mer/matrix/{p}/Z.npy",
    output:
        x = "results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse.npz",
    run:
        import numpy as np
        import scipy.sparse as sp

        x_train = np.load(input.x)
        x_new = sp.csr_matrix(x_train)

        sp.save_npz(output.x,x_new,compressed=True)


rule test_npz_load:
    input:
        x = "results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse.npz",
    output:
        x = "results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse_blah.npz"
    run:
        #import scipy.sparse as sp

        print("---starting load---", flush=True)

        X_path = input.x

        print(X_path)

        X = sp.load_npz(X_path)

        print(type(X))

        sys.exit()




################################################################################
# Model Training
################################################################################


rule model_training:
    group: 'param_tuning'
    input:
        matrix = "results/wgs_standard/{k}mer/matrix/{p}/matrix.npy",
        rows = "results/wgs_standard/{k}mer/matrix/{p}/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/{p}/labels.npy",
        #cols = "results/wgs_standard/{k}mer/matrix/cols.npy",
        #x = "results/wgs_standard/{k}mer/matrix/{p}/X.npy",
        x = "results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse.npz",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols_vt.npy",
        #x_fal_filt = "results/wgs_standard/{k}mer/matrix/{p}/validation_filtered/X_filtered.npy",
    output:
        final_results = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/final_model_results_micro.tsv",
        final_results2 = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/final_model_results_macro.tsv",
        #feats = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/feature_importance.npy"
    benchmark:
        "benchmarks/param_tuning/{k}mer/{model_name}/{p}_{f}_{model_type}.benchmark.txt"
    params:
        matrix = "results/wgs_standard/{k}mer/matrix/{p}/matrix.npy",
        rows = "results/wgs_standard/{k}mer/matrix/{p}/rows.npy",
        labels = "results/wgs_standard/{k}mer/matrix/{p}/labels.npy",
        #cols = "results/wgs_standard/{k}mer/matrix/cols.npy",
        #x = "results/wgs_standard/{k}mer/matrix/{p}/X.npy",
        x = "results/wgs_standard/{k}mer/matrix/{p}/sparse_test/X_sparse.npz",
        cols = "results/wgs_standard/{k}mer/matrix/{p}/cols_vt.npy",
        #x_fal_filt = "results/wgs_standard/{k}mer/matrix/{p}/validation_filtered/X_filtered.npy",
        param_grid_file = config['param_grid_configfile'],
        param_fixed_file = config['param_fixed_configfile'],
        model_config = config['model_composition_configfile'],
        outdir = "results/wgs_standard/{k}mer/models/{model_name}/",
        data_column = lambda wcs: model_config[wcs.model_name]['data_column'],
        data_criteria = lambda wcs: model_config[wcs.model_name]['data_criteria'],
        param_tuning = lambda wcs: model_config[wcs.model_name]['param_tuning'],
        verbose = 0,
        tf_inter_threads=2,
        tf_intra_threads=lambda wcs: param_config['resource_settings'][wcs.model_type]['n_jobs'],
        python_warn='ignore' # 'ignore' or 'default'
    shell:
        "python scripts/model_creation/25mer_train.py  \
                        --modelname {wildcards.model_name} \
                        --paramconfig {params.param_fixed_file} \
                        --matrix {params.matrix} \
                        --rows {params.rows} \
                        --cols {params.cols} \
                        --labels {params.labels} \
                        --feats {wildcards.f} \
                        --verbose {params.verbose} \
                        --modeltype {wildcards.model_type} \
                        --metadata {master_path} \
                        --predfor {wildcards.p} \
                        --outdir {params.outdir} \
                        --forcecol {params.data_column} \
                        --criteria {params.data_criteria}"


# ------------------------------------------------------------------------------
# Model Combine Results
# ------------------------------------------------------------------------------

# combine the results into one file
rule model_cat_results:
    group:'model_cat_results'
    input:
        #model_cat_results_input_micro
        expand("results/wgs_standard/{{k}}mer/models/{{model_name}}/{p}/nested_cv/{f}feats/{mt}/final_model_results_{{type}}.tsv",
            p = sir_drugs_we_care_about,
            f = nfeats,
            mt = model_types)
    output:
        all_results = "results/wgs_standard/{k}mer/models/{model_name}/all_nested_cv_results_{type}.tsv"
    run:
        # Create one df of all results
        df = pd.read_csv(input[0], sep='\t',low_memory=False)
        for file in input[1:]:
            tmp = pd.read_csv(file, sep='\t',low_memory=False)
            df = pd.concat([df,tmp])
        df.to_csv(output.all_results, sep='\t', index=False)
        df.to_csv( output.all_results, index=False, sep='\t' )

################################################################################
# Feature Extraction
################################################################################

# get top features for each model trained
rule impt_feats:
    group: 'xgb_extract'
    input:
        #flagfile = "results/wgs_standard/{k}mer/models/{model_name}/all_nested_cv_results_macro.tsv",
        #flag = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/final_model_results_micro.tsv",
    output:
        out_tsv = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/top_{n_top_feats}_feats.tsv",
        out_npy = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/top_{n_top_feats}_feats.npy"
    params:
        xgb_feats = "results/wgs_standard/{k}mer/models/{model_name}/{p}/nested_cv/{f}feats/{model_type}/feature_importance.npy"
    run:
        import numpy as np
        #"python feature_annotation/get_top_impt_feats.py -i {params.imptfeat} -o {output} -n {wildcards.num_top_feats}"
        # Load the feature importances of the saved xgb model
        feat_imps = np.load(params.xgb_feats, allow_pickle=True)
        # Turn the feature importances into a datframe; name columns appropriately
        df = pd.DataFrame(feat_imps.T)
        df = df.rename(columns={0:'feat',1:'impt'})
        # Drop features with an importance of 0 (not important features)
        df = df[df.impt != '0.0']
        # Convert the df to float
        df['impt'] = df['impt'].astype(str).astype(float)
        # Extract the top X most important features
        top = df.nlargest(int(wildcards.n_top_feats), 'impt', keep='first')
        # Add column with drug name, for easier stitching together later
        top['drug'] = wildcards.p
        # Add column with num feats, for easier stitching together later
        top['numfeats'] = wildcards.f
        # Save a tsv of the top feats and their importances
        top.to_csv(output.out_tsv,sep='\t',index=False)
        # Save a numpy file of the top feats
        top_npy = top['feat'].to_numpy()
        np.save(output.out_npy, top_npy)

# make a single summary table
rule dense:
    input:
        top_feat_files = expand("results/wgs_standard/{{k}}mer/models/{{model_name}}/{p}/nested_cv/{f}feats/{model_type}/top_{n_top_feats}_feats.tsv",
            p=sir_drugs_we_care_about,
            f=nfeats,
            model_type=model_types,
            n_top_feats=n_top_feats)
    output:
        out_tsv = "results/wgs_standard/{k}mer/top_feats/{model_name}/top_feats_dense.tsv",
        out_npy = "results/wgs_standard/{k}mer/top_feats/{model_name}/top_feats_order.npy",
    #"tmp.txt"
    run:
        # read each file and stick it into 1 dataframe
        df = pd.concat( [ pd.read_csv( f, sep='\t' ) for f in input.top_feat_files ] )
        # save it
        df.to_csv( output.out_tsv, index=False, sep='\t' )
        # also make an np
        df_npy = df['feat'].to_numpy()
        np.save( output.out_npy, df_npy )


# make a frequency matrix of these top kmers
# query the kmers in each .jf file
rule query_top_feats:
    input:
        jelly_file = "results/wgs_standard/{k}mer/jellyfish/{id}.jf",
        feat_order = "results/wgs_standard/{k}mer/top_feats/{model_name}/top_feats_order.npy"
    output:
        "results/wgs_standard/{k}mer/top_feats/{model_name}/jellyfish_query/{id}.txt"
    run:
        feat_list = np.load(input.feat_order, allow_pickle=True)
        feat_list_formatted = " ".join(feat_list)
        print(feat_list_formatted)
        shell("jellyfish query {input.jelly_file} {feat_list_formatted} > {output}")


# make a frequency matrix of these top kmers
# faster than trimming the master matrix ?
rule top_feat_table:
    input:
        query_files = expand("results/wgs_standard/{{k}}mer/top_feats/{{model_name}}/jellyfish_query/{id}.txt", id=ids),
        feat_order = "results/wgs_standard/{k}mer/top_feats/{model_name}/top_feats_order.npy"
    output:
        outfile = "results/wgs_standard/{k}mer/top_feats/{model_name}/feat_count_table.tsv"
    run:
        # goal is a table of the format:

        #  genome     kmer1     kmer2    kmer 3
        #  ERR123     11        0         2

        # order of the columns
        #feat_order = np.load(input.feat_order,allow_pickle=True)

        final_df = pd.DataFrame()

        # for each file
        for f in input.query_files:
            # get id name from the file path
            id = f.split("/")[-1].split(".txt")[0]
            #id = id.split(".txt")[0]
            #print("id",id)
            #print('-----------')

            # load the tsv, make the index the first column (kmers)
            jf = pd.read_csv(f,header=None,sep=' ',index_col=0)
            #jf.set_index()
            #print(jf)
            #print('-----------')

            # transpose it so that kmers are the columns
            #print('-----------transpose')
            jf = jf.T
            #print(jf)

            # add the name to the row
            #print('-----------add id')
            new_col = [id]
            jf.insert(loc=0,column="id",value=new_col)
            #print(jf)

            #print('-----------concat')
            final_df = pd.concat([final_df,jf])
            #print(final_df)

        print('-----------final df')
        print(final_df)
        final_df.to_csv(output.outfile,sep="\t",index=False)




################################################################################
# Feature Annotation
################################################################################

# ------------------------------------------------------------------------------
# annotate assemblies
# ------------------------------------------------------------------------------

# Annotate every genome
rule prokka:
    group: 'prokka'
    input:
        "data/genomes/asmbl/{id}/scaffolds.fasta"
    output:
        ffn = "results/wgs_standard/annotated_assemblies/{id}/{id}.ffn",
        gff = "results/wgs_standard/annotated_assemblies/{id}/{id}.gff"
    params:
        outdir = "results/wgs_standard/annotated_assemblies/{id}",
        cpus = 8#workflow.cores
    run:
        shell("export OMP_NUM_THREADS={params.cpus} && \
        export OPENBLAS_NUM_THREADS={params.cpus} && \
        prokka {input} \
        --outdir {params.outdir} \
        --prefix {wildcards.id} \
        --cpus {params.cpus} \
        --force \
        --compliant")

# Turn the results into a pandas df and save with pickle
rule prokka_to_df:
    group: 'prokka'
    input:
        "results/wgs_standard/annotated_assemblies/{id}/{id}.gff"
    output:
        "results/wgs_standard/annotated_assemblies/{id}_df.pkl"
    run:
        shell("export OPENBLAS_NUM_THREADS=1")
        import gffpandas.gffpandas as gffpd
        anno = gffpd.read_gff3(input[0])
        ### consider cleaning dataframe here ###
        anno.df.to_pickle(str(output),protocol=4)


# ------------------------------------------------------------------------------
# blast
# ------------------------------------------------------------------------------
rule makeblastdb_setup:
    group: 'blastdb_setup'
    input:
        expand("data/genomes/asmbl/{id}/scaffolds.fasta", id=ids)
    output:
        "results/wgs_standard/blast_db_input/all_multifasta.fasta"
    run:
        from Bio import Seq, SeqIO
        fasta_out = str(output[0])
        for file in input:
            fasta_in = str(file)
            id = fasta_in.rsplit("/",1)[0].split("/")[-1]
            #print(id)
            with open(fasta_in) as fin, open(fasta_out, 'a') as fout:
                for record in SeqIO.parse(fin, "fasta"):
                    contig_seq = record.seq
                    contig_id = record.id
                    contig_desc = record.description
                    contig_num = contig_desc.split("NODE_")[-1].split("_length")[0]
                    record.id = "{}_{}".format(id,contig_num)
                    #print(contig_seq)
                    #print(contig_id)
                    #print(contig_desc)
                    #print(contig_num)
                    #print(record.id)
                    SeqIO.write(record, fout, "fasta")

# make a blast db out of the assemblies
rule makeblastdb:
    group: 'makeblastdb'
    input:
        "results/wgs_standard/blast_db_input/all_multifasta.fasta"
    output:
        "results/wgs_standard/asmbl_blastdb/asmbl_blastdb.nal"
    #"results/wgs_standard/asmbl_blastdb/asmbl_blastdb.nhr"
    params:
        db = "results/wgs_standard/asmbl_blastdb/asmbl_blastdb"
    shell:
        "makeblastdb -in {input} -parse_seqids -dbtype nucl -out {params.db}"

# needs the top x features extracted from the xgb model along with the blastdb
rule blast_top_feats:
    group: 'blast'
    input:
        #imptfeats = model("results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/xgb/top_{num_top_feats}_feats.npy"),
        imptfeats = "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/xgb/top_{num_top_feats}_feats.npy",
        blastdb = "results/wgs_standard/asmbl_blastdb/asmbl_blastdb.nal"
    #blastdb = "results/wgs_standard/asmbl_blastdb/asmbl_blastdb.nhr"
    output:
        "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/top_{num_top_feats}_feats_blast.txt"
    params:
        db = "results/wgs_standard/asmbl_blastdb/asmbl_blastdb"
    run:
        alt_out = str(output).split('.')[0]+'.query'
        top_feats = np.load(input.imptfeats, allow_pickle=True)
        if(isinstance(top_feats[0],bytes)):
            top_feats = [i.decode('utf-8') for i in top_feats]
        assert(len(top_feats[0])==int(wildcards.kmer_len))
        with open(alt_out,'a') as fh:
            for feat in top_feats:
                fh.write(">{}\n".format(feat))
                fh.write(feat+"\n")
        shell("blastn -task blastn-short -db {params.db} -query {alt_out} -ungapped -perc_identity 100 -dust no -word_size {kmer_len} -max_target_seqs 50000 -evalue 100000 -outfmt 6 -out {output}")

#if kmer_len == "11":
#    alt_out = str(output).split('.')[0]+'.query'
#    top_feats = np.load(input.imptfeats, allow_pickle=True)
#    if(isinstance(top_feats[0],bytes)):
#        top_feats = [i.decode('utf-8') for i in top_feats]
#    assert(len(top_feats[0])==11)
#    with open(alt_out,'a') as fh:
#        for feat in top_feats:
#            fh.write(">{}\n".format(feat))
#            fh.write(feat+"\n")
#    shell("blastn -task blastn-short -db {params.db} -query {alt_out} -ungapped -perc_identity 100 -dust no -word_size 11 -max_target_seqs 50000 -evalue 100000 -outfmt 6 -out {output}")
#shell("blastn -task blastn-short -db data/master.db -query {alt_out} -ungapped -perc_identity 100 -dust no -word_size 11 -max_target_seqs 50000 -evalue 100000 -outfmt 6 -out {output}")


# ------------------------------------------------------------------------------
# find hits & condense the table for better reading
# ------------------------------------------------------------------------------

# if expanding, use {{ }} so it doesnt expand on stuff we dont need ?
# need to check input, and the script itself
rule find_hits:
    group: 'hit'
    input:
        blast_path = "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/top_{num_top_feats}_feats_blast.txt",
        prokka_files = expand("results/wgs_standard/annotated_assemblies/{id}_df.pkl", id=ids),
        #top_f = model("results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/xgb/top_{num_top_feats}_feats.npy"),
        top_f = "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/xgb/top_{num_top_feats}_feats.npy"
    output:
        p = "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/top_{num_top_feats}_feats_hits.pkl",
        t = "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/top_{num_top_feats}_feats_hits.tsv"
    params:
        prokka_loc = "results/wgs_standard/annotated_assemblies/"
    shell:
        "python scripts/feature_annotation/find_hits.py \
        -blast_path {input.blast_path} \
        -top_f {input.top_f} \
        -prokka_loc {params.prokka_loc} \
        -kmer_len {kmer_len} \
        -out {output.p}"

# condense identical entries and add a count of them
rule summarize_hits_brief:
    group: 'brief_sum'
    input:
        "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/top_{num_top_feats}_feats_hits.pkl"
    output:
        "results/wgs_standard/{kmer_len}mer/models/{group}/{pred_for}/nested_cv/{num_feats}feats/annotation/brief_top_{num_top_feats}_feats_hits_{pred_for}.tsv"
    run:
        indf = pd.read_pickle(input[0])
        tmp = indf[['kmer', 'gene_up', 'genome_id', 'dist_up', 'gene_down', 'dist_down']].groupby(['kmer', 'gene_up', 'dist_up', 'gene_down', 'dist_down']).agg(['count'])
        tmp.columns = tmp.columns.droplevel(0)
        tmp.to_csv(output[0],sep='\t')


# ------------------------------------------------------------------------------
# CARD
# ------------------------------------------------------------------------------

rule card:
    input:

    output:

    shell:
