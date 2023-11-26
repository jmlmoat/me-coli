
import pandas as pd
import numpy as np
import pickle as pkl
import json


# stop numpy from being a greedy bastard
import os
os.environ["OMP_NUM_THREADS"] = "1"

configfile: "config/config.yaml"

# Load the metadata master
master = pd.read_csv(config['metadata_master'],sep='\t')

# Settings
wgs_type = config['wgs_type'] # standard or kraken-filtered
kmer_len = config['kmer_size']

'''
# New master file - create new file for wgs type so no overwriting
new_master = config['metadata_master'].split('.tsv')[0]+'_'+wgs_type+'.tsv'
'''

# All non-reference ids
ids = master['id'].tolist()

# Path to resfinder location; will be downloaded if not present
resfinder_path = config['resfinder_path']
if '~' in resfinder_path:
    resfinder_path = os.path.expanduser(resfinder_path)
    # the file to run resfinder
    resfinder = resfinder_path+"run_resfinder.py"
# Path to resfinder db location; will be downloaded if not present
#resfinder_db_path = config['resfinder_db_path']
#if '~' in resfinder_db_path:
#    resfinder_db_path = os.path.expanduser(resfinder_db_path)
    #print(resfinder_db_path)

# Paths
OUT = "results/wgs_{}/external_amr_tools/".format(wgs_type)
ABR = OUT+"abricate/"
RGI = OUT+"rgi/"
AMF = OUT+"amrfinder/"
RES = OUT+"resfinder/"
RAW = "data/genomes/raw/" # raw genomes

# Flag to show if amrfinder was updated
AMF_FLAG = OUT+"update_amrfinder.done".format(wgs_type)
#ABR_FLAG = OUT+"update_{db}.done"


# ALL dbs used; must update the abr_amr_dbs list as well below
#dbs = ['rgi','amrfinder','card', 'ncbi', 'resfinder', 'cgeresfinder']
dbs = ['rgi','amrfinder','card', 'ncbi', 'resfinder', 'rawcgeresfinder','asmblcgeresfinder']
#*******************************************************************************
# All DBs used by ABRicate; can just spefify additional dbs like vfdb if want
abr_dbs = ['card', 'ncbi', 'resfinder', 'plasmidfinder', 'vfdb']
# ABRicate dbs used for AMR (card, ncbi, resfinder)
abr_amr_dbs = ['card', 'ncbi', 'resfinder']
#abr_dbs = ['argannot']
#*******************************************************************************

# ******************************************************************************
# Testing
#dbs = ['rawcgeresfinder','asmblcgeresfinder']
#dbs = ['cgeresfinder']
#dbs = ['card']
#ids = ['W98_2']
# ResFinder 4.0 paper's 390 genome val set

#predgroup = 'resf_val_data'
#resf_prj = "metadata/public/PRJNA6164_metadata.tsv"
#resf_df = pd.read_csv(resf_prj, sep='\t')
#ids = resf_df['id'].tolist()

#ids = [ 'EC2013AMR1568' ]

#******************************************************************************


rule all:
    input:
    # ABRicate all dbs incl vfdb etc
        #expand(ABR+"{id}/{db}.tsv",id=ids, db=abr_dbs)
    # --------------------------------------------------------------------------
    ## Resistance Lists
    # RGI resistance list
        #expand(rgi = RGI+"{id}/{id}.tsv",id=ids)
        #expand(RGI+"{id}/rgi_resistance.pkl",id=ids),
    # AMRFinder resistance list
        #expand(AMF+"{id}/{id}.tsv", id=ids),
        #expand(AMF+"{id}/amrfinder_resistance.pkl",id=ids),
    # ResFinder
        #expand(RES+"{state}/{id}/pheno_table_escherichia_coli.txt",
        #    id=ids,state=['raw','asmbl']), #raw, asmbl
        #expand(RES+"{id}/cgeresfinder_resistance.pkl", id=ids),
    # ABR resistance list
        #expand(ABR+"{id}/{db}_resistance.pkl", id=ids, db=abr_amr_dbs)
    # --------------------------------------------------------------------------
    # ABRicate presence/absence table
        #expand(ABR+"pres_abs_{abr_db}.tsv",abr_db=abr_dbs),#abr_db=['card','ncbi','resfinder']),
    # --------------------------------------------------------------------------
    ## Summaries
    # Summary by sample
        #expand(OUT+"summary_by_sample/{id}/{db}.tsv",id=ids,db=dbs)
    # Summary by sample for all
        #expand(OUT+"{db}_summary_by_sample_all.tsv", db=dbs),
    # Summary for each db
        #expand(OUT+"{db}_acc.tsv", db=dbs)
        #expand(OUT+"summary_tables/{db}_predictions.tsv",db=dbs)
        expand(OUT+"summary_tables/{db}_acc.tsv", db=['card','ncbi','resfinder'])
    # Combined summary file - rule not complete yet
        #expand(OUT+"comparison.tsv")
    #run:
    #    shell("rm {AMF_FLAG}")
        # remove the database flags so if it's rerun, it will update the dbs
        #shell("rm {general_out}*.done")

subworkflow download:
    workdir: "../"
    snakefile: "download.smk"

subworkflow assemble:
    workdir: "../"
    snakefile: "assemble.smk"

subworkflow filter_taxa:
    workdir: "../"
    snakefile: "filter_taxa.smk"

def get_asmbl_input(wcs):
    """
    Returns a list of input files to be used by Snakemake rules depending on
    the wgs_type passed (currently either standard or kraken).
    Path to the assemblies differ if using normal assemblies or kraken-filtered
    assemblies.

    :wgs_type: currently either standard or kraken; specified in config file
    :return: list of input files for use with snakefiles
    """
    # Try this with subworkflows?
    if wgs_type == 'standard':
        #return [assemble("data/genomes/asmbl/{id}/scaffolds.fasta")]
        return ["data/genomes/asmbl/{id}/scaffolds.fasta"]
    elif wgs_type == 'kraken':
        return [filter_taxa("data/genomes/kraken/{id}/filtered_ecoli.fa")]
    else:
        raise Exception("{} is an invalid wgs_type".format(wgs_type))

def get_raw_input_1(wcs):
    # if the id is in the master sheet, check if it's downloaded externally
    # otherwise pass to download subworkflow
    if wcs.id in master['id'].values.tolist():
        # look up file_id
        fid = master[master.id == wcs.id]['file_id'].values[0]
        # look up dataset
        ds = master[master.id == wcs.id]['dataset'].values[0]
        # If dataset is external
        if config['external_datasets'] and ds in config['external_datasets']:
            #path = config['external_datasets'][ds] + fid + '/'
            path = config['external_datasets'][ds]['raw_path'] + fid + '/'
            ext = config['external_datasets'][ds]['file_ext']
            # If the path doesn't exist, raise an error
            if not os.path.isdir(path):
                raise Exception("check path to external dataset {}".format(ds))
            # Walk through the path and look at the files
            for root, indirs, files in os.walk(path):
                # If BCRC or CIPARS then deal with the super fun naming scheme
                if ds == 'BCRC' or ds == 'CIPARS':
                    # Look at each file
                    for file in files:
                        # this one specific ID needs extra help because it sucks
                        if "0067K" in file:
                            if"_S2_" in file:
                                if "_R1_" in file:
                                    #print(path+file)
                                    return [path+file]
                        # if it's file 1
                        else:
                            if "_R1_" in file:
                                #print(path+file)
                                return [path+file]
                # If not BCRC or CIPARS assume normal paired end naming
                else:
                    # Look at each file
                    for file in files:
                        # if it's file 1
                        #if "_1.fastq.gz" in file:
                        if "_1{}".format(file_ext) in file:
                            #print(path+file)
                            return [path+file]
        # if its a downloaded genome, call the subworkflow
        else:
            #return [download("data/genomes/raw/{id}/{id}_1.fastq.gz")]
            return ["data/genomes/raw/{id}/{id}_1.fastq.gz"]
    else:
        #return [download("data/genomes/raw/{id}/{id}_1.fastq.gz")]
        return ["data/genomes/raw/{id}/{id}_1.fastq.gz"]

def get_raw_input_2(wcs):
    # if the id is in the master sheet, check if it's downloaded externally
    # otherwise pass to download subworkflow
    if wcs.id in master['id'].values.tolist():
        # look up file_id
        fid = master[master.id == wcs.id]['file_id'].values[0]
        # look up dataset
        ds = master[master.id == wcs.id]['dataset'].values[0]
        # If dataset is external
        if config['external_datasets'] and ds in config['external_datasets']:
            #path = config['external_datasets'][ds] + fid + '/'
            path = config['external_datasets'][ds]['raw_path'] + fid + '/'
            ext = config['external_datasets'][ds]['file_ext']
            if not os.path.isdir(path):
                raise Exception("check path to external dataset {}".format(ds))
            # Walk through the path
            for root, indirs, files in os.walk(path):
                # If BCRC or CIPARS then deal with the super fun naming scheme
                if ds == 'BCRC' or ds == 'CIPARS':
                    # Look at each file
                    for file in files:
                        # this one specific ID needs extra help because it sucks
                        if "0067K" in file:
                            if"_S2_" in file:
                                if "_R2_" in file:
                                    #print(path+file)
                                    return [path+file]
                        # if it's file 1
                        else:
                            if "_R2_" in file:
                                #print(path+file)
                                return [path+file]
                else:
                    # Look at each file
                    for file in files:
                        # if it's file 2
                        #if "_2.fastq.gz" in file:
                        if "_2{}".format(file_ext) in file:
                            return [path+file]
        # if its a downloaded genome, call the subworkflow
        else:
            #return [download("data/genomes/raw/{id}/{id}_1.fastq.gz")]
            return ["data/genomes/raw/{id}/{id}_1.fastq.gz"]
    # if its a downloaded genome, call the subworkflow
    else:
        #return [download("data/genomes/raw/{id}/{id}_2.fastq.gz")]
        return ["data/genomes/raw/{id}/{id}_2.fastq.gz"]

def get_raw_input(wcs):
    f1 = get_raw_input_1(wcs)
    f2 = get_raw_input_2(wcs)
    f3 = f1+f2
    return f3

################################################################################
## ABRicate - identify known antimicrobial reistance & virulence genes
################################################################################

# Update the abricate databases
rule update_abricate:
    group: 'update_abricate'
    output:
        touch(OUT+"update_{db}.done")
    shell:
        "abricate-get_db --db {wildcards.db} --force"

# Run ABRicate on assembly, compare with selected databases
rule abricate:
    group: 'abricate'
    input:
        #infile = get_asmbl_input,
        infile = "data/genomes/asmbl/{id}/scaffolds.fasta",
        update_flag = OUT+"update_{db}.done"
    output:
        ABR+"{id}/{db}.tsv",
    shell:
        "abricate --db {wildcards.db} {input.infile} > {output}"

# Error may be due to character limit.
# Worked when path was ABR = "results/wgs_{}/abricate/".format(wgs_type)
# works if using * SO
# THIS MAKES A FILE FOR ALL RESULTS NOT JUST IDS SPECIFIED.
# MUST FILTER IT AFTERWARDS IF ONLY WANT CERTAIN IDS.
# Create gene pres/abs file for the group
rule abricate_summary:
    group: 'abricate_summary'
    input:
        expand(ABR+"{id}/{{abr_db}}.tsv", id=ids)
    output:
        ABR+"pres_abs_{abr_db}.tsv"
    params:
        in_all = ABR+"*/{abr_db}.tsv"
    shell:
        "abricate --summary {params.in_all} > {output}"
        #"abricate --summary {input} > {output}"

# List of resistances for each id
rule abricate_resistance_list:
    group: 'abr_res_list'
    input:
        infile = ABR+"{id}/{db}.tsv"
    output:
        out = ABR+"{id}/{db}_resistance.pkl"
    run:
        # load the input
        df = pd.read_csv(input.infile,sep='\t')
        # Get the resistance column
        res_col = df['RESISTANCE'].tolist()
        # abricate output is already lowercase
        #res_col = [ i.lower() for i in res_col]
        # Remove duplicates
        res_col = list(set(res_col))
        # Some entries in the list have multiple items listed with /
        # Split them into two entries
        new_dc = []
        for item in res_col:
            # item must be a string; nans are not included
            if isinstance(item, str):
                if ';' in item:
                    # if the item is a list, split into separate items
                    new_dc_items = item.split(';')[:]
                    # strip whitespace
                    new_dc_items = [ i.strip() for i in new_dc_items ]
                    # append to new list
                    new_dc = new_dc + new_dc_items
                else:
                    new_dc = new_dc + [item]
        # Remove duplicates
        new_dc = list(set(new_dc))
        # ensure lowercase
        new_dc = [ i.lower() for i in new_dc ]
        # Save the list of drug classes the id is resistant to
        with open(output.out, 'wb') as f:
            pkl.dump(new_dc,f)


'''
# Condense the abricate output for each database into a single file
rule abricaknit:
    input:
        # {{}} expands only on the one given id instead of all ids
        # needed to prevent speed problems with recursion
        expand(ABR+"{{id}}_raw_out/{abr_db}.tsv", abr_db=abr_dbs)
    output:
        ABR+"{id}_abr.tsv"
    run:
        shell("python src/abricaknit.py -i {wildcards.id} -o {output}")
'''

''' # Make a gene-presence-absence file
rule abr_pres_abs:
    input:
        expand(ABR+"{id}/condensed.tsv", id=ids)
    output:
        ABR+DATASET+"summary.tsv"
    shell:
        "abricate --summary {input} > {output}"
'''


###############################################################################
# RGI / CARD
###############################################################################
"""
Run the Resistance Gene Identifier (RGI) on an assembly.
Also find the set of drug classes that the id is resistant to and save it.
https://card.mcmaster.ca/home
https://card.mcmaster.ca/analyze/
conda install -c bioconda rgi=5.1.0
"""

rule rgi:
    group: 'rgi'
    input:
        asmbl = get_asmbl_input
    output:
        RGI+"{id}/rgi.tsv",
    params:
        dir = RGI+ "{id}/",
        intermediate = RGI+"{id}/rgi.txt"
    threads:
        4
    run:
        shell("rgi main \
        --num_threads {threads} \
        --input_type contig \
        --alignment_tool BLAST \
        --data wgs \
        --clean \
        --input_sequence {input.asmbl} \
        --output_file {params.dir}rgi")
        shell("mv {params.intermediate} {output}")
# --clean removes temp files

# Save the list of drug classes that the id is resistant to
rule rgi_resistant_list:
    group: 'rgi_res_list'
    input:
        rgi = RGI+"{id}/rgi.tsv",
    output:
        rlist = RGI+"{id}/rgi_resistance.pkl",
    run:
        # load the input
        rgi = pd.read_csv(input.rgi,sep='\t')

        # mepA gene is generic MDR, and has Drug Class nan
        # dropping these nans for now
        rgi = rgi.dropna(subset=['Drug Class'])
        '''
        find_nan = rgi[rgi['Drug Class'].isnull()]
        print("********************")
        print(find_nan)
        print("********************")
        print(find_nan['Drug Class'].values[0])
        print(find_nan['Resistance Mechanism'].values[0])
        print(find_nan['AMR Gene Family'].values[0])
        print(find_nan['Best_Hit_ARO'].values[0])
        print("********************")
        sys.exit()
        '''

        # We care about 'Drug Class' column at the moment
        dc = rgi['Drug Class'].tolist()
        # Remove immediate duplicates
        dc = list(set(dc))
        # Some entries in the list have multiple items listed with ;
        # Split them into two entries
        new_dc = []
        for item in dc:
            if ';' in item:
                # if the item is a list, split into separate items
                new_dc_items = item.split(';')[:]
                # strip whitespace
                new_dc_items = [ i.strip() for i in new_dc_items ]
                # append to new list
                new_dc = new_dc + new_dc_items
            else:
                new_dc = new_dc + [item]
        # Remove the word ' antibiotic' it's redundant
        new_dc = [ i.replace(' antibiotic', '') for i in new_dc ]
        # Remove duplicates
        new_dc = list(set(new_dc))
        # ensure lowercase
        new_dc = [ i.lower() for i in new_dc ]
        # Save the list of drug classes the id is resistant to
        with open(output.rlist, 'wb') as f:
            pkl.dump(new_dc,f)


###############################################################################
# AMRFinder
###############################################################################
"""
Run AMRFinder on an assembly.
Also find the set of drug classes that the id is resistant to and save it.
https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
conda install -c bioconda ncbi-amrfinderplus
"""

rule update_amf:
    group: 'update_amf'
    output:
        touch(AMF_FLAG)
    shell:
        "amrfinder --update"

rule amf:
    group: 'amf'
    input:
        asmbl = get_asmbl_input,
        amrflag = AMF_FLAG
    output:
        AMF+"{id}/amrfinder.tsv"
    shell:
        "amrfinder \
        --nucleotide {input.asmbl} \
        --organism Escherichia \
        > {output}"


rule amf_plus:
    group: 'amf'
    input:
        asmbl = get_asmbl_input,
        amrflag = AMF_FLAG
    output:
        AMF+"{id}/amrfinderPLUS.tsv"
    shell:
        "amrfinder \
        --nucleotide {input.asmbl} \
        --organism Escherichia \
        --plus \
        > {output}"

# Save the list of drug classes that the id is resistant to
rule amf_resistant_list:
    group: 'amf_res_list'
    input:
        amf = AMF+"{id}/amrfinder.tsv"
    output:
        rlist = AMF+"{id}/amrfinder_resistance.pkl"
    run:
        # load the input
        amf = pd.read_csv(input.amf,sep='\t')
        # We care about 'Subclass' column at the moment
        # Don't need 'Class', because 'Subclass' is either more specific or is
        # the same as 'Class'
        dc = amf['Subclass'].tolist()
        # Make all lowercase
        dc = [ i.lower() for i in dc]
        # Remove duplicates
        dc = list(set(dc))
        # Some entries in the list have multiple items listed with /
        # Split them into two entries
        new_dc = []
        for item in dc:
            if '/' in item:
                # if the item is a list, split into separate items
                new_dc_items = item.split('/')[:]
                # strip whitespace
                new_dc_items = [ i.strip() for i in new_dc_items ]
                # append to new list
                new_dc = new_dc + new_dc_items
            else:
                new_dc = new_dc + [item]
        # Remove duplicates
        new_dc = list(set(new_dc))
        # ensure lowercase
        new_dc = [ i.lower() for i in new_dc ]
        # Save the list of drug classes the id is resistant to
        with open(output.rlist, 'wb') as f:
            pkl.dump(new_dc,f)


###############################################################################
# ResFinder
###############################################################################
"""
Full documentation:
https://bitbucket.org/genomicepidemiology/resfinder/src/master/

Using ResFinder 4.0, but had to make minor modification.
4.0 I had an error for some genomes, so I grabbed the master instead of 4.0,
then reverted a small piece of code to the 4.0 version because of a different
error.
"""
'''
# If resfinder not present, get it
rule resfinder_clone:
    output:
        resfinder_path+"README.md"
    shell:
        "git clone \
        -b 4.0 \
        https://git@bitbucket.org/genomicepidemiology/resfinder.git \
        {resfinder_path}"
'''


# Run ResFinder on each sample
rule resfinder:
    group: 'resf'
    input:
        asmbl = get_asmbl_input,
        raw = get_raw_input
    output:
        o1 = RES+"{state}/{id}/pheno_table.txt",
        o2 = RES+"{state}/{id}/pheno_table_escherichia_coli.txt"
    params:
        outdir = RES+"{state}/{id}/"
    run:
        print(input.raw)
        if wildcards.state=='raw':
            shell("{resfinder} \
            --inputfastq {input.raw} \
            --outputPath {params.outdir} \
            -s 'Escherichia coli' \
            --acquired \
            --point")
        elif wildcards.state=='asmbl':
            shell("{resfinder} \
            --inputfasta {input.asmbl} \
            --outputPath {params.outdir} \
            -s 'Escherichia coli' \
            --acquired \
            --point")
        else:
            raise Exception("Wrong state {}".format(wildcards.state))

# the resf output phenotype tables are text files with a table in the middle
rule resf_extract_table:
    group: 'resf_extract_table'
    input:
        i1 = RES+"{state}/{id}/pheno_table.txt",
        i2 = RES+"{state}/{id}/pheno_table_escherichia_coli.txt"
    output:
        o1 = RES+"{state}/{id}/pheno_table.tsv",
        o2 = RES+"{state}/{id}/pheno_table_escherichia_coli.tsv"
    run:
        shell("sed -e '/^$/,/^$/!d' -e '/^$/d' -e 's/^# *//' {input.i1} > {output.o1}")
        shell("sed -e '/^$/,/^$/!d' -e '/^$/d' -e 's/^# *//' {input.i2} > {output.o2}")

# The e-coli only output doesnt have co-amoxiclav anymore?
# So use the regular pheno output file
rule resf_resistant_list:
    group: 'resf_res_list'
    input:
        res = RES+"{state}/{id}/pheno_table.tsv",
        res_ec = RES+"{state}/{id}/pheno_table_escherichia_coli.tsv",
    output:
        rlist = RES+"{state}/{id}/cgeresfinder_resistance.pkl",
        rlist_ec = RES+"{state}/{id}/cgeresfinder_resistance_ec.pkl",
    run:
        # load the input phenotype file
        res = pd.read_csv(input.res,sep='\t')
        # Keep all rows with "Resistant" and drop all rows with "No resistance"
        # in the "WGS-predicted phenotype" column
        #********
        # Note: considering all matches, level 1-3
        #********
        res = res.loc[res['WGS-predicted phenotype'] == "Resistant"]
        # turn the drug names into a list and pickle it
        reslist = res['Antimicrobial'].tolist()
        # Save the list of drug classes the id is resistant to
        with open(output.rlist, 'wb') as f:
            pkl.dump(reslist,f)
        ########################################################################
        # E coli version
        # load the input phenotype file
        res_ec = pd.read_csv(input.res_ec,sep='\t')
        # Keep all rows with "Resistant" and drop all rows with "No resistance"
        # in the "WGS-predicted phenotype" column
        #********
        # Note: considering all matches, level 1-3
        #********
        res_ec = res_ec.loc[res_ec['WGS-predicted phenotype'] == "Resistant"]
        # turn the drug names into a list and pickle it
        reslist_ec = res_ec['Antimicrobial'].tolist()
        # Save the list of drug classes the id is resistant to
        with open(output.rlist_ec, 'wb') as f_ec:
            pkl.dump(reslist_ec,f_ec)


'''
# Run ResFinder on each sample
rule resfinder:
    group: 'resf'
    input:
        asmbl = get_asmbl_input,
    output:
        o1 = RES+"{id}/pheno_table.txt",
        o2 = RES+"{id}/pheno_table_escherichia_coli.txt"
    params:
        outdir = RES+"{id}"
    shell:
        "{resfinder} \
        --inputfasta {input.asmbl} \
        --outputPath {params.outdir} \
        --db_path_res {resfinder_db_path} \
        -s 'Escherichia coli' \
        --acquired \
        --point"
'''
'''
# the resf output phenotype tables are text files with a table in the middle
rule resf_extract_table:
    input:
        i1 = RES+"{id}/pheno_table.txt",
        i2 = RES+"{id}/pheno_table_escherichia_coli.txt"
    output:
        o1 = RES+"{id}/pheno_table.tsv",
        o2 = RES+"{id}/pheno_table_escherichia_coli.tsv"
    run:
        shell("sed -e '/^$/,/^$/!d' -e '/^$/d' -e 's/^# *//' {input.i1} > {output.o1}")
        shell("sed -e '/^$/,/^$/!d' -e '/^$/d' -e 's/^# *//' {input.i2} > {output.o2}")
'''

'''
rule resf_resistant_list:
    group: 'resf_res_list'
    input:
        res = RES+"{id}/pheno_table_escherichia_coli.tsv"
    output:
        rlist = RES+"{id}/cgeresfinder_resistance.pkl"
    run:
        # load the input phenotype file
        res = pd.read_csv(input.res,sep='\t')
        # Keep all rows with "Resistant" and drop all rows with "No resistance"
        # in the "WGS-predicted phenotype" column
        #********
        # Note: considering all matches, level 1-3
        #********
        res = res.loc[res['WGS-predicted phenotype'] == "Resistant"]
        # turn the drug names into a list and pickle it
        reslist = res['Antimicrobial'].tolist()
        # Save the list of drug classes the id is resistant to
        with open(output.rlist, 'wb') as f:
            pkl.dump(reslist,f)
'''

''' ResF 3.2
shell("resfinder.py \
--inputfile {input.asmbl} \
--outputPath {params.outdir} \
--databasePath {resfinder_db_path} \
--quiet")
shell("mv {params.intermediate} {output}")
'''
'''
ResFinder 3.2 output a json file; 4 noe has tsv
#find all keys in a json file
def find_values(id, json_repr):
    results = []
    def _decode_dict(a_dict):
        try:
            results.append(a_dict[id])
        except KeyError:
            pass
        return a_dict
    json.loads(json_repr, object_hook=_decode_dict) # Return value ignored.
    return results

rule resf_resistant_list:
    group: 'resf_res_list'
    input:
        j_file = RES+"{id}/cgeresfinder.json"
    output:
        rlist = RES+"{id}/cgeresfinder_resistance.pkl"
    run:
        import json

        with open(input.j_file, 'rb') as jf:
            # open the json file
            j = jf.read()
            # find all phenotype predictions
            r = find_values('predicted_phenotype',j)
            # enforce lowercase
            r = [i.lower() for i in r]
            # remove warning: uncurated genes
            r = [ i for i in r if 'warning:' not in i ]
            # list is in the format ['sulfonamide resistance' etc]
            # split each entry to just get the class
            r = [i.split(' resistance')[0] for i in r]
            # remove duplicates
            r = list(set(r))
            # Save the list of drug classes the id is resistant to
            with open(output.rlist, 'wb') as f:
                pkl.dump(r,f)
'''

###############################################################################
# Summary Tables
###############################################################################
def get_all_resistance_lists(wcs):
    if wcs.db == 'rgi':
        return expand(RGI+"{id}/rgi_resistance.pkl",id=ids)
    elif wcs.db == 'amrfinder':
        return expand(AMF+"{id}/amrfinder_resistance.pkl",id=ids)
    elif wcs.db == 'asmblcgeresfinder':
        return expand(RES+"asmbl/{id}/cgeresfinder_resistance.pkl",id=ids)
    elif wcs.db == 'rawcgeresfinder':
        return expand(RES+"raw/{id}/cgeresfinder_resistance.pkl",id=ids)
    elif wcs.db in ['card', 'ncbi', 'resfinder']:
        return expand(ABR+"{id}/{db}_resistance.pkl",id=ids,db=wcs.db)
    else:
        raise Exception("Wrong db")


rule fast_table:
    group: 'fast_table'
    input:
        get_all_resistance_lists
    output:
        out = OUT+"summary_tables/{db}_predictions.tsv"
    run:
        # Database name
        db = wildcards.db

        # Load the metadata master
        master_df = pd.read_csv(config['metadata_master'],sep='\t')

        # Get list of only drugs we made models for
        drugs = config['drugs_we_care_about']
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
                amr_lab = metadata_row['SIR_'+drug].values[0]
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

                # SXT is a special case
                if drug == 'SXT':
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


rule fast_table_compare:
    group:' fast_table_compare'
    input:
        pred = OUT+"summary_tables/{db}_predictions.tsv"
    output:
        out = OUT+"summary_tables/{db}_acc.tsv"
    run:
        # Load the file
        all_df = pd.read_csv(input.pred,sep='\t')
        print(all_df)

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
        drugs = config['drugs_we_care_about']
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
            n_lab_R = df[lab_col].value_counts()['R']

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



'''
def get_resistance_list(wcs):
    if wcs.db == 'rgi':
        return RGI+"{id}/rgi_resistance.pkl"
    elif wcs.db == 'amf':
        return AMF+"{id}/amrfinder_resistance.pkl"
    elif wcs.db == 'cgeresfinder':
        return RES+"{id}/cgeresfinder_resistance.pkl"
    elif wcs.db in ['card', 'ncbi', 'resfinder']:
        return ABR+"{id}/{db}_resistance.pkl"
    else:
        raise Exception("Wrong db")
'''
'''
def get_resistance_list(wcs):
    if wcs.db == 'rgi':
        return RGI+"{id}/rgi_resistance.pkl"
    elif wcs.db == 'amf':
        return AMF+"{id}/amrfinder_resistance.pkl"
    elif wcs.db == 'asmblcgeresfinder':
        return RES+"asmbl/{id}/cgeresfinder_resistance.pkl"
    elif wcs.db == 'rawcgeresfinder':
        return RES+"raw/{id}/cgeresfinder_resistance.pkl"
    elif wcs.db in ['card', 'ncbi', 'resfinder']:
        return ABR+"{id}/{db}_resistance.pkl"
    else:
        raise Exception("Wrong db")
'''
'''
# make a per-drug summary for each sample
rule db_summary_by_sample:
    group: 'db_summary_by_samp'
    input:
        res_list = get_resistance_list
    output:
        summary = OUT+"summary_by_sample/{id}/{db}.tsv"
    run:
        # Genome's id
        id = wildcards.id
        db = wildcards.db

        # Get list of only drugs we made models for
        drugs = config['drugs_we_care_about']

        # Load drug information
        df = pd.read_csv(config['antibiotic_master'],sep='\t')
        # Filter for only the columns we want
        df = df[['code', 'name', 'name_alt', 'group_broad', 'group', 'class', 'subclass', 'card_class', 'combination']]
        # Filter for the drugs in the list
        df = df[df['code'].isin(drugs)]
        df = df.reset_index(drop=True)
        # Make names lowercase
        df['name'] = df['name'].astype(str).apply(lambda x: x.lower())
        df['name_alt'] = df['name_alt'].astype(str).apply(lambda x: x.lower())
        # Make the drugs list match the order of the df
        drugs = df['code'].tolist()

        # Get the id's lab AMR stats from the master df
        id_row = master.loc[master['id'] == id]
        # Make a list of the AMR stats in the same order the names appear
        lab_amr = []
        for drug in drugs:
            sir = id_row['SIR_'+drug].values[0]
            lab_amr = lab_amr + [sir]
        # Add the Lab AMR to the df
        df['lab_amr'] = lab_amr

        # Load the list of resistances
        with open(input.res_list, 'rb') as r:
            res_list = pkl.load(r)

        # Init the new columns
        db_amr = []
        db_class = []

        # Create a row for each drug
        for drug in drugs:
            # get the row of the df for this drug
            drug_row = df.loc[df['code'] == drug]
            # Is the drug a combination drug eg co-amoxiclav
            combo_flag = drug_row['combination'].values[0]
            # Create a list of the drug's names, ensure no nans
            ori_drug_name = drug_row['name'].values.tolist()[0]
            drug_names = drug_row['name'].values.tolist()
            # the alts cell has items delimited by commas but is not actually in
            # list form yet
            alts = drug_row['name_alt'].values.tolist()[0].split(',')
            # strip whitepaces
            alts = [a.strip(' ') for a in alts]
            # make one big list of drug names
            drug_names.extend(alts)

            # no nans!
            drug_names = list(filter(lambda a: a != 'nan', drug_names))
            if np.nan in drug_names or 'nan' in drug_names or 'NaN' in drug_names:
                raise Exception("nans!")

            #print("drug names",drug_names)


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

            # SXT is a special case
            if drug == 'SXT':
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
                tmp_amr = np.nan
                tmp_class = np.nan
                if db == 'cgeresfinder':
                    if 'sulphonamide' in res_list and 'trimethoprim' in res_list:
                        tmp_amr = 'R'
                        tmp_class = 'sulphonamide, trimethoprim'
                elif db == 'resfinder':
                    if 'sulfamethoxazole' in res_list and 'trimethoprim' in res_list:
                        tmp_amr = 'R'
                        tmp_class = 'sulfamethoxazole, trimethoprim'
                elif db == 'ncbi' or db == 'amf':
                    if 'sulfonamide' in res_list and 'trimethoprim' in res_list:
                        tmp_amr = 'R'
                        tmp_class = 'sulfonamide, trimethoprim'
                elif db == 'card' or db == 'rgi':
                    if 'sulfonamide' in res_list and 'diaminopyrimidine' in res_list:
                        tmp_amr = 'R'
                        tmp_class = 'sulfonamide, diaminopyrimidine'
                db_amr += [tmp_amr]
                db_class += [tmp_class]

            # Special cases where we only want to look at direct drug
            # for resfinder
            # dont do this - worse results
            # also, dont do for AMP, CIP for ncbi (doesnt work)
            #elif db in ['amf','ncbi'] and drug in ['AMP','CIP','GEN','TBM']:
            #    if any(d in drug_names for d in res_list):
            #        db_amr += ['R']
            #        # same drug name as it appears in res_list
            #        #db_class += [i for i in drug_names if i in res_list]
            #        # what names appeared in both lists
            #        overlap = list(set(drug_names).intersection(res_list))
            #        db_class += overlap
            #    else:
            #        db_amr += [np.nan]
            #        db_class += [np.nan]

            # Special cases where we only want to look at cephalosporin
            # but not beta-lactam (large accuracy boost)
            elif db in ['amf','ncbi'] and drug in ['CET','CRO','CTX','CTZ','CXM','TIO']:
                if any(d in drug_names for d in res_list):
                    db_amr += ['R']
                    # same drug name as it appears in res_list
                    #db_class += [i for i in drug_names if i in res_list]
                    # what names appeared in both lists
                    overlap = list(set(drug_names).intersection(res_list))
                    db_class += [','.join(overlap)]
                # if the drug class is in the res list
                # class is a list for [sulfonamide, sulphonamide] because
                # multiple spellings
                #elif drug_class in res_list:
                elif any(d in drug_class for d in res_list):
                    db_amr += ['R']
                    #db_class += [drug_class]
                    overlap = list(set(drug_class).intersection(res_list))
                    db_class += [','.join(overlap)]
                else:
                    db_amr += [np.nan]
                    db_class += [np.nan]

            # rest of the drugs proceed the same
            else:
                # if one of the drug names is in the resistant list
                if any(d in drug_names for d in res_list):
                    db_amr += ['R']
                    # same drug name as it appears in res_list
                    #db_class += [i for i in drug_names if i in res_list]
                    # what names appeared in both lists
                    overlap = list(set(drug_names).intersection(res_list))
                    db_class += [','.join(overlap)]
                # if the drug class is in the res list
                # class is a list for [sulfonamide, sulphonamide] because
                # multiple spellings
                #elif drug_class in res_list:
                elif any(d in drug_class for d in res_list):
                    db_amr += ['R']
                    #db_class += [drug_class]
                    overlap = list(set(drug_class).intersection(res_list))
                    db_class += [','.join(overlap)]
                # if the card class is in the res list
                # CARD considers nalidixic acid fluoroquinolone while others dont
                elif (db in ['card','rgi']) and card_class in res_list:
                    db_amr += ['R']
                    db_class += [card_class]
                # group eg cephem
                elif group in res_list:
                    db_amr += ['R']
                    db_class += [group]
                # group_broad group eg beta-lactam
                elif group_broad in res_list:
                    db_amr += ['R']
                    db_class += [group_broad]
                else:
                    db_amr += [np.nan]
                    db_class += [np.nan]

        #turn res list into string
        if len(res_list)>1:
            res_list = ', '.join(res_list)
        elif len(res_list)==0:
            res_list = np.nan
        else:
            res_list = res_list[0]


        # Add the AMR to the df
        df[db] = db_amr
        df[db+'_class'] = db_class
        df[db+'_full_res_list'] = [res_list]*len(db_amr)

        # Add id column
        df['id'] = [id]*len(drugs)


        df.to_csv(output.summary, sep='\t', index=False)
'''
'''
# cat all of the sample summaries for each db
rule db_summary_combined:
    group: 'db_summary_combined'
    input:
        expand(OUT+"summary_by_sample/{id}/{{db}}.tsv", id=ids)
    output:
        out = OUT+"cat_summary_by_sample/cat_{db}.tsv"
    run:
        cat = pd.concat([pd.read_csv(f, sep='\t') for f in input])
        #print(cat)
        cat.to_csv(output.out, index=False, sep='\t')
'''
'''
# create a condensed summary of accuracy etc for each db
rule db_summary_accuracy:
    group: 'db_summary_accuracy'
    input:
        infile = OUT+"cat_summary_by_sample/cat_{db}.tsv"
    output:
        out = OUT+"{db}_acc.tsv"
    run:
        db = wildcards.db

        all_df = pd.read_csv(input.infile, sep='\t')

        drugs = config['drugs_we_care_about'] # drugs to predict on

        # Load the antibiotic drug master file to get subclass etc info
        dmdf = pd.read_csv(config['antibiotic_master'],sep='\t')
        # Ensure all drug names and classes are lowercase
        dmdf['name'] = dmdf['name'].astype(str).apply(lambda x: x.lower())
        dmdf['name_alt'] = dmdf['name_alt'].astype(str).apply(lambda x: x.lower())
        dmdf['group_broad'] = dmdf['group_broad'].astype(str).apply(lambda x: x.lower())
        dmdf['group'] = dmdf['group'].astype(str).apply(lambda x: x.lower())
        dmdf['class'] = dmdf['class'].astype(str).apply(lambda x: x.lower())
        dmdf['subclass'] = dmdf['subclass'].astype(str).apply(lambda x: x.lower())

        # Initialize a list to hold one dict for each drug; to be made into a df
        row_list = []

        for drug in drugs:

            # Get the drug info from the drug master
            drug_name = dmdf.loc[dmdf['code'] == drug,'name'].item()
            drug_name_alt = dmdf.loc[dmdf['code'] == drug,'name_alt'].item()
            drug_class = dmdf.loc[dmdf['code'] == drug,'class'].item()
            drug_subclass = dmdf.loc[dmdf['code'] == drug,'subclass'].item()

            # Filter dataframe for current drug
            df = all_df[all_df['code']==drug]

            # Filter to keep only the samples which have
            # experimental/lab values for SIR
            df = df[df['lab_amr'].notnull()]

            # When RGI or AMRFinder doesn't find R / resistance, then it is S
            df[[db]] = df[[db]].fillna(value='S')

            # Number of samples with lab values for this drug
            n_samples = len(df['id'].tolist())

            # Number of samples that have a lab value of Resistant
            n_lab_R = df.lab_amr.value_counts()['R']

            # Make a column to indicate if the result matched the lab amr
            df['correct'] = np.where(df['lab_amr']==df[db], 1, 0)
            # Count the correct predictions by summing the column
            correct_count = df['correct'].sum()

            # Make a col to indicate if the result was R
            df['predicted_R'] = np.where(df[db]=='R', 1, 0)
            # Count the number of R predictions
            # In the case of the penicillins, they're all R
            predicted_R_count = df['predicted_R'].sum()

            # What did the db classify the drug as? eg is db using penicillin,
            # penam, or beta-lactam to make the choice of R?
            db_class = df[db+'_class'].tolist()
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
                'class': drug_class,
                'subclass': drug_subclass,
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
        #print(new_df)

        new_df.to_csv(output.out, index=False, sep='\t')
'''
'''
rule full_comparison_table:
    group: 'full_compare_table'
    input:
        oth_tools = expand(OUT+"{db}_acc.tsv",db=dbs),
        xgb = "results/wgs_standard/11mer/models/all/best_xgb_nested_cv_results.tsv",
        svm = "results/wgs_standard/11mer/models/all/best_svm_nested_cv_results.tsv"
    output:
        out = OUT+"comparison.tsv"
    run:
        # Load list of all drug acronyms
        drugs = config['drugs_we_care_about'] # drugs to predict on

        # Load the model accuracy, just need the drug and acc cols
        df = pd.read_csv(input.xgb, sep='\t', usecols=['predfor','mean'])
        # Strip the "SIR" from the front of the predfor column for
        # easier indexing
        df['predfor'] = df['predfor'].str.replace(r'SIR', '')
        # rename cols
        df.rename(columns={'predfor':'drug','mean':'XGBoost'})

        print(df)

        # Update the dataframe with the result of all of the other amr tools
        for file in input.oth_tools:
            # db name
            #db_name = file.split()
            # load the file
            df = pd.read_csv(file, sep='\t')
            # Want correct_perc, R_perc, and rename by db, and update main df


        # Load each input


'''
