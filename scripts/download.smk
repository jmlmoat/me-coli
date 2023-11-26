"""
Download data
"""

import pandas as pd

configfile: "config/config.yaml"


#-------------------------------------------------------------------------------
# Load files
# Load the metadata master
master = pd.read_csv(config['metadata_master'],sep='\t',low_memory=False)

# Load the reference genome file
ref_master = pd.read_csv(config['reference_master'],sep='\t')

# Metagenome metadata file
meta_master = pd.read_csv(config['metagenome_master'],sep='\t')


#-------------------------------------------------------------------------------
# Paths - output
RAW = "data/genomes/raw/" # raw genomes
REF_ASMBL = "data/reference/" # reference genome assemblies
ASMBL = "data/genomes/asmbl/" # assembled genomes


#-------------------------------------------------------------------------------
# WGS - raw

# ENA genomes to be downloaded
ena_ids = master[master['collection'].isin(['ENA'])]['id'].tolist()

# NCBI/PATRIC isolates that can be downloaded via SRA accession
#ncbi_downloadable_df = master.loc[(master['collection'] != 'CIPARS') & (master['collection'] != 'BCRC') & (master['collection'] != 'ENA')]
srr_df = master.loc[ (master['download_key'] == 'srr') ]
srr_ids = srr_df['id'].tolist()

# ResFinder 4.0 Paper Validation Data
# To download from the resfinder project: annoyingly the folders and files
# and ids annoyingly have different names. Not like the other ENA download.
## PRJNA6164
PRJNA6164 = pd.read_csv(config['PRJNA6164'], sep='\t')
PRJNA6164_ids = PRJNA6164['id'].tolist()
# dict for renaming appropriately
PRJNA6164_dict = dict(zip(PRJNA6164.id, PRJNA6164.experiment_accession))
## PRJEB22091
PRJEB22091 = pd.read_csv(config['PRJEB22091'], sep='\t')
PRJEB22091_ids = PRJEB22091['id'].tolist()
# dict for renaming appropriately
PRJEB22091_dict = dict(zip(PRJEB22091.id, PRJEB22091.experiment_accession))


#-------------------------------------------------------------------------------
# WGS - assemblies

# The download_key for assemblies is either for refseq or genbank
asmbl_df = master.loc[ (master['download_key'] == 'asmbl_refseq') | (master['download_key'] == 'asmbl_genbank') ]
asmbl_ids = asmbl_df['id'].tolist()

'''
refseq_df = master.loc[ (master['download_key'] == 'asmbl_refseq') ]
refseq_ids = refseq_df['id'].tolist()
# filter for the ones with GenBank
genbank_df = master.loc[ (master['download_key'] == 'asmbl_genbank') ]
genbank_ids = genbank_df['id'].tolist()
'''


#-------------------------------------------------------------------------------
# Reference samples

# List of reference genomes to download
#ref_ids = ref_master['id'].tolist()
# filter for REF-M only
ref_meta_only = ref_master[ (ref_master['dataset'] == 'REF-M') ]
ref_meta_ids = ref_meta_only['biosample'].tolist()
#print(refm_ids)


# Reference samples for tree
ref_tree_only = ref_master[ (ref_master['dataset'] == 'REF-TREE') ]
ref_tree_ids = ref_tree_only['biosample'].tolist()


#-------------------------------------------------------------------------------
# Metagenome stuff

# List of metagenomic data to download
m_ids = meta_master['id'].tolist()
# Metagenomic samples from Rahat's data; N = 4 + 13 + 6 + 28 = 51


#-------------------------------------------------------------------------------
# When running on a cluster, these need to be local so they don't flood and fail
# Could also try messing with --dependency singleton (but can run 6 at time fine
# so this would be a slower method) or profiles
#localrules: ENA_download, ref_ncbi_fna_list, ref_ncbi_dl
localrules: raw_wgs_dl, ref_ncbi_fna_list, ref_ncbi_dl



rule all:
    input:
    # WGS raw ------------------------------------------------------------------
    # ENA raw
        expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=ena_ids, pe=[1,2]),
    # NCBI raw
        expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=srr_ids, pe=[1,2]),
    # ResFinder Validation set (ENA) raw
        expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=PRJNA6164_ids, pe=[1,2]),
        expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=PRJEB22091_ids, pe=[1,2]),
    # WGS assemblies -----------------------------------------------------------
    # NCBI assemblies
        expand(ASMBL+"{id}/scaffolds.fasta", id=asmbl_ids),
    # Reference Genome Assemblies
        expand(REF_ASMBL+"{ref_id}/{ref_id}.fasta", ref_id=ref_meta_ids),
        expand(REF_ASMBL+"{ref_id}/{ref_id}.fasta", ref_id=ref_tree_ids),

     # WGS - ENA raw -----------------------------------------------------------
     # ENA - download raw genomes; pe is paried end
        #expand(RAW+"{ena_id}/{ena_id}_{pe}.fastq.gz", ena_id=ena_ids, pe=[1,2]),
     # ENA - download the resfinder comparison dataset
        #expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=PRJNA6164_ids, pe=[1,2]),
        #expand(RAW+"{id}/{id}_{pe}.fastq.gz", id=PRJEB22091_ids, pe=[1,2]),
    # WGS - NCBI raw -----------------------------------------------------------
    # WGS - NCBI assemblies ----------------------------------------------------
    # WGS - NCBI reference -----------------------------------------------------
    # NCBI - Download reference genomes, rename as needed...
        ##expand(REF_ASMBL+"{ref_id}/{ref_id}.fna.gz", ref_id=ref_ids),
        #expand(REF_ASMBL+"{ref_id}/{ref_id}.fasta", ref_id=ref_meta_ids),
    # Meta - NCBI -------- -----------------------------------------------------
    # Meta (NCBI) - Download raw metagenomic data
        #expand("data/metagenomes/raw/{m_id}/{m_id}_pass_1.fastq.gz", m_id=m_ids),



################################################################################
""" replaced by raw_wgs_dl
Download raw WGS from the European Nucleotide Archive (ENA) using
enaBrowserTools
"""
'''
rule ENA_download:
    output:
        o1 = RAW+"{id}/{id}_1.fastq.gz",
        o2 = RAW+"{id}/{id}_2.fastq.gz"
    params:
        outdir = RAW
    # Force all download jobs to be submitted as one group when in cluster mode
    #group:
    #    'ena-dl'
    # This group will use maximum 6 cores, so the request limit isn't exceeded
    #threads:
    #    6
    #resources: cpus=1
    run:
        # resf ids require extra crap because it makes you use the experiment
        # accession, not sure why
        # Flag indicating if current id is resfinder 4 data or not
        resfinder4_id = False
        # if its a resfinder id get the download accession based on bioproject
        if wildcards.id in PRJNA6164_ids:
            resfinder4_id = True
            id = wildcards.id
            dl_accession = PRJNA6164_dict[id]
        if wildcards.id in PRJEB22091_ids:
            id = wildcards.id
            dl_accession = PRJEB22091_dict[id]

        # if its a resfinder 4 id, download it in this annoying way
        if resfinder4_id:
        #if wildcards.id in PRJNA6164_ids:
            #id = wildcards.id
            #dl_accession = PRJNA6164_dict[id]
            dl_folder = "{}{}/".format(params.outdir,dl_accession)
            interm_f1_gz = "{0}/{0}_1.fastq.gz".format(id)
            interm_f2_gz = "{0}/{0}_2.fastq.gz".format(id)
            # enaDataGet https://github.com/enasequence/enaBrowserTools
            shell("enaDataGet -w -e -f fastq -d {params.outdir} {dl_accession}")
            #print("Unzipping")
            #shell("echo {params.outdir} {dl_accession} {dl_folder}{interm_f1_gz}")
            #shell("gunzip {dl_folder}{interm_f1_gz}")
            #shell("gunzip {dl_folder}{interm_f2_gz}")
            shell("mv {dl_folder}{id} {params.outdir}/")
            shell("rm -rf {dl_folder}")
        # standard download
        else:
            # enaDataGet https://github.com/enasequence/enaBrowserTools
            shell("enaDataGet -w -e -f fastq -d {params.outdir} {wildcards.id}")
            #print("Unzipping")
            #shell("gunzip {output.o1}.gz")
            #shell("gunzip {output.o2}.gz")
'''



################################################################################
"""
Download raw WGS
- download from the European Nucleotide Archive (ENA) using enaBrowserTools
- download from the NCBI using fastq-dump
"""

rule raw_wgs_dl:
    output:
        o1 = RAW+"{id}/{id}_1.fastq.gz",
        o2 = RAW+"{id}/{id}_2.fastq.gz"
    params:
        i1 = RAW+"{id}/{id}_pass_1.fastq.gz",
        i2 = RAW+"{id}/{id}_pass_2.fastq.gz",
        outdir = RAW,
        collection = lambda wildcards: master.loc[master['id']==wildcards.id]['collection'].values[0]
    run:
        #--------------------------------------
        # Download procedure for ENA samples
        if params.collection == 'ENA':
            # resf ids require extra crap because it makes you use the experiment
            # accession, not sure why
            # Flag indicating if current id is resfinder 4 data or not
            resfinder4_id = False
            # if its a resfinder id get the download accession based on bioproject
            if wildcards.id in PRJNA6164_ids:
                resfinder4_id = True
                id = wildcards.id
                dl_accession = PRJNA6164_dict[id]
            if wildcards.id in PRJEB22091_ids:
                id = wildcards.id
                dl_accession = PRJEB22091_dict[id]

            # if its a resfinder 4 id, download it in this annoying way
            if resfinder4_id:
            #if wildcards.id in PRJNA6164_ids:
                #id = wildcards.id
                #dl_accession = PRJNA6164_dict[id]
                dl_folder = "{}{}/".format(params.outdir,dl_accession)
                interm_f1_gz = "{0}/{0}_1.fastq.gz".format(id)
                interm_f2_gz = "{0}/{0}_2.fastq.gz".format(id)
                # enaDataGet https://github.com/enasequence/enaBrowserTools
                shell("enaDataGet -w -e -f fastq -d {params.outdir} {dl_accession}")
                #print("Unzipping")
                #shell("echo {params.outdir} {dl_accession} {dl_folder}{interm_f1_gz}")
                #shell("gunzip {dl_folder}{interm_f1_gz}")
                #shell("gunzip {dl_folder}{interm_f2_gz}")
                shell("mv {dl_folder}{id} {params.outdir}/")
                shell("rm -rf {dl_folder}")
            # standard download
            else:
                # enaDataGet https://github.com/enasequence/enaBrowserTools
                shell("enaDataGet -w -e -f fastq -d {params.outdir} {wildcards.id}")
        #--------------------------------------
        # Download procedure for NCBI-available samples
        elif params.collection in ['NCBI','PATRIC','AAFC']:
            shell("fastq-dump \
            --outdir {params.outdir}{wildcards.id} \
            --gzip \
            --skip-technical \
            --read-filter pass \
            --dumpbase \
            --split-3 \
            --clip \
            {wildcards.id}")
            shell("mv {params.i1} {output.o1}")
            shell("mv {params.i2} {output.o2}")
#  --readids \
# do not use this because it messes up spades in assembly


################################################################################
"""
Download assembled WGS
- download from the NCBI using fastq-dump
"""

rule asmbl_wgs_dl:
    output:
        o = ASMBL+"{id}/scaffolds.fasta",
    params:
        intermediate1 = ASMBL+"{id}/scaffolds.fasta.gz",
        #intermediate2 = ASMBL+"{id}/{id}.fna",
        outdir = ASMBL,
        ref_or_gen = lambda wildcards: asmbl_df.loc[asmbl_df['id']==wildcards.id]['download_key'].values[0],
        refseq_ftppath = lambda wildcards: asmbl_df.loc[asmbl_df['id']==wildcards.id]['refseq_ftppath'].values[0],
        genbank_ftppath = lambda wildcards: asmbl_df.loc[asmbl_df['id']==wildcards.id]['genbank_ftppath'].values[0],
    run:
        ftp_path = ''
        if params.ref_or_gen == 'asmbl_refseq':
            ftp_path = params.refseq_ftppath
        elif params.ref_or_gen == 'asmbl_genbank':
            ftp_path = params.genbank_ftppath
        else:
            raise Exception("invalid assembly download key")
        shell("""
        curl -o {params.intermediate1} {ftp_path}
        gunzip {params.intermediate1}
        """)


#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/770/275/GCF_000770275.1_ASM77027v1/GCF_000770275.1_ASM77027v1_genomic.fna.gz
#wget -O SAMN03075588/SAMN03075588.fna.gz -i ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/770/275/GCF_000770275.1_ASM77027v1/GCF_000770275.1_ASM77027v1_genomic.fna.gz

# before using curl, tried this
#wget -O {params.intermediate1} -i {ftp_path}
#gunzip {params.intermediate1}
#mv {params.intermediate2} {output}







################################################################################
"""
Download reference genomes from the NCBI
Uses list of ids
"""


rule ref_ncbi_fna_list:
    output:
        temp(REF_ASMBL+"{id}_esearch.txt")
    run:
        shell("""
        esearch -db biosample -query {wildcards.id} \
          | elink -target assembly \
          | esummary \
          | grep "FtpPath_RefSeq" \
          | sed -r "s|.+>(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)<.+|\\1\\2/\\2_genomic.fna.gz|" \
          > {output}
        """)

rule ref_ncbi_dl:
    input:
        REF_ASMBL+"{id}_esearch.txt"
    output:
        REF_ASMBL+"{id}/{id}.fasta"
    params:
        intermediate1 = REF_ASMBL+"{id}/{id}.fna.gz",
        intermediate2 = REF_ASMBL+"{id}/{id}.fna"
    shell:
        """
        wget -O {params.intermediate1} -i {input}
        gunzip {params.intermediate1}
        mv {params.intermediate2} {output}
        """
###!!!!!!!!!!! swith to curl with following file instead
# I would switch it right now but i dont want to test it
'''
rule ref_ncbi_dl:
    input:
        REF_ASMBL+"{id}_esearch.txt"
    output:
        REF_ASMBL+"{id}/{id}.fasta"
    params:
        intermediate1 = REF_ASMBL+"{id}/{id}.fasta.gz",
        #intermediate2 = REF_ASMBL+"{id}/{id}.fna"
    shell:
        """
        curl -o {params.intermediate1} {input}
        gunzip {params.intermediate1}
        """
'''


### OLD VERSION THAT DOWNLOADED WHOLE BIOPROJECT AND TOOK 1 GENOME FROM IT
'''
# get input file list
rule ref_ncbi_fna_list:
    output:
        temp(REF_ASMBL+"{id}_esearch.txt")
    #group:
    #    'ncbi'
    #threads:
    #    6
    params:
        prjna = lambda wildcards: ref_master.loc[ref_master['id']==wildcards.id]['bioproject'].values[0]
    run:
        shell("""
        esearch -db bioproject -query {params.prjna} \
          | elink -target nucleotide \
          | esummary \
          | grep "FtpPath_RefSeq" \
          | sed -r "s|.+>(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)<.+|\\1\\2/\\2_genomic.fna.gz|" \
          > {output}
        """)
        # If the bioproject has multiple components, we only want one, so need
        # to trim off the excess
        # Open the file to read it
        with open(output[0], "r") as f:
            lines = f.readlines()
        # Open output file for overwriting
        with open(output[0], "w") as f:
            # if there's only one line, keep it, otherwise need to trim the file
            if len(lines) > 1:
                # grab the refseq_assembly we want to download
                accession = ref_master.loc[ref_master['id']==wildcards.id]['refseq_accession'].values[0]
                for line in lines:
                    # if the accession is in the line, keep it
                    if accession in line.strip("\n"):
                        #print(accession)
                        f.write(line)
            else:
                for line in lines:
                    f.write(line)
'''
'''
# download the files from the list, renaming to id in the process
# unzip into .fasta
rule ref_ncbi_dl:
    input:
        REF_ASMBL+"{id}_esearch.txt"
    output:
        REF_ASMBL+"{id}/{id}.fna.gz"
    #params:
    #    intermediate1 = REF_ASMBL+"{id}/{id}.fna.gz",
    #    intermediate2 = REF_ASMBL+"{id}/{id}.fna"
    #group:
    #    'ncbi'
    #threads:
    #    6
    shell:
        "wget -O {output} -i {input}"
    #shell:
    #    """
    #    wget -O {params.intermediate1} -i {input}
    #    gunzip {params.intermediate1}
    #    mv {params.intermediate2} {output}
    #    """
'''

################################################################################
"""
Download all metagenome samples from Rahat's papers
"""

# Download
rule meta_ncbi_dl:
    output:
        o1 = "data/metagenomes/raw/{m_id}/{m_id}_pass_1.fastq.gz",
        o2 = "data/metagenomes/raw/{m_id}/{m_id}_pass_2.fastq.gz"
    shell:
        "fastq-dump \
        --outdir data/metagenomes_raw/{wildcards.m_id} \
        --gzip \
        --skip-technical \
        --readids \
        --read-filter pass \
        --dumpbase \
        --split-3 \
        --clip \
        {wildcards.m_id}"

'''
# Convert the fastq into fasta files
rule meta_fastq_to_fasta:
    input:
        i1 = "data/metagenomes_raw/{m_id}_pass_1.fastq",
        i2 = "data/metagenomes_raw/{m_id}_pass_2.fastq"
    output:
        o1 = "data/metagenomes_raw/{m_id}_1.fasta",
        o2 = "data/metagenomes_raw/{m_id}_2.fasta"
    params:
        awk_fastq_to_fasta = '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}'
    run:
        shell("seqret -sequence {input.i1} -outseq {output.o1}")
        shell("seqret -sequence {input.i2} -outseq {output.o2}")
'''


''' Downloading Biosample for Rylan

ids = ['SAMN12648957']

REF_ASMBL = "whateverpath/"

rule all:
     input:
        expand(REF_ASMBL+"{id}/{id}.fasta", id=ids)


# get input file list
rule ncbi_fna_list:
    output:
        temp(REF_ASMBL+"{id}_esearch.txt")
    shell:
        """
        esearch -db biosample -query {wildcards.id} \
          | elink -target assembly \
          | esummary \
          | grep "FtpPath_RefSeq" \
          | sed -r "s|.+>(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)<.+|\\1\\2/\\2_genomic.fna.gz|" \
          > {output}
        """

# download the files from the list, renaming to id in the process
# unzip into .fasta
rule ncbi_dl:
    input:
        REF_ASMBL+"{id}_esearch.txt"
    output:
        REF_ASMBL+"{id}/{id}.fasta"
    params:
        intermediate1 = REF_ASMBL+"{id}/{id}.fna.gz",
        intermediate2 = REF_ASMBL+"{id}/{id}.fna"
    shell:
        """
        wget -O {params.intermediate1} -i {input}
        gunzip {params.intermediate1}
        mv {params.intermediate2} {output}
        """
'''
