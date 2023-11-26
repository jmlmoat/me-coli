"""
Clean & assemble raw WGS data.
"""

import pandas as pd

configfile: "config/config.yaml"

# Load the metadata master
master = pd.read_csv(config['metadata_master'],sep='\t',low_memory=False)
# filter our the ones that are already assembled
master = master.loc[ (master['download_key'] != 'asmbl_refseq') & (master['download_key'] != 'asmbl_genbank') ]
# get the id list for the ones that need to be assembled
ids = master['id'].tolist()


def get_file_id(id):
    cipars_strains[cipars_strains.id == 'id']['file_id'].values[0]
    return file_id



# Paths - output
CLEAN = "data/genomes/clean/" # cleaned/trimmed genomes
ASMBL = "data/genomes/asmbl/" # assembled genomes


#**************************************************
# ResFinder 4.0 paper's 390 genome val set
#predgroup = 'resf_val_data'
#resf_prj = "metadata/public/PRJNA6164_metadata.tsv"
#resf_df = pd.read_csv(resf_prj, sep='\t')
#ids = resf_df['id'].tolist()





rule all:
    input:
        # FastQC on raw reads
        #expand("data/genomes/raw_fastqc/{ds}/{id}/{id}_1_fastqc.zip", zip, ds=dss, id=ids),
        # Trim-galore
        #expand(CLEAN+"{id}/trim-galore/{id}_val_1.fq.gz", id=ids)
        # Flash
        #expand(CLEAN+"{id}/flash/{id}.extendedFrags.fastq.gz", id=ids)
        # Assembly
        expand("data/genomes/asmbl/{id}/scaffolds.fasta", id=ids)

subworkflow download:
    workdir: "../"
    snakefile: "download.smk"


################################################################################
# Helper Functions
################################################################################
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
            #return [download("data/genomes/raw/{id}/{id}_2.fastq.gz")]
            return ["data/genomes/raw/{id}/{id}_2.fastq.gz"]
    # if its a downloaded genome, call the subworkflow
    else:
        #return [download("data/genomes/raw/{id}/{id}_2.fastq.gz")]
        return ["data/genomes/raw/{id}/{id}_2.fastq.gz"]


'''
def zip_min_max(zip_path, fname):
    """
    Get min and max frag length from the fastqc output
    """
    from zipfile import ZipFile
    # Open the fastqc data file to find seqeunce length
    with ZipFile(zip_path, 'r') as zip:
        with zip.open(fname) as fp:
            # Line 8 has "Sequence Length #-#"
            # So grab this as the range and return it
            for i, line in enumerate(fp):
                if i == 8:
                    range = line.split()[2].decode("utf-8")
                    min = range.split('-')[0]
                    max = range.split('-')[1]
                    #print(min,max)
                elif i > 8:
                    break
    return(min,max)
'''

################################################################################
# Trim & QC
################################################################################

# fastqc on raw reads
'''
rule raw_fastqc:
    input:
        f1 = get_raw_input_1,
        f2 = get_raw_input_2
    output:
        "data/genomes/raw_fastqc/{id}/{id}_1_fastqc.zip",
        "data/genomes/raw_fastqc/{id}/{id}_2_fastqc.zip",
    params:
        outdir = "data/genomes/raw_fastqc/{id}"
    shell:
        "fastqc --outdir {params.outdir} {input.f1} {input.f2}"
'''

# trim-galore on raw reads
rule trim_galore:
    group: 'trim_galore'
    input:
        f1 = get_raw_input_1,
        f2 = get_raw_input_2
    output:
        CLEAN+"{id}/trim-galore/{id}_val_1.fq.gz",
        CLEAN+"{id}/trim-galore/{id}_val_2.fq.gz"
    params:
        outdir = CLEAN+"{id}/trim-galore/",
        cores = 4,
        quality = 10
    shell:
        "trim_galore \
        --quality {params.quality} \
        --cores {params.cores} \
        --paired \
        --gzip \
        --basename {wildcards.id} \
        -o {params.outdir} \
        {input.f1} \
        {input.f2}"
        #"--cores 4 could be a sweet spot, anything above has diminishing returns"
        #"--cores 4 is: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15"
# had to remove fastqc because of a stupid font error

# flash on trim-galore results
# using FLASH2 https://github.com/dstreett/FLASH2
rule flash:
    group: 'flash'
    input:
        f1 = CLEAN+"{id}/trim-galore/{id}_val_1.fq.gz",
        f2 = CLEAN+"{id}/trim-galore/{id}_val_2.fq.gz"
    output:
        CLEAN+"{id}/flash/{id}.extendedFrags.fastq.gz",
        CLEAN+"{id}/flash/{id}.notCombined_1.fastq.gz",
        CLEAN+"{id}/flash/{id}.notCombined_2.fastq.gz"
    params:
        outdir = CLEAN+"{id}/flash",
        zip_path = CLEAN+"{id}/trim-galore/{id}_val_1_fastqc.zip",
        fname = "{id}_val_1_fastqc/fastqc_data.txt",
        max_overlap = 400 # see below for explanation
    threads:
        4
    shell:
        "flash2 \
        --threads {threads} \
        --max-overlap={params.max_overlap} \
        --output-prefix={wildcards.id} \
        --output-directory={params.outdir} \
        --compress {input.f1} {input.f2}"

'''
Docs
    http://ccb.jhu.edu/software/FLASH/MANUAL
    https://github.com/dstreett/FLASH2

Using the FLASH2 defailt -M=65 gives this error
    WARNING: An unexpectedly high proportion of combined pairs (85.56%)
    overlapped by more than 65 bp, the --max-overlap (-M) parameter.  Consider
    increasing this parameter.  (As-is, FLASH is penalizing overlaps longer than
    65 bp when considering them for possible combining!)


-M or --max-overlap:
maxOverlap is the maximum overlap length expected in approximately 90% of read
pairs. It is by default set to 70bp, which works well for 100bp reads generated from
180bp library (normal distribution of fragment lengths is assumed). Overlaps longer
than maxOverlap are still considered as good overlaps, but the mismatch ratio
(explained below) is calculated over the maxOverlap rather than the true overlap
length. If you enter a value for maxOverlap, then the read length, fragmetn length
and standard deviaiton of fragment lengths that you enter will be ignored for
calculation of maxOverlap parameter. Default: 70bp.

-m or --min-overlap:
minOverlap is the minimum required overlap length between two reads to provide
a confident overlap. Default: 10bp.

-t or --threads:
The default number of combiner threads is the number of processors.
When multiple combiner threads are used, the order of the combined and
uncombined reads in the output files will be nondeterministic.  If you need to
enforce that the output reads appear in the same order as the input, you must
specify --threads=1.
'''



'''
run:
    min,max = zip_min_max(params.zip_path,params.fname)
    shell("flash2 \
    --max-overlap={max} \
    --output-prefix={wildcards.id} \
    --output-directory={params.outdir} \
    --compress {input.f1} {input.f2}")
'''

# fastqc on flash output - not used currently
rule flash_qc:
    group: 'flash_qc'
    input:
        CLEAN+"{id}/flash/{id}.extendedFrags.fastq.gz"
    output:
        CLEAN+"{id}/flash_qc/{id}.extendedFrags_fastqc.zip",
    params:
        outdir = CLEAN+"{id}/flash_qc/"
    shell:
        "fastqc --outdir {params.outdir} {input}"

################################################################################
# Assembly
################################################################################

# assemble with spades
rule spades:
    group: 'spades'
    input:
        f1 = CLEAN+"{id}/flash/{id}.notCombined_1.fastq.gz",
        f2 = CLEAN+"{id}/flash/{id}.notCombined_2.fastq.gz",
        f3 = CLEAN+"{id}/flash/{id}.extendedFrags.fastq.gz"
    params:
        outdir = ASMBL+"{id}/",
        max_mem = 48,
    output:
        ASMBL+"{id}/scaffolds.fasta"
    threads:
        16
    shell:
        "export OMP_NUM_THREADS={threads} && \
        spades.py \
        --careful \
        -m {params.max_mem} \
        --threads {threads} \
        -o {params.outdir} \
        --pe1-1 {input.f1} \
        --pe1-2 {input.f2} \
        --s2 {input.f3}"
        # -m <> is max memory in GB
        # -k
        # --tmp-dir /tmp
        # --careful I used this
        #    > as of dec30 2019, now there is --isolate option but dont want to
        #      rerun everything right now
        #    > --careful and --isolate are incomaptable
'''
def rename_header(fasta_in, fasta_out, prefix):
    from Bio import Seq, SeqIO
    fasta_in = str(fasta_in)
    fasta_out = str(fasta_out)
    contig_num = 1
    for record in SeqIO.parse(fasta_in, "fasta"):
        contig_seq = record.seq
        contig_seq = contig_seq._get_seq_str_and_check_alphabet(contig_seq)
        prefix = str(prefix)
        contig_header = ">lcl|{}|contig{}".format(prefix, contig_num)
        with open(fasta_out,'a') as fout:
            fout.write(contig_header)
            fout.write("\n")
            fout.write(contig_seq)
            fout.write("\n")
        contig_num += 1
    return 0

# add id to the head of every contig
rule add_header_id:
    input:
        ASMBL+DATASET+"/{id}/scaffolds.fasta"
    output:
        ASMBL+DATASET+"/{id}/{id}.fasta"
    params:
        rel_path_for_link = "../../ALL/",
        path = ASMBL + "ALL/"
    run:
        pre = wildcards.id
        rename_header(input,output,pre)
        shell("if [ ! -d {params.path} ]; then mkdir -p {params.path}; fi")
        #shell("cp {output} {params.outdir}")
        shell("ln -s ../../../{output} data/assemblies/ALL/")
        #shell("ln -s {output} {params.rel_path_for_link}")
'''
