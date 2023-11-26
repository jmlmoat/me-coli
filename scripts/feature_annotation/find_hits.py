import pandas as pd
import gffpandas.gffpandas as gffpd
import os,sys
import skbio.io
import numpy as np
import argparse
#import pickle5 as pickle

def find_locs(kmer, blast_df):
    # where to store all the location hits
    locs = []

    #print(blast_df)
    #sys.exit()

    # filter blast results to just our kmer of interest
    kmer_df = blast_df[blast_df['qseqid'] == kmer]
    # make sure it's not an empty df
    assert(kmer_df.shape[0]>0)
    # for every kmer hit in the genome, append the location
    for i in range(kmer_df.shape[0]):
        send = kmer_df['send'].values[i]
        sstart = kmer_df['sstart'].values[i]
        # note that for my data, the genome_id and contig_name are the same,
        # and are just the sseqid
        genome_id = kmer_df['sseqid'].values[i].rsplit('_',1)[0].split('_NZ')[0]
        contig_name = kmer_df['sseqid'].values[i]
        locs.append([send,sstart,genome_id,contig_name])

    # locs is 2d list, each row is (start, stop, genome_id, contig_name)
    return locs

def find_gene(start, stop, genome_id, contig_name, prokka_loc):
    gene_up = ''
    dist_up = -1
    gene_down = ''
    dist_down = -1

    # Prokka renames everything with {locustag}_{contig_number}.
    # The contig numbers match the original values.
    # Extract the contig number from the contig name which is in the format
    #  {id}_{contig_num} eg ERR121_1
    contig_num = contig_name.split('_')[-1]

    # gff location is what was passed into the function
    gff_loc = prokka_loc

    # Scan through contigs until the correct number is found
    with open("{0}{1}/{1}.gff".format(gff_loc,genome_id)) as file:
        for line in file:
            # if the matching contig number is found
            if("_{} ".format(contig_num) in line):
                # store the contig name
                prokka_contig = line.split(" ")[1]
                #print(prokka_contig)
                break
    # If the contig number wasnt found
    if('prokka_contig' not in locals()):
        print("Contig number {2} and contig name {3} not located in {0}{1}/{1}.gff".format(gff_loc,genome_id, contig_num, contig_name))
        return [gene_up, dist_up, gene_down, dist_down]


    # Load the prokka annotated genome dataframe
    # columns are:
    #  ['seq_id','source','type','start','end',
    #   'score','strand','phase','attributes']
    gff_df = pd.read_pickle(prokka_loc+genome_id+'_df.pkl')

    # keep only the contig the kmer was found in and only show coding sequences (Prodigal)
    contig_df = gff_df[gff_df['seq_id']==prokka_contig]
    contig_df = contig_df[contig_df['type']=='CDS']


    start = int(start)
    stop = int(stop)

    df_length = contig_df.values.shape[0]

    # find the nearest gene/genes
    # for every gene found by prokka, does it contain the kmer or is it near?
    for gene_num, gene_anno in enumerate(contig_df.values):
        gene_start = int(gene_anno[3])
        gene_stop = int(gene_anno[4])

        try:

            if(start > gene_stop):
                if(gene_num==df_length-1):
                    # we are after the last gene
                    gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                    dist_up = start - gene_stop
                    if(gene_dict['product']=='hypothetical protein'):
                        gene_up = "HypoProt:hypothetical protein"
                    else:
                        gene_up = gene_dict['gene']+':'+gene_dict['product']
                    break

                # we are not at the correct gene yet
                continue
            elif(stop < gene_start):
                # the kmer is before the current gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_down = gene_start-stop
                if(gene_dict['product']=='hypothetical protein'):
                    gene_down = "HypoProt:hypothetical protein"
                else:
                    try:
                        gene_down = gene_dict['gene']+':'+gene_dict['product']
                    except KeyError:
                        gene_down = 'none:'+ dict(i.split('=') for i in gene_anno[8].split(';'))['product']

                prev_gene_anno = contig_df.values[gene_num-1]

                gene_dict = dict(i.split('=') for i in prev_gene_anno[8].split(';'))
                dist_up = start - prev_gene_anno[4]
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start >= gene_start and stop <= gene_stop):
                # the kmer is inside of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = 0
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start <= gene_stop <= stop):
                # kmer hanging over right end of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = stop-gene_stop
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            elif(start <= gene_start <= stop):
                # kmer hanging over the left end of a gene
                gene_dict = dict(i.split('=') for i in gene_anno[8].split(';'))
                dist_up = gene_start-start
                if(gene_dict['product']=='hypothetical protein'):
                    gene_up = "HypoProt:hypothetical protein"
                else:
                    gene_up = gene_dict['gene']+':'+gene_dict['product']
                break

            else:
                raise Exception("Unexpected kmer start: {} stop: {} in relation to gene start: {} stop: {}".format(start, stop, gene_start, gene_stop))
        except KeyError:
            gene_up = 'none:'+ dict(i.split('=') for i in gene_anno[8].split(';'))['product']
            break

    return [gene_up, dist_up, gene_down, dist_down]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-blast_path', help='', required=True)
    parser.add_argument('-top_f', help='', required=True)
    parser.add_argument('-prokka_loc', help='', required=True)
    parser.add_argument('-kmer_len', help='', required=True)
    parser.add_argument('-out', help='', required=True)
    args = parser.parse_args()

    blast_path = args.blast_path
    top_f = args.top_f
    prokka_loc = args.prokka_loc
    kmer_length = args.kmer_len
    outfile = args.out

    top_feats = np.load(top_f, allow_pickle = True)

    with open(blast_path) as fh:
        blast_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True)

    ############################################################################
    #print(blast_df)
    #sys.exit()
    ############################################################################

    # each row in gene hits will be [kmer, gene_up, dist_up, gene_down, dist_down, start, stop, genome_id, contig_name]
    gene_hits = []
    for kmer in top_feats:
        #if kmer_length != '11':
        #    kmer = kmer[0]
        # locs is 2d list, each row is (start, stop, genome_id, contig_name)
        locs = find_locs(kmer, blast_df)

        for loc in locs:
            # gene_info is 1D list: gene_up, dist_up, gene_down, dist_down
            gene_info = find_gene(*loc, prokka_loc)
            gene_hits.append([kmer]+gene_info+loc)



    # save gene hits as a dataframe
    hits_df = pd.DataFrame(data = np.asarray(gene_hits),columns = ['kmer', 'gene_up', 'dist_up', 'gene_down', 'dist_down', 'start', 'stop', 'genome_id', 'contig_name'])
    hits_df.to_pickle(outfile)
    # also save as a .tsv
    alt_outfile = outfile.split('.')[0]+'.tsv'
    hits_df.to_csv(alt_outfile, sep='\t', index=False)
