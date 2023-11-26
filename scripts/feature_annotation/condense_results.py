import pandas as pd
import os,sys
import numpy as np
import argparse






def merge_df(df_path, drug):
    from statistics import stdev
    """
    Takes path to pandas df with the columns:
    [kmer, gene_up, dist_up, gene_down, dist_down, start, stop, genome_id, contig_name]
    and returns a df with the columns:
    [drug, kmer, gene_up, gene_up_count, avg_dist_up, gene_down, gene_down_count, avg_dist_down]
    """
    print(df_path)
    df = pd.read_pickle(df_path)
    hit_summary = []

    for kmer in set(df['kmer']):
        # filter for only a single kmer
        kmer_df = df[df['kmer']==kmer]

        for gene_direc, gene_dist in [['gene_up','dist_up'],['gene_down','dist_down']]:
            for gene in set(kmer_df[gene_direc]):
                if(len(gene)==0):
                    continue
                # filter for only a single gene
                gene_df = kmer_df[kmer_df[gene_direc]==gene]

                total_dist = 0

                for dist in gene_df[gene_dist]:
                        total_dist += abs(float(dist))

                count = gene_df.shape[0]
                average_dist = total_dist/count
                if(len(gene_df[gene_dist])<2):
                    std_dev = 0
                else:
                    std_dev  = stdev([abs(int(float(i))) for i in gene_df[gene_dist]])

                try:
                    gene, name = gene.split(':')
                except:
                    print("Gene: {}".format(gene))
                    print("{}".format(drug))
                    gene, carb, name = gene.split(':')

                hit_summary.append([drug, kmer, gene, name, count, average_dist, std_dev])

    return hit_summary



def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def add_value(df,col,value,id,cn,start,stop):
    df.loc[
        (df['genome_id']==id) &
        (df['contig_num']==cn) &
        (df['start']==start) &
        (df['stop']==stop),
        col
    ] = value
    return df



if __name__ == "__main__":
    '''
    python scripts/feature_annotation/condense_results.py -i feats/ori.pkl -o SIRAMC.tsv -p SIRAMC -m metadata/metadata_master.tsv

    sbatch -c 1 --mem 4G --wrap="python \
    scripts/feature_annotation/condense_results.py \
    -i results/wgs_standard/11mer/models/all/SIRAMC/nested_cv/1000feats/annotation/top_10_feats_hits.pkl \
    -o feats/SIRAMC_dense.tsv \
    -p SIRAMC \
    -m metadata/metadata_master.tsv"

    '''
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', help='input file', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-p', help='predict for', required=True)
    parser.add_argument('-m', help='path to metadata master', required=True)
    args = parser.parse_args()

    infile = args.i
    outfile = args.o
    pred = args.p
    master = args.m

    wgs_type = 'standard'

    # Load input file
    df = pd.read_pickle(infile)
    #df = pd.read_csv(infile,sep='\t')

    ############################################################################
    # rylan's stuff
    '''
    data = merge_df(infile,pred)
    data = np.asarray(data)
    what_amg = np.zeros((data.shape[0]), dtype = object)
    all_df = pd.DataFrame(
        data=np.concatenate(
                (data,np.asarray([what_amg]).T),axis=1),
                columns=[
                    'drug','kmer','gene','name',
                    'count','average_dist', 'std_dev' ,"AMG"
                ]
             )

    print(all_df)
    '''
    ############################################################################

    # Make a column containing just the contig number
    df['contig_num'] = df['contig_name'].str.split('_').str[-1]
    # Make sure the start and stop columns are in the right dtype
    df['start'] = df['start'].astype(float).astype(int)
    df['stop'] = df['stop'].astype(float).astype(int)
    df['dist_up'] = df['dist_up'].astype(float).astype(int)

    df = df.replace({'dist_down': {'-1': 0}})
    # replace the -1's with 0's
    df['dist_down'] = df['dist_down'].astype(float).astype(int)

    df[pred] = ''

    # Load the metadata master too look up SIR class
    # look up the SIR classification in master
    mdf = pd.read_csv(master,sep='\t')
    mdf = mdf.set_index('id')
    mdf.columns = mdf.columns.str.replace("_", "")
    # get list of all unique genome_ids
    unique_genome_ids = df.genome_id.unique()
    # for each genome in the list
    for id in unique_genome_ids:
        # add the SIR classification to the df
        df.loc[ (df['genome_id']==id), pred ] = mdf.at[id,pred]


    # Drop the rows with genomes that don't have a SIR designation
    # These weren't used in training the models, so just scrap them
    df.dropna(subset=[pred], inplace=True)
    # Reset index just so i can see the updated count
    df = df.reset_index(drop=True)

    #############################################################################
    # removing the _1 _2 etc to condense prokka more
    df['gene_up'] = df['gene_up'].str.replace('_1','')
    df['gene_up'] = df['gene_up'].str.replace('_2','')
    df['gene_up'] = df['gene_up'].str.replace('_3','')
    ############################################################################

    pred_df = df.copy()
    #print(df)


    df['card_gene'] = ''
    df['card_product'] = ''
    df['card_resistance'] = ''
    df['card_start'] = ''
    df['card_Stop'] = ''

    # get list of all unique genome_ids
    unique_genome_ids = df.genome_id.unique()
    # for each genome in the list
    for id in unique_genome_ids:
        # filter the indf for current genome
        id_df = df[df['genome_id']==id]

        # open the card.tsv and filter for current genome
        cardfile = "results/wgs_{}/abricate/{}/card.tsv".format(
            wgs_type, id )
        card = pd.read_csv(cardfile,sep='\t')
        card['contig_num'] = card['SEQUENCE'].str.split('_').str[1]
        # go through each row and update it with card data
        for i, row in id_df.iterrows():
            # get the start and stop points of the current kmer
            up = 0#row['dist_up']
            down = 0#row['dist_down']
            start = row['start']+up
            stop = row['stop']+down

            # get the number of the contig of this row
            cn = row['contig_num']
            # filter the card df for this number
            card_cn = card[card['contig_num']==cn]
            # keep track if we've already found an overlap
            # currently can only deal with 1 match
            found_overlap=False
            # Look at each entry in the filtered card df and check if
            # the start/stop overlap with the other df
            for j,card_row in card_cn.iterrows():
                # get the gene stop and start points
                card_start = card_row['START']
                card_stop = card_row['END']
                # Check if the ranges overlap with the other df
                overlap = get_overlap([start,stop],[card_start,card_stop])
                # If the regions overlap, then there's a match
                # and it needs to be added to the df
                if overlap != 0:
                    if found_overlap==True:
                        raise Exception("More than one hit found in CARD, need to handle this")
                    found_overlap = True
                    # update the df
                    df = add_value(
                        df,'card_gene',card_row['GENE'],
                        id,cn,start,stop)
                    df = add_value(
                        df,'card_product',card_row['PRODUCT'],
                        id,cn,start,stop)
                    df = add_value(
                        df,'card_resistance',card_row['RESISTANCE'],
                        id,cn,start,stop)
                    df = add_value(
                        df,'card_start',card_row['START'],
                        id,cn,start,stop)
                    df = add_value(
                        df,'card_stop',card_row['END'],
                        id,cn,start,stop)


    df = df.drop(columns=[pred])

    # need to insert temp values so that the nans arent disregarded
    # by the groupby & aggregate
    df['card_gene'] = df['card_gene'].fillna('tmp')
    df['card_product'] = df['card_product'].fillna('tmp')

    df = df.groupby(
            ['kmer', 'gene_up', 'card_gene', 'card_product']
         )['genome_id'].apply(list).reset_index(name='genome_id_list_rep')

    # remove duplicates
    #df['genome_id_list_rep'] = df['genome_id_list_rep'].fillna('tmp')
    df['genome_id_list'] = df['genome_id_list_rep'].apply(lambda x: list(set(x)))
    df = df.replace({'genome_id_list_rep':{'tmp': np.nan}})
    df['genome_count'] = df['genome_id_list'].str.len()
    #df['uniq_genome_count'] = df['genome_id_list'].str.nunique()

    # get rid of the row multiindexing
    #df = df.reset_index()
    # revert the tmp values
    df = df.replace({'card_gene':{'tmp': np.nan},'card_product':{'tmp': np.nan}})
    # nename columns appropriately
    df = df.rename(columns={'gene_up':'prokka_gene', 'genome_id':'genome_count'})
    # get rid of the column multiindexing
    #df.columns = df.columns.droplevel(1)
    # reorder the columns so count is second
    cols = df.columns.tolist()
    cols.insert(1, cols.pop(cols.index('genome_count')))
    df = df.reindex(columns= cols)


    # add count of # R, I, S

    #pred_df = pd.read_csv(infile,sep='\t')
    #pred_df = pd.read_pickle(infile)
    if 'sir' in pred.lower():
        pred_df = pred_df[['genome_id', pred]].drop_duplicates(keep='first').set_index('genome_id')
        for i, row in df.iterrows():
            id_list = row['genome_id_list']
            S = 0
            I = 0
            R = 0
            U = 0
            for id in id_list:
                sir = pred_df.at[id,pred]
                if sir=='S': S+=1
                elif sir=='I': I+=1
                elif sir=='R': R+=1
                else: U+=1
            df.at[i,'S'] = S
            df.at[i,'I'] = I
            df.at[i,'R'] = R
            df.at[i,'U'] = U

        cols = df.columns.tolist()
        cols.insert(2, cols.pop(cols.index('U')))
        cols.insert(2, cols.pop(cols.index('I')))
        cols.insert(2, cols.pop(cols.index('S')))
        cols.insert(2, cols.pop(cols.index('R')))
        df = df.reindex(columns= cols)

    elif 'clade' in pred.lower():
        pred_df = pred_df[['genome_id', pred]].drop_duplicates(keep='first').set_index('genome_id')

        for i, row in df.iterrows():
            id_list = row['genome_id_list']

            #print(df)

            #remove duplicates
            #id_list = list(dict.fromkeys(id_list))

            #clade_list = []
            #for id in id_list:
            #    clade_list = clade_list.append(pred_df.at[id,pred])

            clade_list = [ pred_df.at[id,pred] for id in id_list ]

            #from collections import Counter
            #clade_counts = Counter(clade_list)
            counts = dict()
            for j in clade_list:
              counts[j] = counts.get(j, 0) + 1

            print(i)
            print(id_list)
            print(clade_list)
            print(counts)
            print()


            for clade,clade_count in counts.items():
                df.at[i, str(clade)] = clade_count


        movetoback = [ 'genome_id_list', 'genome_id_list_rep' ]
        df = df[ [ col for col in df.columns if col not in movetoback ] + movetoback ]

        #cols = df.columns.tolist()
        #cols.insert(-1, cols.pop(cols.index('genome_id_list')))
        #cols.insert(-1, cols.pop(cols.index('genome_id_list_rep')))
        #df = df.reindex(columns= cols)

    df.to_csv(outfile,sep='\t', index=False)
