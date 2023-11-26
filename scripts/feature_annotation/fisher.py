
import pandas as pd
import os,sys
import numpy as np
import argparse
import scipy.stats as stats




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', help='the tsv from find_hits', required=True)
    parser.add_argument('-m', help='metadata master', required=True)
    parser.add_argument('-p', help='predictfor, eg SIRAMP', required=True)
    parser.add_argument('-o', help='output folder', required=True)
    parser.add_argument('-t', help='num top feats', required=True)
    args = parser.parse_args()

    infile = args.i
    metadata_file = args.m
    predfor = args.p
    outpath = args.o
    n_top_feats = args.t

    #sir_col = 'SIR_'+predfor.upper()
    sir_col = predfor.split('SIR')[-1]
    sir_col = 'SIR_'+sir_col

    # print whole dataframes
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    # Load the metadata master, only id and SIR column
    df = pd.read_csv(metadata_file,sep='\t',usecols=['id',sir_col])

    # Load the input file, only need two cols
    kmerdf = pd.read_csv(infile,sep='\t',usecols=['kmer','genome_id'])

    # Get the names of the top 10 kmers
    kmers = kmerdf.kmer.unique()

    # For each kmer, get the list of genome ids
    # and mark it as present or absent in the main df

    def pres_abs_kmer(d):
        if d['id'] in ids:
            return 1
        else:
            return 0

    # For each kmer
    for kmer in kmers:
        # Take the chunk of the df for this kmer
        tmp = kmerdf[kmerdf['kmer']==kmer]
        # Get the list of unique genome_ids which contain this kmer
        ids = tmp.genome_id.unique()
        # Mark the kmer as present or absent for each genome in df
        #df[kmer] = np.where(df['genome_id'] in ids)
        df[kmer] = df.apply(pres_abs_kmer,axis=1)

        # Save the PA table
        pa_file = "{}top_{}_feats_pres_abs.tsv".format(outpath,n_top_feats)
        df.to_csv(pa_file, sep='\t', index=False)

        # total num genomes with kmer present; not currently used
        # NOTE THAT THIS NUMBER INCLUDES GENOMES WITH NO SIR DATA
        #total = df[kmer].sum()
        #print (total)

    #print(df)

    # to make a dataframe after (faster than appending to a df)
    predfor_list = []
    kmer_list = []
    pvalue_list = []
    oddsratio_list = []
    con_table_list = []
    for kmer in kmers:
        predfor_list.append(predfor)
        kmer_list.append(kmer)

        # Groupby to get the counts for the contingency table
        gb = df.groupby([sir_col, kmer]).size().reset_index(name='count')
        #print(gb)

        # Extract the values for the contingency table
        # If there is no entry, then there are 0 entries for that combination

        # Number that have kmer & are S
        try:
            ks = gb.loc[ (gb[sir_col]=='S') & (gb[kmer]==1), 'count'].values[0]
        except:
            ks = 0
        # Number that have kmer & are R
        try:
            kr = gb.loc[ (gb[sir_col]=='R') & (gb[kmer]==1), 'count'].values[0]
        except:
            kr = 0
        # Number that dont have kmer & are S
        try:
            nks = gb.loc[ (gb[sir_col]=='S') & (gb[kmer]==0), 'count'].values[0]
        except:
            nks = 0
        # Number that dont have kmer & are R
        try:
            nkr = gb.loc[ (gb[sir_col]=='R') & (gb[kmer]==0), 'count'].values[0]
        except:
            nkr = 0

        '''
        num of genomes that....
                                   | are S | are R
        have kmer hit (blast)      | ks    | kr
        dont have kmer hit (blast) | nks   | nkr

        same result as transposed table

        contingency_table = np.array( [ [ks, kr], [nks, nkr] ] )
        '''

        # 2 tailed fisher's test
        oddsratio, pvalue = stats.fisher_exact([[ks, kr], [nks, nkr]], alternative='two-sided')

        pvalue_list.append(pvalue)
        oddsratio_list.append(oddsratio)
        con_table_list.append([[ks,kr],[nks,nkr]])

        #print("-------------------------------------")
        #print(kmer)
        #print(ks,kr,nks,nkr)
        #print(oddsratio)
        #print(pvalue)

        # Save the groupby table
        gb_file = "{}top_{}_feats_contingency.tsv".format(outpath,n_top_feats)
        gb.to_csv(gb_file, sep='\t', index=False)


    '''
    predfor | kmer | pvalue | oddsratio | contingency_table
    SIR_AMP | AAAA | 0.007  |  ...      | [ks,kr,nks,nkr]
    SIR_AMP | TTTT |
    '''


    # make a results df
    output_df_dict = {
           'predfor': predfor_list,
           'kmer': kmer_list,
           'pvalue': pvalue_list,
           'oddsratio': oddsratio_list,
           'contingency_table': con_table_list
           }
    output_df = pd.DataFrame.from_dict(output_df_dict)

    print(output_df)

    # Save the fisher table
    fisher_file = "{}top_{}_feats_fisher.tsv".format(outpath,n_top_feats)
    output_df.to_csv(fisher_file, sep='\t', index=False)
