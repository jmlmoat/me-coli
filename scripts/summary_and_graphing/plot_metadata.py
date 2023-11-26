import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ruamel.yaml
import os, sys
import argparse
import warnings


def plot_sir_class_freq(out_prefix,sheet,outdir="figures/metadata_summary/sir_freq/"):
    """
    Plots the SIR class frequencies for the specified sheet and drug(s)
    ____________________________________________________________________
    ____________________________________________________________________
    ____________________________________________________________________
    :dataset: name of the dataset (str)
    :sheet: pandas dataframe of the metadata sheet; should match the dataset
    :outdir: directory to save the plot in
    :return:
    """

    # Use default seaborn formatting
    sns.reset_orig()
    sns.set()

    # Order we want the graph bars to be in
    bar_order = ['S','I','R']

    '''
    # Make one graph for each drug
    for drug in drugs:
        sns.catplot(x=drug, kind='count', data=sheet, order=bar_order)
        #plt.show()
    '''

    ## Make a single graph with all the drugs, with catplot
    # Slice off just the SIR columns from the dataframe
    sir_df = sheet.filter(like='sir_')
    # Make uppercase again so the labels look nice
    sir_df.columns = map(str.upper, sir_df.columns)
    # Sum the number of S, I, and R occurences for each drug
    sir_df = sir_df.apply(pd.Series.value_counts)
    # Rename columns from SIR_DRUG to DRUG
    sir_df.rename(columns=lambda x: x.split('_')[1], inplace= True)
    # Manipulate data for graphing with catplot:
    #   transpose the df
    sir_df = sir_df.transpose()
    #   give the column of drugs the title 'drug'
    sir_df = sir_df.rename_axis('Drug')
    #   move the index (drug) column to the first column
    sir_df = sir_df.reset_index()

    #print(sir_df)

    # Reshape/melt sir_df for catplot
    sir_df = pd.melt(sir_df, id_vars='Drug', value_vars=bar_order, \
                     var_name='SIR', value_name='Count')

    # Plot & format
    g = sns.catplot(x='Drug', y='Count', hue='SIR', hue_order=bar_order, \
                    data=sir_df, kind='bar', height=5, aspect=2)

    g.set_xticklabels(rotation=90)
    #g.set(ylim=(0, 2000))

    #g.fig.suptitle(dataset+' '+' SIR Frequencies')

    ############################################################################
    ## Testing options to make the graph prettier
    #g.map(plt.subplots_adjust(hspace=20, wspace=20))
    #g.set_axis_labels('Drug','Count')
    #g.set_xticklabels(['S','I','R'])
    #g.despine(left=True)
    '''
    test = sheet.filter(like='SIR_')
    test.rename(columns=lambda x: x.split('_')[1], inplace= True)
    test = test.transpose()
    test = test.rename_axis('Drug')
    test = test.reset_index()
    test = pd.melt(test, id_vars='Drug', value_vars=bar_order, \
                     var_name='SIR', value_name='Count')
    print(test)
    t = sns.countplot(x='Drug', data=test)
    #t.set_xticklabels(bar_order)
    '''
    ############################################################################

    outfile = outdir+out_prefix+"_sir_freq.png"
    outdir = outdir.lower()
    outfile = outfile.lower()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outfile, dpi=300)

    # clear plots
    plt.clf()

    return

def plot_mic_class_freq(dataset,sheet,outdir="figures/metadata_summary/mic_freq/"):
    """
    ____________________________________________________________________
    ____________________________________________________________________
    :dataset: name of the dataset (str)
    :sheet: pandas dataframe of the metadata sheet; should match the dataset
    :outdir: directory to save the plot in
    :return:
    """

    # Use default seaborn formatting
    sns.set()
    sns.reset_orig()

    # Slice off just the MIC columns from the dataframe
    mic_df = sheet.filter(like='mic_')
    # Make uppercase again so the labels look nice
    mic_df.columns = map(str.upper, mic_df.columns)

    # load the config file to get the class list for the drugs
    configfile = "config/config.yaml"
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(configfile))
    drug_classes = config['drug_class_ranges']

    # For each drug in the dataframe
    mic_drugs = list(mic_df.columns.values)
    for mic_drug in mic_drugs:
        drug = mic_drug.split('_')[-1]

        if not mic_df[mic_drug].isnull().all():
            # If the drug class range is not None
            if drug_classes[drug]:
                # Get the class range from the config file to specify the order of
                # the bars on the graph. Need to convert to strings.
                class_order = drug_classes[drug]
                class_order = [str(x) for x in class_order]

                # Plot the mic frequencies
                ax = sns.countplot(x=mic_drug, data=mic_df, order=class_order)
                # Rename the x and y axis
                ax.set(xlabel='MIC class', ylabel='Count')
                # Add a figure title
                ax.set_title(dataset+' '+drug+' MIC Frequencies')
                # Set y axis limit to 500
                #ax.set(ylim=(0, 500))

                ## Add counts to the top of each bar
                # For each bar in the graph
                for i,p in enumerate(ax.patches):
                    # get the height of tha bar
                    height = p.get_height()
                    # if height is nan set it to 0
                    if np.isnan(height): height = 0
                    # how high above the bar to put the count text
                    _, ymax = plt.ylim()
                    bump = ymax * 0.01
                    # put the height above the bar
                    ax.text(p.get_x() + p.get_width()/2., p.get_height() + bump, str(int(height)), ha="center")

            # Save figure
            outfile = outdir+dataset+"_"+drug+"_mic_freq.png"
            outdir = outdir.lower()
            outfile = outfile.lower()
            os.makedirs(outdir, exist_ok=True)
            plt.savefig(outfile, dpi=300)

            # Clear plot
            plt.clf()
    return

def plot_sero_freq(dataset,sheet,column,maxbars=50,outdir="figures/metadata_summary/serotype_freq/"):
    """
    ____________________________________________________________________
    ____________________________________________________________________
    :dataset: name of the dataset (str)
    :sheet: pandas dataframe of the metadata sheet; should match the dataset
    :outdir: directory to save the plot in
    :return:
    """

    # Use default seaborn formatting
    sns.set()
    sns.reset_orig()

    # Slice off just the serotype column from the dataframe
    sero_df = sheet[column]

    '''
    ax = sns.countplot(x=sero_df, data=sero_df) # plot the frequencies
    ax.set(xlabel='Serotype', ylabel='Count') # rename x and y axes
    ax.set_title(dataset+' Serotype Frequencies') # add figure title
    plt.xticks(rotation=90, fontsize=10) # rotate x axis labels
    '''

    serovar_counts = sheet[column].value_counts()
    serovar_counts = serovar_counts[:maxbars,]
    i1_indx = 100
    for i, val in enumerate(serovar_counts.index):
        if val == 'I 1':
            i1_indx = i


    plt.figure(figsize=(7,4))

    ax = sns.barplot([x for i,x in enumerate(serovar_counts.index) if i!=i1_indx], [x for i,x in enumerate(serovar_counts.values) if i!=i1_indx])
    plt.xticks(rotation=90, fontsize=10)
    plt.tight_layout()

    ''' skipping numbering the top of columns for now; looks cluttered
    ## Add counts to the top of each bar
    # For each bar in the graph
    for i,p in enumerate(ax.patches):
        # get the height of tha bar
        height = p.get_height()
        # if height is nan set it to 0
        if np.isnan(height): height = 0
        # how high above the bar to put the count text
        _, ymax = plt.ylim()
        bump = ymax * 0.01
        # put the height above the bar
        ax.text(p.get_x() + p.get_width()/2., p.get_height() + bump, str(int(height)), ha="center")
    '''

    # Save figure
    outfile = outdir+dataset+"_"+column+"_freq.png"
    outdir = outdir.lower()
    outfile = outfile.lower()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outfile, dpi=300)

    # Clear plot
    plt.clf()

    return

def plot_freq_by_col(dataset,sheet,column,outdir="figures/metadata_summary/"):
    """
    ____________________________________________________________________
    ____________________________________________________________________
    :dataset: name of the dataset (str)
    :sheet: pandas dataframe of the metadata sheet; should match the dataset
    :outdir: directory to save the plot in
    :return:
    """

    # Use default seaborn formatting
    sns.set()
    sns.reset_orig()

    # Slice off just the one column from the dataframe
    df = sheet[column].astype(str)

    # remove nan if it's not wanted
    #df = pd.Series(np.asarray(df)[[i !='nan' for i in df]])

    # this combination needs more space
    if column == 'source' and dataset.upper() == 'ALL':
        plt.figure(figsize=(9,4))
    if column == 'country' and dataset.upper() == 'ALL':
        plt.figure(figsize=(9,9))
        plt.xticks(rotation=90, fontsize=8)

    if column == 'origin':
        plt.figure(figsize=(9,9))
        plt.xticks(rotation=90, fontsize=8)

    if column == 'year':
        #plt.figure(figsize=(9,9))
        plt.xticks(rotation=90, fontsize=8)

    # Plot the frequencies
    #ax = sns.countplot(x=df)
    #in highest to lowest count:
    ax = sns.countplot(x=df, order=df.value_counts().index)

    # Rename the x and y axis
    ax.set(xlabel=column.capitalize(), ylabel='Count')

    # need to change this later
    if column=='esbl': ax.set(xlabel=column.upper(), ylabel='Count')

    # Add a figure title
    ax.set_title(dataset+' '+column.capitalize()+' Frequencies')

    # need to change this later
    if column=='esbl': ax.set_title(dataset+' '+column.upper()+' Frequencies')

    ## Add counts to the top of each bar
    # For each bar in the graph
    for i,p in enumerate(ax.patches):
        # get the height of tha bar
        height = p.get_height()
        # if height is nan set it to 0
        if np.isnan(height): height = 0
        # how high above the bar to put the count text
        _, ymax = plt.ylim()
        bump = ymax * 0.01
        # put the height above the bar
        ax.text(p.get_x() + p.get_width()/2., p.get_height() + bump, str(int(height)), ha="center")



    # Save figure
    out = outdir+column+"_freq/"
    outfile = out+dataset+"_"+column+"_freq.png"
    outdir = outdir.lower()
    outfile = outfile.lower()
    os.makedirs(out, exist_ok=True)
    plt.savefig(outfile, dpi=300)

    # Clear plot
    plt.clf()
    return

if __name__ == "__main__":
    """
    $ python src/plot_metadata -ds DATASET -plt PLOT_TYPE
    $ python scripts/summary_and_graphing/plot_metadata.py -c country -p sir_amc -i canada thailand


     python scripts/summary_and_graphing/plot_metadata.py -m metadata/metadata_master_standard.tsv -c collection -i CIPARS BCRC ENA -p serotype_ec

     python scripts/summary_and_graphing/plot_metadata.py -m metadata/metadata_master_standard.tsv -c collection -i all -p serotype_ec -b 50

     python scripts/summary_and_graphing/plot_metadata.py -m metadata/metadata_master.tsv -c collection -i all -p country

    python scripts/summary_and_graphing/plot_metadata.py -m metadata/metadata_master.tsv -c collection -i all -p sir



    """
    # Get dataset and plot from the cmd line
    parser = argparse.ArgumentParser(description="Plots metadata information")
    parser.add_argument(
        '-c',
        help="column title in the tsv for the data filtration eg country, year, etc",
        required=True)
    parser.add_argument(
        '-p',
        help='plot column name exactly OR mic OR sir',
        required=True)
    parser.add_argument(
        '-i',
        help="the items to pull out eg canada etc; all for use all items", nargs='+',
        required=True)
    parser.add_argument(
        '-b',
        help="number of bars to display; only use for serotype etc", required=False)
    parser.add_argument(
        '-m',
        help="metadata tsv",
        required=False,
        default="metadata/metadata_master.tsv")
    args = parser.parse_args()

    data_col_name = args.c.lower()
    metadata_master_file = args.m
    data_col_items = args.i
    data_col_items = [i.lower() for i in data_col_items]
    plt_type = args.p.lower()
    if args.b:
        max_bars = int(args.b)
    else:
        max_bars = 50

    # file name prefix
    out_prefix = "_".join(data_col_items)


    # Load configfile for the list of valid datasets
    configfile = "config/config.yaml"
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(configfile))


    # Load the metadata master file
    #df = pd.read_csv(config['metadata_master'],sep='\t')
    df = pd.read_csv(metadata_master_file,sep='\t')

    # lowercase col names
    df.columns = map(str.lower, df.columns)

    # Filter the file for the desired column and items if not plotting all
    if data_col_items[0] != 'all':
        df = df[df[data_col_name].str.lower().isin(data_col_items)]


    # valdi plot options
    #valid_plt = ['sir','mic','serotype', 'country', 'year', 'source', 'esbl']
    valid_plts = ['sir','mic','serotype', 'collection', 'year', 'source_bin',
                  'country','region','phylogroup','esbl',
                  'serotype_ec','otype_ec','htype_ec',
                  'phylogroup_ct','esbl_abr','origin']
    # Check if the plot option is valid
    #plt_type = args.plt.lower()
    if plt_type not in valid_plts:
        raise Exception('{} not valid plot type'.format(plt_type))

    # Load the master metadata Excel sheet and all of the sheets
    #metadata_master = pd.ExcelFile('metadata/metadata_master.xlsx')
    #metadata_master = pd.read_csv(config['metadata_master'],sep='\t')


    # plot the metadata depending on the specified type
    if plt_type == 'sir':
        #plot_sir_class_freq(dataset,pd.read_excel(metadata_master, sheet_name=dataset))
        plot_sir_class_freq(out_prefix,df)

    elif plt_type == 'mic':
        # currently there's no MIC data for datasets other than CIPARS
        # so if using a diff dataset raise an error
        #if len(data_col_items)>1 or data_col_items[0] != 'cipars':
        #    raise Exception('{} has no MIC data. Currently only CIPARS has MIC data.'.format(data_col_items[0]))
        #plot_mic_class_freq(dataset,pd.read_excel(metadata_master, sheet_name=dataset))
        plot_mic_class_freq(out_prefix,df)
    elif plt_type in ['serotype','serotype_ec','otype_ec','htype_ec']:
        #plot_sero_freq(dataset,pd.read_excel(metadata_master, sheet_name=dataset))
        plot_sero_freq(out_prefix,df,plt_type,max_bars,outdir="figures/metadata_summary/"+plt_type+"_freq/")
    else:
        #plot_freq_by_col(dataset,pd.read_excel(metadata_master, sheet_name=dataset),plt_type)
        #plot_freq_by_col(dataset,metadata_master,plt_type)
        plot_freq_by_col(out_prefix,df,plt_type)
