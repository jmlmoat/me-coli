import pandas as pd
import os, sys
import seaborn as sns
import numpy as np
import ete3
from ete3 import TreeStyle, TextFace, faces, AttrFace, NodeStyle, ClusterTree
from PIL import Image, ImageDraw, ImageFont
from math import ceil
import argparse



def remove_path_from_name(tree):
    """
    removes filepath from the name of the genome
    """
    for node in tree.traverse("levelorder"):
        if node.name!='':
            node.name=node.name.split('/')[-1]
    return t


def make_legend(color_dict, col_name='nan', num_per_col=10, font_size=40, column_width=250):
    """
    Make a separate legend because ETE legends are the worst thing I've ever
    seen in my life oh my god
    """

    # Load a font for the legend
    cwd = os.getcwd()
    font_path = cwd+"/calibri.ttf"
    font = ImageFont.truetype(font_path,font_size)

    # If making a legend for serotype, need much larger col widths
    if 'sero' in col_name: column_width = 500

    # Number of unique entries in the dictionary
    num_unique = len(color_dict)

    # If all the entries will fit in one column, make the num_per_col equal
    # to the number of entries so we can strip tailing whitespace
    if num_unique <= num_per_col: num_per_col = num_unique

    # Based on number of entries per column, calculate the number of columns
    # we need for the legend
    num_cols = ceil(num_unique/num_per_col)
    # If the number of unique entries is less than our limit, use one column
    if num_cols == 0:
        width = column_width
    # If the number of unique entries over our limit, split into multiple cols
    else:
        width = column_width*num_cols

    # Whitespace padding around the border of the image
    horizontal_padding = font_size
    vertical_padding = font_size/4
    # Height of the image
    height = ((num_per_col+1)*font_size)+(vertical_padding*4)
    height = ceil(height)
    #   Each entry has the height of the font, so the total height is the
    #   number of entries multiplied by the font size, plus a top margin and
    #   bottom margin which are both equal to the font size

    # Create an image to write to, with the background colour as white
    img = Image.new('RGB', (width,height), color = (255,255,255))
    d = ImageDraw.Draw(img)

    # Distance from the left side for the current column; initially just the pad
    x_margin = horizontal_padding
    # Distance from the top of the img; initially 0 for the first entry in a col
    y_margin = vertical_padding
    # Current column number; initially on column 0
    curr_col = 0
    # For each entry in the dictionary, write the key to the image
    for index,key in enumerate(color_dict):
        # If we're on the last entry of the column
        if index % num_per_col == 0:
            # Reset the y index to write to the top of the next column
            y_margin = vertical_padding
            # Move horizonatally to the next column
            x_margin = horizontal_padding + column_width*curr_col
            # Increment the current column number
            curr_col = curr_col + 1
        # Move down by one unit to enter the next key
        y_margin = y_margin+font_size
        # Write the key to the image in the appropriate colour
        d.text((x_margin,y_margin), key, font=font, fill=color_dict[key])

    return img


def display_tree(tree,out_path,color_dict,num_per_col=10,ts_mode="c"):
    # Output directory

    if outgroup_name:
        outfile = "{}{}_outgroup_{}.png".format(outpath,plt,outgroup_name)
        outleg = "{}{}_outgroup_{}_legend.png".format(outpath,plt,outgroup_name)
    else:
        outfile = "{}{}.png".format(outpath,plt)
        outleg = "{}{}_legend.png".format(outpath,plt)

    #print(outdir,outtree,outleg)

    ts.mode = ts_mode
    if ts.mode == "c":
        ts.arc_start = -180
        ts.arc_span = 180

    #ts.legend.add_face(TextFace("STUFF", fsize=5000), column=1)

    #outleg = dfname+"_"+column+"_legend.png"
    # Save tree
    #t.show(tree_style=ts)
    t.render(outfile, dpi=int(600), w=int(5000), tree_style=ts)
    # Make legend separately because ETE sucks
    legend = make_legend(color_dict, col_name='clade', num_per_col=num_per_col)
    legend.save(outleg)
    #legend = make_legend(color_dict, col_name=column, num_per_col=num_per_legend_cols)
    #legend.save(outdir+outleg)
    return


def color_by_column(t, df, column, label_ref=False):
    """
    color tree by column in the metadata
    """

    # Change datatype to str for consistency
    df = df.astype({column: str})

    # Get the unique values in the column
    unique = list(df[column].unique())
    num_unique = len(unique)

    # Generate a colour palette with seaborn; one colour for each unique value
    pal = sns.color_palette("hls", num_unique)
    # Translate the colour palette into hex
    pal = pal.as_hex()
    # Make a dictionary mapping each unique name to a colour
    color_dict = dict(zip(unique, pal))
    # Make a dictionary with id as key, and the column info as the value
    id_dict = dict(zip(df['id'], df[column]))

    # Traverse through the tree
    for leaf in t.get_leaves():
        # Grab the id of the node
        id = leaf.name.split('/')[-1]
        # Grab the column info associated with the id; 'None' when id not present
        set = id_dict.get(id, None)
        # if it's a reference sample dont colour it (for now)
        # typo means k12 is labelled as k21 at the moment, oops
        ref_names = ['o2h6', 'o25h4_2', 'o25h4_1', 'o157h7_2', 'o157h7_1',
                     'k12', 'k21', 'atcc117755', 'Salmonella']
        if id in ref_names:
            print(id)
            if label_ref == True:
                ref_names = ['o2h6', 'o25h4_2', 'o25h4_1', 'o157h7_2', 'o157h7_1', 'k12', 'atcc117755']
                if leaf.name in ref_names:
                    name_face = TextFace(leaf.name, fgcolor='red', fsize=1000)
                    leaf.add_face(name_face, column=0, position='float')
        else:
            # Colour the background by value, as according to the color_dict
            ns = ete3.NodeStyle()
            ns['bgcolor'] = color_dict[set]
            leaf.set_style(ns)

    return t,color_dict


if __name__ == "__main__":
    """
    python tree_color.py -p collection
    python tree_color.py -p collection -g salmonella

    python tree_color.py -p collection -t binary_kmer_alignment_oneline.phy.treefile -m metadata_master_copy.tsv

    python tree_color.py -p collection -t binary_kmer_alignment_oneline.phy.treefile -m metadata_master_copy.tsv -g midpoint

    conda create -n figs-tree ete3 seaborn pandas pillow -c bioconda -c conda-forge -c anaconda

    python tree_color.py -p collection -t final_alignment.phy.treefile -m metadata_master_copy.tsv -g Salmonella -o test/

    python tree_color.py -p collection -t final_alignment.phy.treefile -m AHH.tsv -g Salmonella -o test/



    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--plot', help='plot type ie serotype', required = True)
    parser.add_argument('-t', '--treepath', help='path to tree', default="tree.parsimony.tre",required=False)
    parser.add_argument('-m', '--metadata', help='path to metadata master', default="metadata_master.tsv", required=False)
    parser.add_argument('-o', '--outpath', help='folder path for saving output', default="figures/", required=False)
    parser.add_argument('-g', '--outgroup', help='outgroup; default none', required=False)
    parser.add_argument('-c', '--numcol', help='number of entries per column in the legend', default=10, required=False)
    parser.add_argument('-l', '--labelref', help='label reference genomes, True or False', default=False, required=False)


    args = parser.parse_args()

    plt = args.plot
    tpath = args.treepath
    master_path = args.metadata
    outgroup_name = args.outgroup
    outpath = args.outpath
    num_per_col = int(args.numcol)
    label_ref = args.labelref

    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)


    # Load metadata
    master = pd.read_csv(master_path,sep='\t', low_memory=False)


    # define ts as the tree's style
    ts = ete3.TreeStyle()
    # load the tree
    t = ete3.Tree(tpath)
    # remove filepath from genome names
    t = remove_path_from_name(t)


    # set the outgroup
    if outgroup_name:
        #new_outgroup = t&outgroup_name
        #t.set_outgroup(new_outgroup)

        if outgroup_name == "midpoint":
            midpt = t.get_midpoint_outgroup()
            t.set_outgroup(midpt)
            print(midpt.name)
            new_outgroup = t&midpt.name
        else:
            # set the outgroup to the desired name
            new_outgroup = t&outgroup_name
            t.set_outgroup(new_outgroup)
            new_outgroup = t&outgroup_name


    t, color_dict = color_by_column(t,master,plt,label_ref)
    display_tree(t,outpath,color_dict,num_per_col)
    sys.exit()
