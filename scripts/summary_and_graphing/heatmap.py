import ete3
from ete3 import TreeStyle, TextFace, faces, AttrFace, NodeStyle, Tree, ProfileFace
from ete3 import ClusterTree, faces
#from ete3.treeview.faces import add_face_to_node

import pandas as pd
import numpy as np

import os,sys





if __name__ == "__main__":


    freq_table_file = "feat_count_table.tsv"
    tree_file ="final_alignment.phy.treefile"


    with open(freq_table_file, 'r') as fh:
        matrix = fh.read()

    # make tree with the tree file and table
    t = ClusterTree(tree_file, text_array=matrix)

    circle_style = TreeStyle()
    circle_style.mode = "c"
    circle_style.arc_start = 0 #-180
    circle_style.arc_span = 359 #180
    circle_style.show_leaf_name = False #including the names spreads out the heatmap too much
    #circle_style.force_topology=True
    #circle_style.scale = 20 # does weird format
    #circle_style.optimal_scale_level = "full"
    #circle_style.complete_branch_lines_when_necessary = True

    t.show("heatmap",tree_style=circle_style)
    sys.exit()



    #----------------------------------------------------------------
    # ProfileFace test - not working

    array =  t.arraytable
    print(array)
    matrix_dist = [i for r in range(len(array.matrix))\
                   for i in array.matrix[r] if np.isfinite(i)]
    matrix_max = np.max(matrix_dist)
    matrix_min = np.min(matrix_dist)
    matrix_avg = matrix_min+((matrix_max-matrix_min)/2)

    print(matrix_max,matrix_min,matrix_avg)

    def mylayout(node):
        if node.is_leaf():
            profileFace  = ProfileFace(matrix_max, matrix_min, matrix_avg, width=200, height=14, style="heatmap")
            faces.add_face_to_node(profileFace, node, 0, aligned=True)

    ts = TreeStyle()
    ts.layout_fn = mylayout
    ts.mode = "c"
    ts.arc_start = 0 #-180
    ts.arc_span = 359 #180
    ts.show_leaf_name = False
    ts.force_topology=True
    #ts.draw_guiding_lines = True
    t.show(tree_style=ts)
