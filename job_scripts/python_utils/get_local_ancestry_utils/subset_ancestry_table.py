###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           subset_ancestry_table.py
###############################################################################


def subset_ancestry_table(ancestry_table, start_node, end_node):
    return ancestry_table[
        (ancestry_table["sample"] >= start_node)
        & (ancestry_table["sample"] <= end_node)
    ].copy()
