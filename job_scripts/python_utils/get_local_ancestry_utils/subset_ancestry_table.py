# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# subset_ancestry_table.py
# ______________________________________________________________________________



def subset_ancestry_table(ancestry_table, start_node, end_node):
    return ancestry_table[
        (ancestry_table["sample"] >= start_node)
        & (ancestry_table["sample"] <= end_node)
    ].copy()
