###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           load_and_subset_ts.py
###############################################################################


import numpy as np
import tskit


def load_and_subset_ts(ts_path, start_ind, end_ind):
    ts = tskit.load(ts_path)
    keep = []
    node_to_individual = np.array(ts.tables.nodes.individual, dtype=int)

    for ind_id in range(start_ind, end_ind + 1):
        keep.extend(ts.individual(ind_id).nodes)

    subset_ts, old_node_to_new_node = ts.simplify(
        keep, keep_unary=True, map_nodes=True
    )

    new_node_to_original_node = np.full(subset_ts.num_nodes, tskit.NULL, dtype=int)
    new_node_to_original_individual = np.full(
        subset_ts.num_nodes, tskit.NULL, dtype=int
    )

    # Invert the simplify node map so output rows can use original IDs.
    for old_node, new_node in enumerate(old_node_to_new_node):
        if new_node == tskit.NULL:
            continue
        new_node_to_original_node[new_node] = old_node
        new_node_to_original_individual[new_node] = node_to_individual[old_node]

    return subset_ts, new_node_to_original_node, new_node_to_original_individual
