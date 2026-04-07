###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_sample_node_table.py
###############################################################################


def build_sample_node_table(ts):
    rows = []
    for ind_id in range(ts.num_individuals):
        nodes = sorted(int(u) for u in ts.individual(ind_id).nodes)
        for hap, node in enumerate(nodes, start=1):
            rows.append((node, ind_id, hap))
    return rows
