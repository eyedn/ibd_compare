###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           phase_of_node.py
###############################################################################


def phase_of_node(u: int, node_to_individual, ind_to_samples) -> int:
    ind = int(node_to_individual[u])
    nodes = ind_to_samples.get(ind, [])
    if len(nodes) >= 2:
        if u == nodes[0]:
            return 1
        if u == nodes[1]:
            return 2
    return -1
