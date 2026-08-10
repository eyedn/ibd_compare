###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_ibd_segments.py
###############################################################################


import numpy as np
from .phase_of_node import phase_of_node


def write_ibd_segments(
    ts, tsv_path, chr_used, min_span_arg, max_time_arg,
    new_node_to_original_node, new_node_to_original_individual,
    population, population_start
):
    if max_time_arg > 0:
        ibd_segments = ts.ibd_segments(
            min_span=min_span_arg, max_time=max_time_arg, store_segments=True
        )
    else:
        ibd_segments = ts.ibd_segments(
            min_span=min_span_arg, store_segments=True
        )

    if len(new_node_to_original_node) != ts.num_nodes:
        raise ValueError("new_node_to_original_node length does not match ts")
    if len(new_node_to_original_individual) != ts.num_nodes:
        raise ValueError(
            "new_node_to_original_individual length does not match ts"
        )

    node_to_individual = np.array(ts.tables.nodes.individual)
    node_times = ts.tables.nodes.time
    samples = np.array(ts.samples(), dtype=int)
    ind_to_samples = {}
    for u in samples:
        ind = int(node_to_individual[u])
        ind_to_samples.setdefault(ind, []).append(int(u))
    for ind in ind_to_samples:
        ind_to_samples[ind].sort()

    with open(tsv_path, "w") as f:
        for (n1, n2), seg_list in ibd_segments.items():
            ind1 = int(new_node_to_original_individual[n1])
            ind2 = int(new_node_to_original_individual[n2])
            sample1 = ind1 - population_start + 1
            sample2 = ind2 - population_start + 1
            if sample1 < 1 or sample2 < 1:
                raise ValueError("tree individuals are outside population bounds")
            hap1 = phase_of_node(n1, node_to_individual, ind_to_samples)
            hap2 = phase_of_node(n2, node_to_individual, ind_to_samples)
            for seg in seg_list:
                left = int(seg.left)
                right = int(seg.right)
                node_mrca = int(new_node_to_original_node[seg.node])
                seg_len = right - left
                time_mrca = float(node_times[seg.node])
                f.write(
                    f"chr{chr_used}\t{left}\t{right}\t"
                    f"{population}_{sample1}\t{hap1}\t"
                    f"{population}_{sample2}\t{hap2}\t"
                    f"{node_mrca}\t{time_mrca}\t{seg_len}\n"
                )
