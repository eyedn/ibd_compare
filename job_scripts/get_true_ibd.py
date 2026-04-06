#!/usr/bin/env python3

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_true_ibd.py
###############################################################################

import tskit
import numpy as np
from datetime import datetime
from sys import argv

def log_msg(msg: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp} | {msg}", flush=True)

def phase_of_node(u: int, node_to_individual, ind_to_samples) -> int:
    ind = int(node_to_individual[u])
    nodes = ind_to_samples.get(ind, [])
    if len(nodes) >= 2:
        if u == nodes[0]:
            return 1
        if u == nodes[1]:
            return 2
    return -1 # unknown/unexpected


def main():
    ts_path = argv[1]
    tsv_path = argv[2]
    chr_used = argv[3]
    min_span_arg = float(argv[4])
    max_time_arg = float(argv[5])

    log_msg("loading ts")
    ts = tskit.load(ts_path)

    log_msg("extending haplotypes")
    ts = ts.extend_haplotypes()

    log_msg("finding ibd segements haplotypes")
    if max_time_arg > 0:
        ibd_segments = ts.ibd_segments(
            min_span=min_span_arg, max_time=max_time_arg, store_segments=True
            )
    else:
        ibd_segments = ts.ibd_segments(
            min_span=min_span_arg, store_segments=True
            )
    node_to_individual = np.array(ts.tables.nodes.individual)
    node_times = ts.tables.nodes.time

    # build robust phase mapping: within each individual, sort its sample nodes.
    samples = np.array(ts.samples(), dtype=int)
    ind_to_samples = {}
    for u in samples:
        ind = int(node_to_individual[u])
        ind_to_samples.setdefault(ind, []).append(int(u))
    for ind in ind_to_samples:
        ind_to_samples[ind].sort()

    # write ibd segments to tsv
    log_msg(f"writing IBD to {tsv_path}")
    
    with open(tsv_path, "w") as f:
        for (n1, n2), seg_list in ibd_segments.items():
            ind1 = int(node_to_individual[n1])
            ind2 = int(node_to_individual[n2])
            hap1 = phase_of_node(n1, node_to_individual, ind_to_samples)
            hap2 = phase_of_node(n2, node_to_individual, ind_to_samples)
            for seg in seg_list:
                left = int(seg.left)
                right = int(seg.right)
                node_mrca = int(seg.node)
                seg_len = right - left
                time_mrca = float(node_times[node_mrca])
                f.write(
                    f"chr{chr_used}\t{left}\t{right}\t"
                    f"tsk_{ind1}\t{hap1}\t"
                    f"tsk_{ind2}\t{hap2}\t"
                    f"{node_mrca}\t{time_mrca}\t{seg_len}\n"
                )
    log_msg("done")

if __name__ == "__main__":
    main()
