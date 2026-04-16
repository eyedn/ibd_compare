###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_true_ibd.py
###############################################################################


from sys import argv

from python_utils import log_msg
from python_utils.get_true_ibd_utils import (
    extend_haplotypes,
    load_and_subset_ts,
    write_ibd_segments
)


def main():
    ts_path = argv[1]
    tsv_path = argv[2]
    chr_used = argv[3]
    min_span_arg = float(argv[4])
    max_time_arg = float(argv[5])
    start_ind = int(argv[6])
    end_ind = int(argv[7])

    log_msg("loading ts")
    ts, new_node_to_original_node, new_node_to_original_individual = load_and_subset_ts(ts_path, start_ind, end_ind)

    log_msg("extending haplotypes")
    ts = extend_haplotypes(ts)

    log_msg(f"writing IBD to {tsv_path}")
    write_ibd_segments(
        ts, tsv_path, chr_used, min_span_arg, max_time_arg, 
        new_node_to_original_node, new_node_to_original_individual
        )
    log_msg("done")

if __name__ == "__main__":
    main()
