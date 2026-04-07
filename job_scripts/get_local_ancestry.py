###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_local_ancestry.py
###############################################################################


from sys import argv

from python_utils import log_msg
from python_utils.get_local_ancestry_utils import (
    build_sample_node_table,
    get_pop_ancestry_table,
    load_ts,
    subset_ancestry_table,
    write_local_ancestry,
)


def main():
    ts_path = argv[1]
    out_tsv = argv[2]
    chr_used = argv[3]
    pop = argv[4]
    start_ind = int(argv[5])
    end_ind = int(argv[6])
    admixture_time = float(argv[7])
    census_time_offset = float(argv[8])
    census_time = admixture_time + census_time_offset

    log_msg(f"loading ts for {pop}")
    ts = load_ts(ts_path)
    ancestry_table = get_pop_ancestry_table(ts, census_time)
    sample_node_rows = build_sample_node_table(ts)
    start_node = 2 * start_ind
    end_node = 2 * end_ind + 1
    ancestry_table = subset_ancestry_table(ancestry_table, start_node, end_node)

    log_msg(f"writing local ancestry to {out_tsv}")
    write_local_ancestry(ancestry_table, sample_node_rows, out_tsv, chr_used)
    log_msg("done")


if __name__ == "__main__":
    main()
