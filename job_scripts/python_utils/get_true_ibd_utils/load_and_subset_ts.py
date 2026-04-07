###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           load_and_subset_ts.py
###############################################################################


import tskit


def load_and_subset_ts(ts_path, start_ind, end_ind):
    ts = tskit.load(ts_path)
    keep = []
    for ind_id in range(start_ind, end_ind + 1):
        keep.extend(ts.individual(ind_id).nodes)
    return ts.simplify(keep, keep_unary=True)
