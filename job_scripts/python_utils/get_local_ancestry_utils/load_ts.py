###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           load_ts.py
###############################################################################


import tskit


def load_ts(ts_path):
    return tskit.load(ts_path)
