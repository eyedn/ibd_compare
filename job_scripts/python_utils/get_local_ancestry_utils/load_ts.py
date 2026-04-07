# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# load_ts.py
# ______________________________________________________________________________


import tskit


def load_ts(ts_path):
    return tskit.load(ts_path)
