# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# get_pop_ancestry_table.py
# ______________________________________________________________________________


import tspop


def get_pop_ancestry_table(ts, census_time):
    return tspop.get_pop_ancestry(ts, census_time=census_time).ancestry_table.copy()
