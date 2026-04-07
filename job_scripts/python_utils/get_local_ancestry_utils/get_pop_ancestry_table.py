###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_pop_ancestry_table.py
###############################################################################


import tspop


def get_pop_ancestry_table(ts, census_time):
    ancestry_table = tspop.get_pop_ancestry(
        ts, census_time=census_time
    ).ancestry_table.copy()
    return ancestry_table
