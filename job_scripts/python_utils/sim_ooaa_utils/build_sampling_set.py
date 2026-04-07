###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_sampling_set.py
###############################################################################


import msprime


def build_sample_sets(sample_size):
    return [
        msprime.SampleSet(sample_size, population="AFR", time=0, ploidy=2),
        msprime.SampleSet(sample_size, population="EUR", time=0, ploidy=2),
        msprime.SampleSet(sample_size, population="ADMIX", time=0, ploidy=2),
    ]
