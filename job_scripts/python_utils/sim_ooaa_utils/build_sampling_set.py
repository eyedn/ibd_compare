# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# build_sampling_set.py
# ______________________________________________________________________________


import msprime


def build_sample_sets(sample_size):
    return [
        msprime.SampleSet(sample_size, population="AFR", time=0, ploidy=2),
        msprime.SampleSet(sample_size, population="EUR", time=0, ploidy=2),
        msprime.SampleSet(sample_size, population="ADMIX", time=0, ploidy=2),
    ]
