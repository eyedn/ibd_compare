# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# simulate_ooaa.py
# ______________________________________________________________________________


import msprime
import tskit


def simulate_ooaa(
    *,
    sample_sets,
    demography,
    contig,
    msprime_model,
    seed,
):
    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        recombination_rate=contig.recombination_map,
        sequence_length=contig.length,
        ploidy=2,
        random_seed=seed,
        model=msprime_model,
    )
    ts = msprime.sim_mutations(
        ts,
        rate=contig.mutation_rate,
        random_seed=seed + 1,
    )
    if not isinstance(ts, tskit.TreeSequence):
        raise TypeError("Expected a tskit.TreeSequence result")
    return ts
