# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# sim_ooaa.R
# ______________________________________________________________________________


import stdpopsim
from python_utils.sim_ooaa_utils import (
    build_ooaa_demography,
    build_sample_sets,
    get_contig,
    parse_args,
    simulate_ooaa,
)


def main():
    args = parse_args()

    species = stdpopsim.get_species("HomSap")
    contig = get_contig(species, args.chromosome, args.genetic_map)
    demography, _ = build_ooaa_demography(
        admixture_time=args.admixture_time,
        afr_proportion=args.afr_proportion,
        eur_proportion=args.eur_proportion,
        admix_initial_size=args.admix_initial_size,
        admix_growth_rate=args.admix_growth_rate,
        census_time_offset=args.census_time_offset,
    )
    sample_sets = build_sample_sets(args.sample_size)
    ts = simulate_ooaa(
        sample_sets=sample_sets,
        demography=demography,
        contig=contig,
        msprime_model=args.msprime_model,
        seed=args.seed,
    )
    ts.dump(args.output)


if __name__ == "__main__":
    main()
