###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_3T.py
###############################################################################


import stdpopsim
from python_utils.sim_3T_utils import (
    build_3T_demography,
    build_sample_sets,
    get_contig,
    parse_args,
    simulate_3T
)


def main():
    args = parse_args()

    species = stdpopsim.get_species("HomSap")
    contig = get_contig(species, args.chromosome, args.genetic_map)
    demography, _ = build_3T_demography(
        generation_time=args.generation_time,
        mutation_rate=args.mutation_rate,
        t_af_years=args.t_af_years,
        t_ooa_years=args.t_ooa_years,
        t_eu0_years=args.t_eu0_years,
        t_eg_years=args.t_eg_years,
        r_eu0=args.r_eu0,
        r_eu=args.r_eu,
        r_af=args.r_af,
        n_a=args.n_a,
        n_af1=args.n_af1,
        n_b=args.n_b,
        n_eu0=args.n_eu0,
        m_af_b=args.m_af_b,
        m_af_eu=args.m_af_eu,
        admixture_time=args.admixture_time,
        afr_proportion=args.afr_proportion,
        eur_proportion=args.eur_proportion,
        admix_founding_size=args.admix_founding_size,
        admix_initial_size=args.admix_initial_size,
        admix_growth_rate_1=args.admix_growth_rate_1,
        admix_growth_rate_2=args.admix_growth_rate_2,
        admix_growth_rate_2_change_time=(
            args.admix_growth_rate_2_change_time
        ),
        census_time_offset=args.census_time_offset
    )
    sample_sets = build_sample_sets(args.sample_size)
    ts = simulate_3T(
        sample_sets=sample_sets,
        demography=demography,
        contig=contig,
        msprime_model=args.msprime_model,
        seed=args.seed
    )
    ts.dump(args.output)


if __name__ == "__main__":
    main()
