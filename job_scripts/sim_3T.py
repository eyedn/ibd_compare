###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_3T.py
###############################################################################


import pickle
from pathlib import Path
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
    output_path = Path(args.output)
    output_stem = output_path.with_suffix("")
    demography_pickle_path = output_stem.with_suffix(".demography.pkl")
    model_metadata_pickle_path = output_stem.with_suffix(
        ".model_metadata.pkl"
    )

    species = stdpopsim.get_species("HomSap")
    contig = get_contig(species, args.chromosome, args.genetic_map)
    demography, model_metadata = build_3T_demography(
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
        admix_generation_count=args.admix_generation_count,
        admix_mixing_generation_count=args.admix_mixing_generation_count,
        admix_ne_by_generation=args.admix_ne_by_generation,
        admix_afr_props_by_generation=args.admix_afr_props_by_generation,
        admix_eur_props_by_generation=args.admix_eur_props_by_generation,
        admix_hybrid_props_by_generation=(
            args.admix_hybrid_props_by_generation
        ),
        admix_modern_growth_rate=args.admix_modern_growth_rate,
        census_time_offset=args.census_time_offset
    )
    with demography_pickle_path.open("wb") as demography_handle:
        pickle.dump(demography, demography_handle)
    with model_metadata_pickle_path.open("wb") as model_metadata_handle:
        pickle.dump(model_metadata, model_metadata_handle)
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
