###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_3T_demography.py
###############################################################################


# pattern: Functional Core

import math
import msprime


def build_3T_demography(
    generation_time=25,
    mutation_rate=2.36e-8,
    t_af_years=148000,
    t_ooa_years=51000,
    t_eu0_years=23000,
    t_eg_years=5115,
    r_eu0=0.00307,
    r_eu=0.0195,
    r_af=0.0166,
    n_a=7310,
    n_af1=14474,
    n_b=1861,
    n_eu0=1032,
    m_af_b=15e-5,
    m_af_eu=2.5e-5,
    admixture_time=15.0,
    admix_generation_count=16,
    admix_mixing_generation_count=11,
    admix_ne_by_generation=None,
    admix_afr_props_by_generation=None,
    admix_eur_props_by_generation=None,
    admix_priorAdmix_props_by_generation=None,
    admix_modern_growth_rate=0.023175,
    census_time_offset=1e-6
):
    if admix_ne_by_generation is None:
        admix_ne_by_generation = [
            187.6667, 2.214286, 276.3148, 2019.588,
            10015.88, 45816.03, 149343.0, 381940.1,
            390895.0, 400059.7, 409439.4, 419039.0,
            428863.6, 438918.6, 449209.3, 459741.3
        ]
    if admix_afr_props_by_generation is None:
        admix_afr_props_by_generation = [
            0.8798077, 0.920000, 0.904820, 0.791384,
            0.786692, 0.719992, 0.494960, 0.130456,
            0.060000, 0.060000, 0.060000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000
        ]
    if admix_eur_props_by_generation is None:
        admix_eur_props_by_generation = [
            0.1201923, 0.080000, 0.080000, 0.080000,
            0.080000, 0.080000, 0.080000, 0.080000,
            0.030000, 0.030000, 0.030000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000
        ]
    if admix_priorAdmix_props_by_generation is None:
        admix_priorAdmix_props_by_generation = [
            0.000000, 0.000000, 0.015180, 0.128616,
            0.133308, 0.200008, 0.425040, 0.789544,
            0.910000, 0.910000, 0.910000, 1.000000,
            1.000000, 1.000000, 1.000000, 1.000000
        ]

    # validate the generation-indexed ADMIX schedule
    admix_vectors = [
        ("admix_ne_by_generation", admix_ne_by_generation),
        ("admix_afr_props_by_generation", admix_afr_props_by_generation),
        ("admix_eur_props_by_generation", admix_eur_props_by_generation),
        ("admix_priorAdmix_props_by_generation", admix_priorAdmix_props_by_generation)
    ]
    if admix_generation_count <= 0:
        raise ValueError("admix_generation_count must be positive")
    if admix_mixing_generation_count <= 0:
        raise ValueError("admix_mixing_generation_count must be positive")
    if admix_mixing_generation_count >= admix_generation_count:
        raise ValueError(
            "admix_mixing_generation_count must be less than "
            "admix_generation_count"
        )
    if admixture_time <= 0:
        raise ValueError("admixture_time must be positive")
    if census_time_offset <= 0:
        raise ValueError("census_time_offset must be positive")
    if not math.isclose(
        admixture_time,
        admix_generation_count - 1,
        rel_tol=0,
        abs_tol=1e-9
    ):
        raise ValueError(
            "admixture_time must equal admix_generation_count - 1"
        )
    if admix_modern_growth_rate < 0:
        raise ValueError("admix_modern_growth_rate must be non-negative")

    for vector_name, vector_value in admix_vectors:
        if len(vector_value) != admix_generation_count:
            raise ValueError(
                f"{vector_name} must have length admix_generation_count"
            )
        if not all(math.isfinite(value) for value in vector_value):
            raise ValueError(f"{vector_name} must contain finite values")

    if not all(value > 0 for value in admix_ne_by_generation):
        raise ValueError("admix_ne_by_generation must contain positive values")

    for idx in range(admix_generation_count):
        props = [
            admix_afr_props_by_generation[idx],
            admix_eur_props_by_generation[idx],
            admix_priorAdmix_props_by_generation[idx]
        ]
        if not all(value >= 0 for value in props):
            raise ValueError("admixture proportions must be non-negative")
        if not math.isclose(sum(props), 1.0, rel_tol=0, abs_tol=1e-6):
            raise ValueError(
                "admixture proportions must sum to 1 for every generation"
            )

    if not math.isclose(
        admix_priorAdmix_props_by_generation[0],
        0.0,
        rel_tol=0,
        abs_tol=1e-9
    ):
        raise ValueError("generation 1 priorAdmix proportion must be 0")

    for idx in range(admix_mixing_generation_count, admix_generation_count):
        if not math.isclose(
            admix_priorAdmix_props_by_generation[idx],
            1.0,
            rel_tol=0,
            abs_tol=1e-9
        ):
            raise ValueError("post-mixing generations must have h = 1")

    post_mixing_duration = (
        admix_generation_count - admix_mixing_generation_count
    )
    expected_post_mixing_ne = (
        admix_ne_by_generation[-1]
        * math.exp(-admix_modern_growth_rate * post_mixing_duration)
    )
    observed_post_mixing_ne = (
        admix_ne_by_generation[admix_mixing_generation_count - 1]
    )
    if not math.isclose(
        expected_post_mixing_ne,
        observed_post_mixing_ne,
        rel_tol=1e-5,
        abs_tol=1e-3
    ):
        raise ValueError(
            "ADMIX present-day size and growth rate must imply the "
            "generation-11 Ne at the mixing boundary"
        )

    # convert Tennessen/Fu timings from years to generations
    t_af = t_af_years / generation_time
    t_ooa = t_ooa_years / generation_time
    t_eu0 = t_eu0_years / generation_time
    t_eg = t_eg_years / generation_time

    # reconstruct Tennessen/Fu present-day population sizes
    n_eu1 = n_eu0 / math.exp(-r_eu0 * (t_eu0 - t_eg))
    n_eu = n_eu1 / math.exp(-r_eu * t_eg)
    n_af = n_af1 / math.exp(-r_af * t_eg)

    # define extant source populations and the present-day ADMIX population
    demography = msprime.Demography()
    demography.add_population(
        name="AFR",
        description="African population from the Tennessen _ooa_2 skeleton",
        initial_size=n_af,
        growth_rate=r_af,
        initially_active=True
    )
    demography.add_population(
        name="EUR",
        description="European population from the Tennessen _ooa_2 skeleton",
        initial_size=n_eu,
        growth_rate=r_eu
    )
    demography.add_population(
        name="ADMIX",
        description="African American admixed population",
        initial_size=admix_ne_by_generation[-1],
        growth_rate=admix_modern_growth_rate
    )

    # one-generation hidden ADMIX populations for G1 through G10
    for gen_number in range(1, admix_mixing_generation_count):
        demography.add_population(
            name=f"ADMIX_G{gen_number}",
            description=(
                "One-generation African American admixed population "
                f"for generation {gen_number}"
            ),
            initial_size=admix_ne_by_generation[gen_number - 1],
            growth_rate=0,
            initially_active=False
        )

    # allow recent symmetric migration between AFR and EUR
    demography.set_symmetric_migration_rate(["AFR", "EUR"], m_af_eu)

    # place a census immediately before the first founder event
    demography.add_census(time=admixture_time + census_time_offset)

    # connect present-day ADMIX back to the final mixing generation, G11
    gen_number = admix_mixing_generation_count
    source_name = f"ADMIX_G{gen_number - 1}"
    sources = ["AFR", "EUR", source_name]
    proportions = [
        admix_afr_props_by_generation[gen_number - 1],
        admix_eur_props_by_generation[gen_number - 1],
        admix_priorAdmix_props_by_generation[gen_number - 1]
    ]
    demography.add_admixture(
        time=admix_generation_count - gen_number,
        derived="ADMIX",
        ancestral=sources,
        proportions=proportions
    )

    # walk backward through G10..G2, each with one generation of constant Ne
    for gen_number in range(admix_mixing_generation_count - 1, 1, -1):
        source_name = f"ADMIX_G{gen_number - 1}"
        raw_sources = ["AFR", "EUR", source_name]
        raw_proportions = [
            admix_afr_props_by_generation[gen_number - 1],
            admix_eur_props_by_generation[gen_number - 1],
            admix_priorAdmix_props_by_generation[gen_number - 1]
        ]
        sources = [
            source
            for source, proportion in zip(raw_sources, raw_proportions)
            if proportion > 0
        ]
        proportions = [
            proportion
            for proportion in raw_proportions
            if proportion > 0
        ]
        demography.add_admixture(
            time=admix_generation_count - gen_number,
            derived=f"ADMIX_G{gen_number}",
            ancestral=sources,
            proportions=proportions
        )

    # generation 1 is founded directly from AFR and EUR ancestry
    demography.add_admixture(
        time=admixture_time,
        derived="ADMIX_G1",
        ancestral=["AFR", "EUR"],
        proportions=[
            admix_afr_props_by_generation[0],
            admix_eur_props_by_generation[0]
        ]
    )

    # Tennessen/Fu recent growth epoch boundary
    demography.add_symmetric_migration_rate_change(
        time=t_eg,
        populations=["AFR", "EUR"],
        rate=m_af_eu
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="EUR",
        growth_rate=r_eu0,
        initial_size=n_eu1
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="AFR",
        growth_rate=0,
        initial_size=n_af1
    )

    # older AFR/EUR migration regime and pre-growth EUR size
    demography.add_symmetric_migration_rate_change(
        time=t_eu0,
        populations=["AFR", "EUR"],
        rate=m_af_b
    )
    demography.add_population_parameters_change(
        time=t_eu0,
        population="EUR",
        initial_size=n_b,
        growth_rate=0
    )

    # out-of-Africa split
    demography.add_population_split(
        time=t_ooa,
        derived=["EUR"],
        ancestral="AFR"
    )
    demography.add_migration_rate_change(time=t_ooa, rate=0)

    # ancestral African size before the OOA split
    demography.add_population_parameters_change(
        time=t_af,
        population="AFR",
        initial_size=n_a
    )
    demography.sort_events()

    model_metadata = {
        "id": "OutOfAfrica_AfricanAmericanGeneration3T",
        "description": (
            "Tennessen _ooa_2 skeleton with generation-indexed "
            "African American admixture"
        ),
        "generation_time": generation_time,
        "mutation_rate": mutation_rate,
        "admix_generation_count": admix_generation_count,
        "admix_mixing_generation_count": admix_mixing_generation_count,
        "admix_present_size": admix_ne_by_generation[-1],
        "admix_size_at_mixing_boundary": observed_post_mixing_ne,
        "admix_modern_growth_rate": admix_modern_growth_rate
    }
    return demography, model_metadata
