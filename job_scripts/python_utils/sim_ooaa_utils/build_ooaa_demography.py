###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_ooaa_demography.py
###############################################################################


import math
import msprime


def build_ooaa_demography(
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
    admixture_time=14.0,
    afr_proportion=0.75,
    eur_proportion=0.25,
    admix_initial_size=1000.0,
    admix_growth_rate_1=0.03,
    admix_growth_rate_2=0.01805,
    admix_growth_rate_2_change_time=7.0,
    census_time_offset=1e-6,
):
    # validate user-specified admixture parameters
    if admixture_time <= 0:
        raise ValueError("admixture_time must be positive")
    if census_time_offset <= 0:
        raise ValueError("census_time_offset must be positive")
    if admix_growth_rate_2_change_time <= 0:
        raise ValueError("admix_growth_rate_2_change_time must be positive")
    if admix_growth_rate_2_change_time >= admixture_time:
        raise ValueError(
            "admix_growth_rate_2_change_time must be less than "
            "admixture_time"
        )
    if not math.isclose(
        afr_proportion + eur_proportion,
        1.0,
        rel_tol=0,
        abs_tol=1e-9,
    ):
        raise ValueError("afr_proportion and eur_proportion must sum to 1")

    # convert Tennessen/Fu timings from years to generations
    t_af = t_af_years / generation_time
    t_ooa = t_ooa_years / generation_time
    t_eu0 = t_eu0_years / generation_time
    t_eg = t_eg_years / generation_time

    # reconstruct Tennessen/Fu present-day population sizes
    n_eu1 = n_eu0 / math.exp(-r_eu0 * (t_eu0 - t_eg))
    n_eu = n_eu1 / math.exp(-r_eu * t_eg)
    n_af = n_af1 / math.exp(-r_af * t_eg)

    # a note on on the paramters for the admix population:
    # 1. the population growth rate set at population creation applies from
    #    the present back to admix_growth_rate_2_change_time
    # 2. the change at admix_growth_rate_2_change_time applies farther back in
    #    time, from admix_growth_rate_2_change_time to admixture_time
    # therefore:
    # - admix_growth_rate_2 for the recent interval
    # - admix_growth_rate_1 for the older post-admixture interval
    #
    # the size update keeps ADMIX continuous at the rate-change boundary
    admix_size_at_growth_rate_2_change = admix_initial_size / math.exp(
        -admix_growth_rate_2 * admix_growth_rate_2_change_time
    )

    # define the three extant populations at the present day
    demography = msprime.Demography()
    demography.add_population(
        name="AFR",
        description="African population from the Tennessen _ooa_2 skeleton",
        initial_size=n_af,
        growth_rate=r_af,
    )
    demography.add_population(
        name="EUR",
        description="European population from the Tennessen _ooa_2 skeleton",
        initial_size=n_eu,
        growth_rate=r_eu,
    )
    demography.add_population(
        name="ADMIX",
        description="African American admixed population",
        initial_size=admix_initial_size,
        growth_rate=admix_growth_rate_2,
    )

    # allow recent symmetric migration between AFR and EUR
    demography.set_symmetric_migration_rate(["AFR", "EUR"], m_af_eu)

    # place a census immediately before the admixture event for ancestry labels
    demography.add_census(time=admixture_time + census_time_offset)

    # move from the recent ADMIX growth phase to the older post-admixture phase
    demography.add_population_parameters_change(
        time=admix_growth_rate_2_change_time,
        population="ADMIX",
        initial_size=admix_size_at_growth_rate_2_change,
        growth_rate=admix_growth_rate_1,
    )

    # create ADMIX by mixing ancestry from AFR and EUR
    demography.add_admixture(
        time=admixture_time,
        derived="ADMIX",
        ancestral=["AFR", "EUR"],
        proportions=[afr_proportion, eur_proportion],
    )

    # Tennessen/Fu recent growth epoch boundary
    demography.add_symmetric_migration_rate_change(
        time=t_eg,
        populations=["AFR", "EUR"],
        rate=m_af_eu,
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="EUR",
        growth_rate=r_eu0,
        initial_size=n_eu1,
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="AFR",
        growth_rate=0,
        initial_size=n_af1,
    )

    # older AFR/EUR migration regime and pre-growth EUR size
    demography.add_symmetric_migration_rate_change(
        time=t_eu0,
        populations=["AFR", "EUR"],
        rate=m_af_b,
    )
    demography.add_population_parameters_change(
        time=t_eu0,
        population="EUR",
        initial_size=n_b,
        growth_rate=0,
    )

    # out-of-Africa split
    demography.add_population_split(
        time=t_ooa,
        derived=["EUR"],
        ancestral="AFR",
    )
    demography.add_migration_rate_change(time=t_ooa, rate=0)

    # ancestral African size before the OOA split
    demography.add_population_parameters_change(
        time=t_af,
        population="AFR",
        initial_size=n_a,
    )
    demography.sort_events()

    model_metadata = {
        "id": "OutOfAfrica_AfricanAmericanSkeleton",
        "description": "Tennessen _ooa_2 skeleton with a recent AFR/EUR admixture event",
        "generation_time": generation_time,
        "mutation_rate": mutation_rate,
    }
    return demography, model_metadata
