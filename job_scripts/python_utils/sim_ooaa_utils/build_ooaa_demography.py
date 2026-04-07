# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# build_ooaa_demography.py
# ______________________________________________________________________________


import math
import msprime


def build_ooaa_demography(
    admixture_time=14.0,
    afr_proportion=0.75,
    eur_proportion=0.25,
    admix_initial_size=1000.0,
    admix_growth_rate=None,
    census_time_offset=1e-6,
):
    if admixture_time <= 0:
        raise ValueError("admixture_time must be positive")
    if census_time_offset <= 0:
        raise ValueError("census_time_offset must be positive")
    if not math.isclose(afr_proportion + eur_proportion, 1.0, rel_tol=0, abs_tol=1e-9):
        raise ValueError("afr_proportion and eur_proportion must sum to 1")

    generation_time = 25
    mutation_rate = 2.36e-8

    t_af = 148e3 / generation_time
    t_ooa = 51e3 / generation_time
    t_eu0 = 23e3 / generation_time
    t_eg = 5115 / generation_time

    r_eu0 = 0.00307
    r_eu = 0.0195
    r_af = 0.0166

    n_a = 7310
    n_af1 = 14474
    n_b = 1861
    n_eu0 = 1032
    n_eu1 = n_eu0 / math.exp(-r_eu0 * (t_eu0 - t_eg))

    m_af_b = 15e-5
    m_af_eu = 2.5e-5

    n_eu = n_eu1 / math.exp(-r_eu * t_eg)
    n_af = n_af1 / math.exp(-r_af * t_eg)

    if admix_growth_rate is None:
        admix_growth_rate = r_eu

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
        growth_rate=admix_growth_rate,
    )

    demography.set_symmetric_migration_rate(["AFR", "EUR"], m_af_eu)
    demography.add_census(time=admixture_time + census_time_offset)
    demography.add_admixture(
        time=admixture_time,
        derived="ADMIX",
        ancestral=["AFR", "EUR"],
        proportions=[afr_proportion, eur_proportion],
    )

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
    demography.add_population_split(
        time=t_ooa,
        derived=["EUR"],
        ancestral="AFR",
    )
    demography.add_migration_rate_change(time=t_ooa, rate=0)
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
