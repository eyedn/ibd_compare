###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_args.py
###############################################################################


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-size", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--msprime-model", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--genetic-map", required=True)
    parser.add_argument("--generation-time", type=float, default=25)
    parser.add_argument("--mutation-rate", type=float, default=2.36e-8)
    parser.add_argument("--t-af-years", type=float, default=148000)
    parser.add_argument("--t-ooa-years", type=float, default=51000)
    parser.add_argument("--t-eu0-years", type=float, default=23000)
    parser.add_argument("--t-eg-years", type=float, default=5115)
    parser.add_argument("--r-eu0", type=float, default=0.00307)
    parser.add_argument("--r-eu", type=float, default=0.0195)
    parser.add_argument("--r-af", type=float, default=0.0166)
    parser.add_argument("--n-a", type=float, default=7310)
    parser.add_argument("--n-af1", type=float, default=14474)
    parser.add_argument("--n-b", type=float, default=1861)
    parser.add_argument("--n-eu0", type=float, default=1032)
    parser.add_argument("--m-af-b", type=float, default=15e-5)
    parser.add_argument("--m-af-eu", type=float, default=2.5e-5)
    parser.add_argument("--admixture-time", type=float, default=15.0)
    parser.add_argument(
        "--admix-generation-count",
        type=int,
        default=16
    )
    parser.add_argument(
        "--admix-mixing-generation-count",
        type=int,
        default=11
    )
    parser.add_argument(
        "--admix-ne-by-generation",
        type=lambda value: [float(item) for item in value.split(",")],
        default=[
            187.6667, 2.214286, 276.3148, 2019.588,
            10015.88, 45816.03, 149343.0, 381940.1,
            390895.0, 400059.7, 409439.4, 419039.0,
            428863.6, 438918.6, 449209.3, 459741.3
        ]
    )
    parser.add_argument(
        "--admix-afr-props-by-generation",
        type=lambda value: [float(item) for item in value.split(",")],
        default=[
            0.8798077, 0.920000, 0.904820, 0.791384,
            0.786692, 0.719992, 0.494960, 0.130456,
            0.060000, 0.060000, 0.060000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000
        ]
    )
    parser.add_argument(
        "--admix-eur-props-by-generation",
        type=lambda value: [float(item) for item in value.split(",")],
        default=[
            0.1201923, 0.080000, 0.080000, 0.080000,
            0.080000, 0.080000, 0.080000, 0.080000,
            0.030000, 0.030000, 0.030000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000
        ]
    )
    parser.add_argument(
        "--admix-priorAdmix-props-by-generation",
        type=lambda value: [float(item) for item in value.split(",")],
        dest="admix_priorAdmix_props_by_generation",
        default=[
            0.000000, 0.000000, 0.015180, 0.128616,
            0.133308, 0.200008, 0.425040, 0.789544,
            0.910000, 0.910000, 0.910000, 1.000000,
            1.000000, 1.000000, 1.000000, 1.000000
        ]
    )
    parser.add_argument(
        "--admix-modern-growth-rate",
        type=float,
        default=0.023175
    )
    parser.add_argument("--census-time-offset", type=float, default=1e-6)
    return parser.parse_args()
