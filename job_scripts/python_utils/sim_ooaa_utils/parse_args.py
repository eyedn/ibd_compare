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
    parser.add_argument("--admixture-time", type=float, default=14.0)
    parser.add_argument("--afr-proportion", type=float, default=0.75)
    parser.add_argument("--eur-proportion", type=float, default=0.25)
    parser.add_argument("--admix-initial-size", type=float, default=1000.0)
    parser.add_argument(
        "--admix-growth-rate-1",
        type=lambda value: None if value == "" else float(value),
        default=0.03
    )
    parser.add_argument(
        "--admix-growth-rate-2",
        type=float,
        default=0.01805
    )
    parser.add_argument(
        "--admix-growth-rate-2-change-time",
        type=float,
        default=7.0
    )
    parser.add_argument("--census-time-offset", type=float, default=1e-6)
    return parser.parse_args()
