# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# parse_args.py
# ______________________________________________________________________________


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-size", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--msprime-model", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--genetic-map", required=True)
    parser.add_argument("--admixture-time", type=float, default=14.0)
    parser.add_argument("--afr-proportion", type=float, default=0.75)
    parser.add_argument("--eur-proportion", type=float, default=0.25)
    parser.add_argument("--admix-initial-size", type=float, default=1000.0)
    parser.add_argument("--admix-growth-rate", type=float, default=None)
    parser.add_argument("--census-time-offset", type=float, default=1e-6)
    return parser.parse_args()
