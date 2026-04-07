# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# get_contig.py
# ______________________________________________________________________________


def get_contig(species, chromosome, genetic_map):
    contig_name = f"chr{chromosome}"
    if "DefaultUniform" in genetic_map:
        return species.get_contig(contig_name)
    try:
        return species.get_contig(contig_name, genetic_map=genetic_map)
    except TypeError:
        return species.get_contig(
            contig_name,
            genetic_map=species.get_genetic_map(genetic_map),
        )
