# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# write_local_ancestry.py
# ______________________________________________________________________________


POP_ID_TO_NAME = {
    0: "AFR",
    1: "EUR",
    2: "ADMIX",
}


def write_local_ancestry(ancestry_table, sample_node_rows, out_tsv, chr_used):
    sample_to_ind = {node: ind_id for node, ind_id, _ in sample_node_rows}
    sample_to_hap = {node: hap for node, _, hap in sample_node_rows}

    ancestry_table["sample_id"] = ancestry_table["sample"].map(
        lambda u: f"tsk_{sample_to_ind[int(u)]}"
    )
    ancestry_table["hap"] = ancestry_table["sample"].map(
        lambda u: sample_to_hap[int(u)]
    )
    ancestry_table["population_name"] = ancestry_table["population"].map(
        lambda p: POP_ID_TO_NAME.get(int(p), f"pop_{int(p)}")
    )
    ancestry_table.insert(0, "chrom", f"chr{chr_used}")
    ancestry_table = ancestry_table[
        ["chrom", "left", "right", "sample_id", "hap", "population_name", "population", "ancestor"]
    ]
    ancestry_table.to_csv(out_tsv, sep="\t", header=False, index=False)
