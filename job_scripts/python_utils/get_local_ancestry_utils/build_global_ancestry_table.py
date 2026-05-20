###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_global_ancestry_table.py
###############################################################################


POP_ID_TO_NAME = {
    0: "AFR",
    1: "EUR",
}


def build_global_ancestry_table(ancestry_table, sample_node_rows):
    sample_to_ind = {node: ind_id for node, ind_id, _ in sample_node_rows}

    ancestry_table = ancestry_table.copy()
    ancestry_table["sample_ind"] = ancestry_table["sample"].map(
        lambda node: sample_to_ind[int(node)]
    )
    ancestry_table["sample_id"] = ancestry_table["sample"].map(
        lambda node: f"tsk_{sample_to_ind[int(node)]}"
    )
    ancestry_table["population_name"] = ancestry_table["population"].map(
        lambda pop: POP_ID_TO_NAME.get(int(pop), f"pop_{int(pop)}")
    )
    ancestry_table["span"] = (
        ancestry_table["right"] - ancestry_table["left"]
    )

    ancestry_by_pop = (
        ancestry_table.groupby(
            ["sample_ind", "sample_id", "population_name"], as_index=False
        )["span"].sum()
    )
    ancestry_wide = ancestry_by_pop.pivot(
        index=["sample_ind", "sample_id"],
        columns="population_name",
        values="span"
    ).fillna(0.0)

    for pop_name in ("AFR", "EUR"):
        if pop_name not in ancestry_wide.columns:
            ancestry_wide[pop_name] = 0.0

    ancestry_wide = ancestry_wide[["AFR", "EUR"]].copy()
    total_span = ancestry_wide.sum(axis=1)
    ancestry_wide["AFR_prop"] = ancestry_wide["AFR"] / total_span
    ancestry_wide["EUR_prop"] = ancestry_wide["EUR"] / total_span
    ancestry_wide = ancestry_wide.reset_index().sort_values("sample_ind")

    return ancestry_wide[["sample_id", "AFR_prop", "EUR_prop"]]
