# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# infer.software.with.params.R
# ______________________________________________________________________________


library(tidyverse)

infer.software.with.params <- function(software.base, software.dir.base) {
    # hapibd directory example:
    # hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2
    if (software.base == "hapibd") {
    m <- stringr::str_match(
        software.dir.base,
        "minseed[^_]+_gap[^_]+_extend[^_]+_out[^_]+_mark[^_]+_mac[^_]+"
    )
    if (!is.na(m[1, 1])) {
        return(paste0("hapibd_", m[1, 1]))
    } else {
        return("hapibd")
    }
    }

    # rapid directory example:
    # rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w300_r10_s2
    if (software.base == "rapid") {
    m <- stringr::str_match(
        software.dir.base,
        "d[^_]+_w[^_]+_r[^_]+_s[^_]+"
    )
    if (!is.na(m[1, 1])) {
        return(paste0("rapid_", m[1, 1]))
    } else {
        return("rapid")
    }
    }

    # fallback for any future software
    software.base
}
