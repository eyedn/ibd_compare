# ---------------------------------------------------------------------------
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# filter.ibd.by.hot.bins.R
# ---------------------------------------------------------------------------

filter.ibd.by.hot.bins <- function(
    calls,
    chr.len.cM,  
    bin.size.cM = 0.25,      
    sd.thresh   = 10,
    trim.frac  = 0.03
) {

    cat("---- hot-bin IBD filter (trimmed stats; full chr) ----\n")

    chr <- unique(calls$chr)
    if (length(chr) != 1) stop("Expected one chromosome in calls")
    chr <- chr[[1]]

    req.cols <- c("chr", "cM.str", "cM.end")
    miss <- setdiff(req.cols, names(calls))
    if (length(miss) > 0) {
        stop("calls missing required columns: ", paste(miss, collapse = ", "))
    }
    if (!is.numeric(chr.len.cM) || length(chr.len.cM) != 1 || !is.finite(chr.len.cM) || chr.len.cM <= 0) {
        stop("chr.len.cM must be a single positive finite number")
    }

    cat("chromosome: ", chr, "\n")
    cat("init. calls: ", nrow(calls), "\n")

    # normalize cM coordinates
    calls.norm <- calls %>%
        mutate(
            cM.str = as.numeric(cM.str),
            cM.end = as.numeric(cM.end),
            cM.lo  = pmin(cM.str, cM.end),
            cM.hi  = pmax(cM.str, cM.end)
        )

    # total bins across FULL chromosome (include empty bins)
    n.bins <- ceiling(chr.len.cM / bin.size.cM)
    cat("bin size (cM): ", bin.size.cM, "\n")
    cat("chr.len.cM: ", chr.len.cM, "\n")
    cat("total bins (full chr): ", n.bins, "\n")

    # map segments to bin indices (clamp to [1, n.bins])
    # treat bins as [0,bin), [bin,2bin), ... by using -eps on hi
    eps <- 1e-12
    bin.lo <- floor(calls.norm$cM.lo / bin.size.cM) + 1L
    bin.hi <- floor((calls.norm$cM.hi - eps) / bin.size.cM) + 1L

    bin.lo <- pmax(bin.lo, 1L)
    bin.hi <- pmax(bin.hi, bin.lo)
    bin.lo <- pmin(bin.lo, n.bins)
    bin.hi <- pmin(bin.hi, n.bins)

    # count overlaps per bin via difference array
    diff.vec <- integer(n.bins + 2L)
    if (nrow(calls.norm) > 0) {
        for (i in seq_along(bin.lo)) {
            diff.vec[bin.lo[i]] <- diff.vec[bin.lo[i]] + 1L
            diff.vec[bin.hi[i] + 1L] <- diff.vec[bin.hi[i] + 1L] - 1L
        }
    }
    counts <- cumsum(diff.vec)[seq_len(n.bins)]

    bin.counts <- tibble(
        chr = chr,
        bin = seq_len(n.bins),
        count = counts,
        cM.bin.str = (seq_len(n.bins) - 1) * bin.size.cM,
        cM.bin.end = pmin(seq_len(n.bins) * bin.size.cM, chr.len.cM)
    )


    # quick summaries (bin counts + segment lengths)
    seg.len <- (calls.norm$cM.hi - calls.norm$cM.lo)
    cat("\n--- Bin overlap counts (segments overlapping each bin) ---\n")
    print(summary(bin.counts$count))  # Min, 1st Qu., Median, Mean, 3rd Qu., Max
    cat("sd(bin counts): ", signif(sd(bin.counts$count), 6), "\n")
    cat("total overlap mass (sum of bin counts): ", sum(bin.counts$count), "\n")

    cat("\n--- Segment length in cM (cM.hi - cM.lo) ---\n")
    print(summary(seg.len))           # Min, 1st Qu., Median, Mean, 3rd Qu., Max
    cat("sd(seg len cM): ", signif(sd(seg.len), 6), "\n")
    cat("n segments: ", length(seg.len), "\n\n")

    # trimmed mean + trimmed SD (3% trimmed) 
    x <- bin.counts$count
    n <- length(x)
    k <- floor(trim.frac * n)

    if (n < 2 || (n - 2*k) < 2) {
        cat("WARNING: too few bins after trimming; no bins flagged.\n")
        cutoff <- Inf
        bin.counts <- bin.counts %>%
            mutate(
                trim.frac = trim.frac,
                mean.trim = NA_real_,
                sd.trim   = NA_real_,
                cutoff    = cutoff,
                flagged   = FALSE
            )
        flagged.bins <- integer(0)
    } else {
        x.sort <- sort(x)
        x.trim <- x.sort[(k + 1):(n - k)]   # drop lowest & highest 3%
        mean.trim <- mean(x.trim)
        sd.trim   <- sd(x.trim)

        cat("trim.frac: ", trim.frac, "\n")
        cat("trimmed mean segments/bin: ", signif(mean.trim, 4), "\n")
        cat("trimmed sd segments/bin: ", signif(sd.trim, 4), "\n")
        max_z <- (max(bin.counts$count) - mean.trim) / sd.trim
        cat("max trimmed z-score across bins: ", signif(max_z, 4), "\n")

        if (!is.finite(sd.trim) || sd.trim <= 0) {
            cat("WARNING: insufficient trimmed variance; no bins flagged.\n")
            cutoff <- Inf
            bin.counts <- bin.counts %>%
                mutate(
                    trim.frac = trim.frac,
                    mean.trim = mean.trim,
                    sd.trim   = sd.trim,
                    cutoff    = cutoff,
                    flagged   = FALSE
                )
            flagged.bins <- integer(0)
        } else {
            cutoff <- mean.trim + sd.thresh * sd.trim
            cat(
                "Hot-bin cutoff (trimmed mean + ", sd.thresh, " trimmed SD): ",
                signif(cutoff, 4), "\n", sep = ""
            )
            bin.counts <- bin.counts %>%
                mutate(
                    trim.frac = trim.frac,
                    mean.trim = mean.trim,
                    sd.trim   = sd.trim,
                    cutoff    = cutoff,
                    flagged   = count > cutoff
                )
            flagged.bins <- bin.counts %>% filter(flagged) %>% pull(bin)
        }
    }

    cat("flagged bins: ", length(flagged.bins), "\n")

    # remove segments overlapping flagged bins
    if (length(flagged.bins) == 0) {
        cat("No hot bins detected; calls unchanged.\n")
        calls.filtered <- calls %>% mutate(cM.str = as.numeric(cM.str), cM.end = as.numeric(cM.end))
    } else {
        flag.vec <- integer(n.bins + 1L)
        flag.vec[flagged.bins] <- 1L
        flag.cum <- cumsum(flag.vec)

        overlap.flagged <-
            (flag.cum[bin.hi] - flag.cum[pmax(bin.lo - 1L, 0L)]) > 0L

        calls.filtered <- calls.norm %>%
            filter(!overlap.flagged) %>%
            select(-cM.lo, -cM.hi)

        cat("Calls removed: ", sum(overlap.flagged), "\n")
    }

    cat("calls after filtering: ", nrow(calls.filtered), "\n")
    cat("---- Done ----\n")

    list(
        calls.filtered = calls.filtered,
        flagged.bins   = bin.counts %>% filter(flagged),
        bin.counts     = bin.counts
    )
}