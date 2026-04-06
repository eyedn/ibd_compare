###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           plot.compare.ibd.by.bins.R
###############################################################################

library(tidyverse)
library(data.table)
library(RColorBrewer)

plot.compare.ibd.by.bins <- function(
    df,
    pop,
    tmap,              # single true.map; must match a column in chr.sizes
    rep,               # single replicate
    unrel.status,      # "all" or "unrel"
    imaps,             # vector of infer.map values
    methods,           # vector of method values, same length as imaps
    use.bp.bin = TRUE, # TRUE = bin in bp, FALSE = bin in cM
    bin.size = 1,      # units: bp or cM depending on use.bp.bin
    chr.sizes,        # 2 x K matrix/df; row1 = cM, row2 = bp
    chr.label,        # e.g., chr1, chr18
    title = NULL,
    x.start = NULL,    # custom x-axis lower bound (same units as binning)
    x.end   = NULL     # custom x-axis upper bound (same units as binning)
) {
  # basic checks
  if (length(imaps) != length(methods)) {
    stop("imaps and methods must be the same length.")
  }
  if (!tmap %in% colnames(chr.sizes)) {
    stop("true.map '", tmap, "' not found as a column in chr.sizes")
  }
  if (!is.null(x.start) || !is.null(x.end)) {
    if (is.null(x.start) || is.null(x.end)) {
      stop("If you provide a custom x-axis interval, both x.start and x.end must be non-NULL.")
    }
    if (!is.numeric(x.start) || !is.numeric(x.end)) {
      stop("x.start and x.end must be numeric.")
    }
    if (x.start >= x.end) {
      stop("x.start must be strictly less than x.end.")
    }
  }
  
  # define the (method, infer.map) combos we want to include
  subsets <- tibble(
    method    = methods,
    infer.map = imaps
  ) %>%
    distinct()
  
  # filter to specific pop, true.map, rep, unrel.status
  # then keep only the desired (method, infer.map) combos
  df.sub <- df %>%
    filter(
      pop          == !!pop,
      true.map     == !!tmap,
      rep          == !!rep,
      unrel.status == !!unrel.status
    ) %>%
    semi_join(subsets, by = c("method", "infer.map"))
  
  if (!nrow(df.sub)) {
    stop("No rows left after filtering by pop, tmap, rep, unrel.status, and (method, infer.map) subsets.")
  }
  
  # construct a sample-haplotype pair ID
  df.sub <- df.sub %>%
    mutate(
      sam.hap.pair.id = sprintf("%d_%d_%d_%d", sam1, hap1, sam2, hap2)
    )
  
  # row order = union of present pairs after filters
  pair.levels <- df.sub %>%
    distinct(sam.hap.pair.id, sam1, hap1, sam2, hap2) %>%
    arrange(sam1, sam2, hap1, hap2) %>%
    pull(sam.hap.pair.id)
  
  if (!length(pair.levels)) {
    stop("No sample-haplotype pairs after filtering.")
  }
  
  ## ---- DYNAMIC OR CUSTOM X-RANGE ------------------------------------------
  if (!is.null(x.start) && !is.null(x.end)) {
    # Use user-specified interval directly
    seg.min <- x.start
    seg.max <- x.end
  } else {
    # Fall back to range of observed segments
    if (use.bp.bin) {
      seg.min <- min(df.sub$bp.str, na.rm = TRUE)
      seg.max <- max(df.sub$bp.end, na.rm = TRUE)
    } else {
      seg.min <- min(df.sub$cM.str, na.rm = TRUE)
      seg.max <- max(df.sub$cM.end, na.rm = TRUE)
    }
  }
  
  # snap to bin boundaries
  cmin <- floor(seg.min / bin.size) * bin.size
  cmax <- ceiling(seg.max / bin.size) * bin.size
  
  n.bin <- ceiling((cmax - cmin) / bin.size)
  bins <- tibble(
    chr        = chr.label,
    bin.idx    = seq_len(n.bin),
    str        = cmin + (bin.idx - 1) * bin.size,
    end        = pmin(cmin + bin.idx * bin.size, cmax),
    global.bin = bin.idx
  )
  ## -------------------------------------------------------------------------
  
  # build bit mask for unique (method, infer.map) pairs
  combo.bits <- df.sub %>%
    distinct(method, infer.map) %>%
    mutate(
      combo = paste(method, infer.map, sep = "_"),
      bit   = bitwShiftL(1L, row_number() - 1L)
    )
  
  # scan for interval overlaps (segments × bins)
  suppressWarnings({
    bins.dt <- as.data.table(bins)
    setnames(bins.dt, c("str", "end"), c("start", "end"))
    setkey(bins.dt, chr, start, end)
    
    segs.dt <- df.sub %>%
      transmute(
        chr = as.character(
          if_else(grepl("^chr", chr), chr, paste0("chr", chr))
        ),
        start = if (use.bp.bin) bp.str else cM.str,
        end   = if (use.bp.bin) bp.end else cM.end,
        sam.hap.pair.id,
        method,
        infer.map
      ) %>%
      as.data.table()
    setkey(segs.dt, chr, start, end)
  })
  
  hits <- foverlaps(segs.dt, bins.dt, type = "any", nomatch = 0L)
  
  hits.tbl <- as_tibble(hits) %>%
    select(sam.hap.pair.id, method, infer.map, global.bin) %>%
    distinct() %>%
    inner_join(combo.bits, by = c("method", "infer.map"))
  
  # OR-combine bits within each (pair, bin)
  codes <- hits.tbl %>%
    group_by(sam.hap.pair.id, global.bin) %>%
    summarise(code = Reduce(bitwOr, bit), .groups = "drop")
  
  # dense n × m table of integer codes (pairs × bins)
  base <- expand_grid(
    sam.hap.pair.id = pair.levels,
    global.bin      = bins$global.bin
  ) %>%
    left_join(codes, by = c("sam.hap.pair.id", "global.bin")) %>%
    mutate(code = coalesce(code, 0L))
  
  mat <- base %>%
    pivot_wider(
      names_from  = global.bin,
      values_from = code,
      values_fill = 0L
    ) %>%
    separate(
      col  = sam.hap.pair.id,
      into = c("sam1", "hap1", "sam2", "hap2"),
      sep  = "_",
      remove = FALSE
    ) %>%
    mutate(
      sam1 = as.numeric(sam1), hap1 = as.numeric(hap1),
      sam2 = as.numeric(sam2), hap2 = as.numeric(hap2)
    ) %>%
    arrange(sam1, sam2, hap1, hap2) %>%
    select(-sam1, -hap1, -sam2, -hap2)
  
  mat.out <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat.out) <- mat$sam.hap.pair.id
  colnames(mat.out) <- paste0("b", bins$global.bin)
  
  # legend for codes (now in terms of full method/map combos)
  uniq.codes <- sort(unique(as.integer(mat.out)))
  uniq.codes <- uniq.codes[uniq.codes > 0L]
  
  decode.combo <- function(code) {
    present <- combo.bits$combo[bitwAnd(code, combo.bits$bit) > 0]
    paste(present, collapse = " + ")
  }
  
  legend.tbl <- tibble(
    code  = c(0L, uniq.codes),
    combo = c("None", vapply(uniq.codes, decode.combo, character(1)))
  )
  
  base.plot.df <- base %>%
    mutate(
      y      = factor(sam.hap.pair.id, levels = rev(pair.levels)),
      code.f = factor(code, levels = legend.tbl$code, labels = legend.tbl$combo)
    )
  
  if (is.null(title)) {
    title <- paste0(
      "pop=", pop,
      " | true.map=", tmap,
      " | rep=", rep,
      " | unrel.status=", unrel.status,
      " | ", if (use.bp.bin) "bp bins" else "cM bins",
      " (bin size=", bin.size, ")"
    )
  }
  
  xlab <- if (use.bp.bin) {
    paste0("chr bins (", format(cmin, big.mark = ","), "–",
           format(cmax, big.mark = ","), " bp)")
  } else {
    paste0("chr bins (", signif(cmin, 4), "–", signif(cmax, 4), " cM)")
  }
  
  ### ---- FORCE "None" TO ALWAYS BE BLACK -----------------------------------
  other.labels <- legend.tbl$combo[legend.tbl$combo != "None"]
  other.colors <- brewer.pal(max(3, length(other.labels)), "Set2")[seq_along(other.labels)]
  
  fill.colors <- c("None" = "black")
  if (length(other.labels) > 0) {
    names(other.colors) <- other.labels
    fill.colors <- c(fill.colors, other.colors)
  }
  ### ------------------------------------------------------------------------
  
  p <- ggplot(
    base.plot.df,
    aes(x = global.bin, y = y, fill = code.f)
  ) +
    geom_tile() +
    scale_fill_manual(
      values = fill.colors,
      drop = FALSE
    ) +
    labs(
      x     = xlab,
      y     = "sample-haplotype pair",
      fill  = "method_map",
      title = title
    ) +
    theme_bw()
  
  list(
    matrix = mat.out,
    plot   = p,
    legend = legend.tbl,
    bins   = bins,
    pairs  = pair.levels,
    bitkey = combo.bits
  )
}
