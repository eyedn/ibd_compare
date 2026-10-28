# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# figure.core.R
# ______________________________________________________________________________



safe.ratio <- function(numerator, denominator) {
  result <- numerator / denominator
  result[!is.finite(result)] <- NA_real_
  result
}

safe.cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  cor(x[keep], y[keep], method = method)
}

organize.figure.plots <- function(pairings) {
  result <- list()
  for (pair.name in names(pairings)) {
    for (family in names(pairings[[pair.name]])) {
      plot <- pairings[[pair.name]][[family]]
      attr(plot, "profiles") <- strsplit(pair.name, "_", fixed = TRUE)[[1]]
      result[[family]][[pair.name]] <- plot
    }
  }
  result
}

profile.metadata <- function(directory, code, param.map) {
  if (identical(unname(code), "true")) {
    return(tibble(software = "true", software.with.params = "true",
                  software.code = "true", param.code = "true", imap = "true"))
  }
  params <- sub("^.*_imap.*?_(minseed|d)", "\\1", basename(directory))
  software <- sub("_.*", "", basename(directory))
  full.params <- paste0(software, "_", params)
  profile <- unname(param.map[full.params])
  map <- c(h = "hapmap", d = "decode", f = "1Mb1cM")[substr(code, 2, 2)]
  if (is.na(profile) || is.na(map) ||
      paste0(substr(code, 1, 1), substr(code, 3, 3)) != profile) {
    stop("Unknown parameterization: ", directory)
  }
  tibble(software = software, software.with.params = full.params,
         software.code = unname(code), param.code = profile, imap = unname(map))
}

segment.metadata <- function(files, metadata) {
  parts <- stringr::str_match(basename(files), paste0(
    "^(true|hapibd|rapid)_([^_]+)_([^_]+)_([^_]+)_(\\d+)",
    "(?:\\.([^.]+))?\\.bed\\.gz$"
  ))
  if (anyNA(parts[, 2:6])) stop("Unrecognized segment filename: ",
                                paste(files[is.na(parts[, 2])], collapse = ", "))
  imap <- if (metadata$param.code == "true") "true" else metadata$imap
  if (any(parts[, 2] != metadata$software) ||
      (metadata$param.code != "true" && any(parts[, 4] != imap))) {
    stop("Segment filename disagrees with profile directory: ", files[1])
  }
  tibble(file = files, pop = parts[, 5], rep = as.integer(parts[, 6]),
         filt = ifelse(is.na(parts[, 7]), "none", parts[, 7])) %>%
    bind_cols(metadata[rep(1L, length(files)), ])
}

ks.figure.summary <- function(calls, datasets, after.centromere = FALSE) {
  # Preserve the original comparisons and the population-wide BH family.
  # Truth includes all segments; inferred calls use unrelated pairs.
  inferred <- calls %>% filter(software.code != "true", is.unrel)
  if (after.centromere) inferred <- inferred %>% filter(bp.str > cen.end.bp)
  keys <- c("pop", "rep", "software.code", "param.code", "imap", "filt")
  inferred <- inferred %>% group_by(across(all_of(keys))) %>%
    summarise(inf = list(cM.len), .groups = "drop")
  truth <- calls %>% filter(software.code == "true") %>%
    group_by(pop, rep) %>% summarise(tru = list(cM.len), .groups = "drop")
  datasets %>% filter(software.code != "true") %>%
    distinct(across(all_of(keys))) %>%
    left_join(inferred, by = keys) %>%
    left_join(truth, by = c("pop", "rep")) %>%
    mutate(p.value = map2_dbl(tru, inf, run.ks),
           n.inf = map_int(inf, length), n.true = map_int(tru, length)) %>%
    select(-inf, -tru) %>% group_by(pop) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           log.p.adj = log10(p.adj), is.sig = p.adj <= 0.05,
           filt = factor(filt, levels = filt.levels),
           imap = factor(imap, levels = c("hapmap", "decode", "1Mb1cM")),
           software.code = factor(software.code), param.code = factor(param.code)) %>%
    ungroup()
}

validate.figure.stats <- function(data, profiles, populations, filters,
                                  focal.rep) {
  required <- tidyr::crossing(
    software.code = profiles, pop = populations,
    imap = c("hapmap", "decode", "1Mb1cM"), filt = filters,
    rep = sort(unique(c(data$rep, focal.rep)))
  )
  missing <- anti_join(required, data,
    by = c("software.code", "pop", "imap", "filt", "rep"))
  if (nrow(missing)) stop("Missing required figure statistics: ",
    paste(unlist(missing[1, ]), collapse = "/"))
  invisible(data)
}
