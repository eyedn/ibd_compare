# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# filter_and_format.R
# ______________________________________________________________________________


library(tidyverse)

# source function ----
source("job_scripts/R_utils/read.hapibd.R")
source("job_scripts/R_utils/read.rapid.R")
source("job_scripts/R_utils/read.true.R")
source("job_scripts/R_utils/filter.ibd.by.hot.bins.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
  stop("incorrect number of parameters")
}

# inputs ----
min.bp.size <- as.numeric(args[1])
min.cm.size <- as.numeric(args[2])
hot.bin.size <- as.numeric(args[3])
hot.bin.sd <- as.numeric(args[4])
ibd.file <- args[5]
gmap.file <- args[6]
ff.out <- args[7]
ff.hot.out <- args[8]
is.hapibd.format <- tolower(args[9])
do.hot.filter <- tolower(args[10])
if (!is.hapibd.format %in% c("true", "false")) {
  stop("is.hapibd.format must be TRUE/true or FALSE/false; got: ", args[9])
}
is.hapibd.format <- is.hapibd.format == "true"
if (!do.hot.filter %in% c("true", "false")) {
  stop("do.hot.filter must be TRUE/true or FALSE/false; got: ", args[10])
}
do.hot.filter <- do.hot.filter == "true"

print("inputs:")
print(paste("min.bp.size = ", min.bp.size))
print(paste("min.cm.size = ", min.cm.size))
print(paste("ibd.file = ", ibd.file))
print(paste("gmap.file = ", gmap.file))
print(paste("ff.out = ", ff.out))
print(paste("ff.hot.out = ", ff.hot.out))
print(paste("is.hapibd.format = ", is.hapibd.format))
print(paste("do.hot.filter = ", do.hot.filter))

empty.calls <- tibble::tibble(
  chr = character(), bp.str = integer(), bp.end = integer(), bp.len = integer(),
  cM.str = numeric(), cM.end = numeric(), cM.len = numeric(),
  sam1 = character(), hap1 = integer(), sam2 = character(), hap2 = integer()
)
if (length(readr::read_lines(ibd.file)) == 0) {
  readr::write_tsv(empty.calls, ff.out, col_names = TRUE)
  if (do.hot.filter) {
    readr::write_tsv(empty.calls, ff.hot.out, col_names = TRUE)
  }
  quit(status = 0)
}

# read in ibd calls ----
if (grepl("true_ibd", ibd.file, ignore.case = TRUE)) {
  segs <- read.true(ibd.file, min.bp.size, min.cm.size, gmap.file)
  
  # print true segments in BED format
  readr::write_tsv(segs, ff.out, col_names = TRUE)
  quit(status = 0)
} else if (is.hapibd.format) {
  calls <- read.hapibd.bed(
    ibd.file,
    min.bp.size,
    min.cm.size,
    gmap.file
  )
  gmap <- readr::read_table(
    file = gmap.file, col_names = FALSE, show_col_types = FALSE
  ) %>%
    rename(chr = X1, snp = X2, cM = X3, bp = X4)
} else {
  calls <- read.rapid.bed(
    ibd.file,
    min.bp.size,
    min.cm.size,
    gmap.file
  )
  gmap <- readr::read_table(
    file = gmap.file, col_names = FALSE, show_col_types = FALSE
  ) %>%
    rename(idx = X1, cM = X2)
}

# print ibd segments in BED format ----
readr::write_tsv(calls, ff.out, col_names = TRUE)

if (do.hot.filter) {
  # filter out hot bins (bins with too great of a number of segemnts) ----
  chr.len.cm <- max(as.numeric(gmap$cM), na.rm = TRUE)
  hot.all <- filter.ibd.by.hot.bins(
    calls = calls, chr.len.cm = chr.len.cm, bin.size.cM = hot.bin.size, 
    sd.thresh = hot.bin.sd, trim.frac = 0.03
    )
  
  # print ibd segments filter by hot bins in BED format ----
  readr::write_tsv(hot.all$calls.filtered, ff.hot.out, col_names = TRUE)
}
