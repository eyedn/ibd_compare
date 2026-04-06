# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# filter_and_format.R
# ______________________________________________________________________________

# libraries ----
library(tidyverse)

# source function ----
source("job_scripts/R_utils/read.hapibd.R")
source("job_scripts/R_utils/read.rapid.R")
source("job_scripts/R_utils/read.true.R")
source("job_scripts/R_utils/filter.ibd.by.hot.bins.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("incorrect number of parameters")
}

# inputs ----
MIN.BP.SIZE <- as.numeric(args[1])
MIN.CM.SIZE <- as.numeric(args[2])
HOT.BIN.SIZE <- as.numeric(args[3])
HOT.BIN.SD <- as.numeric(args[4])
IBD.FILE <- args[5]
GMAP <- args[6]
FF.OUT <- args[7]
FF.HOT.OUT <- args[8]
IS.HAPIBD.FORMAT <- tolower(args[9])
if (!IS.HAPIBD.FORMAT %in% c("true", "false")) {
  stop("IS.HAPIBD.FORMAT must be TRUE/true or FALSE/false; got: ", args[9])
}
IS.HAPIBD.FORMAT <- IS.HAPIBD.FORMAT == "true"

print("inputs:")
print(paste("MIN.BP.SIZE = ", MIN.BP.SIZE))
print(paste("MIN.CM.SIZE = ", MIN.CM.SIZE))
print(paste("IBD.FILE = ", IBD.FILE))
print(paste("GMAP = ", GMAP))
print(paste("FF.OUT = ", FF.OUT))
print(paste("FF.HOT.OUT = ", FF.HOT.OUT))
print(paste("IS.HAPIBD.FORMAT = ", IS.HAPIBD.FORMAT))

# read in ibd calls ----
if (grepl("true_ibd", IBD.FILE, ignore.case = TRUE)) {
  segs <- read.true(IBD.FILE, MIN.BP.SIZE, MIN.CM.SIZE, GMAP)
  
  # print true segments in BED format
  readr::write_tsv(segs, FF.OUT, col_names = TRUE)
  quit(status = 0)
} else if (IS.HAPIBD.FORMAT) {
  calls <- read.hapibd.bed(IBD.FILE, MIN.BP.SIZE, MIN.CM.SIZE, GMAP)
  gmap <- readr::read_table(
    file = GMAP, col_names = FALSE, show_col_types = FALSE
  ) %>%
    rename(chr = X1, snp = X2, cM = X3, bp = X4)
} else {
  calls <- read.rapid.bed(IBD.FILE, MIN.BP.SIZE, MIN.CM.SIZE, GMAP)
  gmap <- readr::read_table(
    file = GMAP, col_names = FALSE, show_col_types = FALSE
  ) %>%
    rename(idx = X1, cM = X2)
}

# print ibd segments in BED format ----
readr::write_tsv(calls, FF.OUT, col_names = TRUE)

# filter out hot bins (bins with too great of a number of segemnts) ----
chr.len.cM <- max(as.numeric(gmap$cM), na.rm = TRUE)
hot.all <- filter.ibd.by.hot.bins(
  calls = calls, chr.len.cM = chr.len.cM, bin.size.cM = HOT.BIN.SIZE, 
  sd.thresh = HOT.BIN.SD, trim.frac = 0.03
  )

# print ibd segments filter by hot bins in BED format ----
readr::write_tsv(hot.all$calls.filtered, FF.HOT.OUT, col_names = TRUE)
