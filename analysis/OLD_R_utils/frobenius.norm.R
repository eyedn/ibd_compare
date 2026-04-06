###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           frobenius.norm.R
###############################################################################

frobenius.norm <- function(A, B) {
  stopifnot(
    is.matrix(A), is.matrix(B),
    all(dim(A) == dim(B)),
    identical(rownames(A), rownames(B)),
    identical(colnames(A), colnames(B))
  )
  sqrt(sum((A - B)^2, na.rm = TRUE))
}
