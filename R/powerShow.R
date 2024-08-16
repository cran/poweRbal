#' Function for displaying TSS distributions
#'
#' \code{showTSSdata} - This function plots histograms of TSS data.
#'
#' @param tss_data Numeric matrix of TSS values (one row per TSS). The row names
#' are used as names for the TSS.
#' @param main Title (default = NULL). A generic title is created by default.
#' @param xlab Label of x-axis (default = NULL). A generic label is created
#' by default.
#' @param sig_lvl Level of significance (default=0.05, must be >0 and <1).
#' @param ... Add further specifications for \code{plot()}.
#'
#' @return \code{showTSSdata} No return value, called for side effects
#' (plotting).
#'
#' @export
#' @rdname powerShow
#'
#' @examples
#' showTSSdata(tss_data = getTSSdata_trees(tss = c("Colless", "Sackin"),
#'                 treeList = lapply(1:20L, function(x) genYuleTree(10))),
#'                 breaks=15)
showTSSdata <- function(tss_data, main = NULL, xlab = NULL,
                        sig_lvl = 0.05, ...) {
  if(is.null(dim(tss_data))){
    # Convert a vector of values to a matrix with one row.
    tss_data <- matrix(tss_data, nrow =1, byrow = TRUE,
                       dimnames = list("TSS", NULL))
  }
  if(is.null(main)){
    # Convert a vector of values to a matrix with one row.
    main <- paste0("Histogram of ", dimnames(tss_data)[[1]][1],
    " values (N=",ncol(tss_data),")")
  }
  if(is.null(xlab)){
    # Convert a vector of values to a matrix with one row.
    xlab <- paste(dimnames(tss_data)[[1]][1],"values")
  }
  # Plot the histogram of first TSS and add quantiles.
  graphics::hist(tss_data[1,], main = main, xlab = xlab, ...)
  graphics::abline(v=stats::quantile(tss_data[1,],
                                     probs = c(0+sig_lvl/2,1-sig_lvl/2),
                    na.rm = TRUE, names = TRUE, type = 1),
         col=c("black","black"),lwd=c(3,3), lty=c(2,2))
  # Plot next.
  if(nrow(tss_data)>1){
    showTSSdata(matrix(tss_data[-1,], nrow = nrow(tss_data)-1, byrow = TRUE,
                       dimnames = list(dimnames(tss_data)[[1]][-1], NULL)),
                main = NULL, xlab = NULL,
                sig_lvl, ...)
  }
}
