#' Functions for displaying TSS distributions and TSS power
#'
#' \code{showTSSdata} - This function plots histograms of TSS data.
#'
#' @author Sophie Kersting
#'
#' @param tss_data Numeric matrix of TSS values (one row per TSS). The row names
#' are used as names for the TSS.
#' @param main Title (default = NULL). A generic title is created by default.
#' @param xlab Label of x-axis (default = NULL). A generic label is created
#' by default.
#' @param test_type Determines the method. Available are: \cr
#' - "two-tailed" (default)
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
                        test_type = "two-tailed", sig_lvl = 0.05, ...) {
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
                    na.rm = TRUE, names = TRUE, type = 7),
         col=c("black","black"),lwd=c(3,3), lty=c(2,2))
  # Plot next.
  if(nrow(tss_data)>1){
    showTSSdata(matrix(tss_data[-1,], nrow = nrow(tss_data)-1, byrow = TRUE,
                       dimnames = list(dimnames(tss_data)[[1]][-1], NULL)),
                main = NULL, xlab = NULL, test_type,
                sig_lvl, ...)
  }
}
#' Functions for displaying TSS distributions and TSS power
#'
#' \code{showPower} - This function displays the power of TSS under
#' different alternative models.
#'
#' @param power_data Numeric matrix, one row per TSS and one column per
#' alternative model. Or alternatively a list containing two such matrices,
#' one for the power values and one for the confidence interval radii. In the
#' latter case, error bars are also depicted.
#' @param tss_names Vector of characters/expression of the TSS names (default
#' = NULL). If none are provided, the row names of \code{power_data} are used
#' as names for the TSS.
#' @param tss_colors Vector of colors for the TSS (default = NULL).
#' @param model_names Vector of characters/expression of the model names
#' (default = NULL). If none are provided, the column names of \code{power_data}
#' are used as names for the models.
#' @param sig_lvl Level of significance (default=0.05, must be >0 and <1). NO
#' horizontal line is depicted if set to NULL.
#'
#' @return \code{showPower} No return value, called for side effects
#' (plotting).
#'
#' @export
#' @rdname powerShow
#'
#' @examples
#' # With error bars:
#' showPower(power_data = powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                                      alt_models = list("pda", "etm"), n=10L,
#'                                      sample_size = 50L,
#'                                      with_CIradius = TRUE),
#'               tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
#'               tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
#'               model_names = c("PDA", "ETM"),
#'               main = "Power (Yule as null model, n = 10, N=50)",
#'               ylim = c(0,1), ylab = "Power (null model rejected)")
#' # Without error bars:
#' showPower(power_data = powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                                      alt_models = list("pda", "etm"), n=10L,
#'                                      sample_size = 50L),
#'               tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
#'               tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
#'               model_names = c("PDA", "ETM"),
#'               main = "Power (Yule as null model, n = 10, N=50)",
#'               ylim = c(0,1), ylab = "Power (null model rejected)")
showPower <- function(power_data, tss_names = NULL, tss_colors = NULL,
                      model_names = NULL, sig_lvl = 0.05, ...){
  radius_data <- NULL
  if(is.list(power_data)){
    radius_data <- power_data[[2]]
    power_data <- power_data[[1]]
  }
  if(is.null(tss_names)){
    tss_names <- rownames(power_data)
  }
  if(is.null(model_names)){
    model_names <- colnames(power_data)
  }
  if(is.null(tss_colors)){
    if(nrow(power_data) <= length(.my_cols)){
      tss_colors <- .my_cols[1:nrow(power_data)]
    } else {
      tss_colors <- c(.my_cols, grDevices::rainbow(nrow(power_data)-
                                                    length(.my_cols)))
    }

  }
  barCenters <- graphics::barplot(height = power_data, beside = TRUE,
                    col = tss_colors,
                    names.arg = model_names,...)
  barCenters
  if(!is.null(radius_data)){
    graphics::segments(barCenters, power_data - radius_data,
                       barCenters, power_data + radius_data, lwd = 1.5)
  }
  if(!is.null(sig_lvl)){
    graphics::abline(h = sig_lvl, col="black", lty=2)
  }
  graphics::legend("topright", legend = tss_names, col = tss_colors, pch = 16,
         bty = "n", pt.cex = 1.5, cex = 1, text.col = "black", horiz = F,
         inset = c(0.05, 0))
}
#' Functions for displaying TSS distributions and TSS power
#'
#' \code{showPower_param} - This function displays the power of TSS under
#' different representatives of a family of tree models which vary only in one
#' parameter.
#'
#' @author Sophie Kersting
#'
#' @param power_data Numeric matrix, one row per TSS and one column per
#' alternative model.
#' @param tss_ltys Vector of line types for the TSS (default = NULL).
#' @param model_name Vector of characters/expressions of the name of the tree
#' model family and of the parameter (default = NULL), e.g.
#' \code{c("Aldous\'", expression(beta))}.
#' If none are provided, the first column name of \code{power_data} is used.
#' @param model_param Numeric vector containing the parameter values of the
#' representatives of the tree model.
#'
#' @return \code{showPower_param} No return value, called for side effects
#' (plotting).
#'
#' @export
#' @rdname powerShow
#'
#' @examples
#' # With confidence bands:
#' showPower_param(power_data = powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                                    alt_models = list(list("aldous", 0.5),
#'                                    list("aldous", 0),list("aldous", -0.5),
#'                                    list("aldous", -1),list("aldous", -1.5)),
#'                                    n=20L, sample_size = 50L,
#'                                    with_CIradius = TRUE),
#'               tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
#'               tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
#'               model_name = c("Aldous-beta splitting model", "beta"),
#'               model_param = c(0.5,0,-0.5,-1,-1.5),
#'               ylim = c(0,1), ylab = "Power (null model rejected)")
#' # Without confidence bands:
#' showPower_param(power_data = powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                                    alt_models = list(list("aldous", 0.5),
#'                                    list("aldous", 0),list("aldous", -0.5),
#'                                    list("aldous", -1),list("aldous", -1.5)),
#'                                    n=20L, sample_size = 50L),
#'               tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
#'               tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
#'               model_name = c("Aldous-beta splitting model", "beta"),
#'               model_param = c(0.5,0,-0.5,-1,-1.5),
#'               ylim = c(0,1), ylab = "Power (null model rejected)")
showPower_param <- function(power_data, tss_names = NULL, tss_colors = NULL,
                            tss_ltys = NULL, model_name = NULL,
                            model_param, sig_lvl = 0.05, ...){
  radius_data <- NULL
  if(is.list(power_data)){
    radius_data <- power_data[[2]]
    power_data <- power_data[[1]]
  }
  if(is.null(tss_names)){
    tss_names <- rownames(power_data)
  }
  if(is.null(model_name)){
    model_name <- c(colnames(power_data)[1],"parameter")
  }
  if(is.null(tss_colors)){
    if(nrow(power_data)<= length(.my_cols)){
      tss_colors <- .my_cols[1:nrow(power_data)]
    } else {
      tss_colors <- c(.my_cols, grDevices::rainbow(nrow(power_data)-
                                                    length(.my_cols)))
    }
  }
  if(is.null(tss_ltys)){
    tss_ltys <- 1:nrow(power_data)
  }
  # Basic plot.
  graphics::plot(0, type="n", main = model_name[1], xlab = model_name[2],
       xlim=c(min(model_param), max(model_param)), xaxt = "n",
       ...)
  graphics::axis(1, at = model_param, labels=model_param)
  if(!is.null(sig_lvl)){
    graphics::abline(h = sig_lvl, col="black", lty=2)
  }
  # Add points and lines to the plot for each TSS
  for(i in 1:nrow(power_data)){
    if(!is.null(radius_data)){
      graphics::polygon(c(rev(model_param), model_param),
                        c(rev(power_data[i,]+radius_data[i,]),
                          power_data[i,]-radius_data[i,]),
                          col = scales::alpha(tss_colors[i], 0.1), border = NA)
    }
    graphics::lines(model_param, power_data[i,], col = tss_colors[i],
          lty=tss_ltys[i],  lwd=2)
    graphics::points(model_param, power_data[i,],
                     col = scales::alpha(tss_colors[i], 0.6),
           pch=16,  cex=2)

  }
  graphics::legend("topleft", legend = tss_names, col = tss_colors, pch = 16,
         bty = "n", pt.cex = 0.5, cex = 0.5, text.col = "black", horiz = F,
         inset = c(0.05, 0))
}
.my_cols <- c("royalblue1", "lightpink3", "lightblue4", "darkblue",
              "deepskyblue4", "darkolivegreen", "mediumorchid", "brown",
              "darkolivegreen2", "tan3", "firebrick1", "darkgoldenrod1",
              "maroon4", "lemonchiffon3", "darkorange4", "chartreuse3",
              "mediumturquoise", "gray23", "deeppink3", "lightgoldenrod" )
