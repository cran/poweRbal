#' Print and summary method for poweRbal_data objects
#'
#' This function prints the contents of an object of class \code{poweRbal_data}.
#' It provides a brief summary of the object structure and its contents.
#'
#' @param x An object of class \code{poweRbal_data}, which is a list containing
#' one mandatory element, \code{power}, and several optional elements:\cr
#' - \code{power}: A numeric matrix containing the power values (one row per
#' TSS and one column per alternative model).\cr
#' - \code{accept_regions}: A numeric matrix containing information on the
#' region of acceptance (one row per TSS and four columns).\cr
#' - \code{CIradius}: A numeric matrix containing the confidence interval
#' radii (one row per TSS and one column per alternative model).\cr
#' - \code{actual_sample_sizes}: A numeric vector containing the actual sample
#' sizes under each alternative model, as some models do not always
#' successfully generate trees.\cr
#' - \code{alt_model_params}: A numeric vector (one element per alternative
#' model) containing the values of a tree model parameter. This is only
#' suitable if the alternative models all belong to the same tree model
#' family and differ only in one parameter.\cr
#' - Other input data from the \code{powerComp()} function, such as \code{tss},
#' \code{null_model}, \code{alt_models}, \code{n}, \code{distribs},
#' \code{N_null}, \code{N_alt}, \code{test_type}, \code{correction}, and
#' \code{sig_lvl}.
#' @param ... Additional arguments passed to the \code{print} or \code{summary}
#' function.
#'
#' @return \code{print.poweRbal_data} No return value, as the primary purpose
#' of this function is the side effect (printing).
#'
#' @rdname poweRbal_data-printsummary
#' @export
#'
#' @examples
#' # Printing a 'poweRbal_data' object:
#' pc1 <- powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                  alt_models = list(list("aldous",-1), "pda", "etm"),
#'                  n = 8L, N_null = 40L, N_alt = 20L)
#' pc1
print.poweRbal_data <- function(x, ...) {
  if(is.null(x$power)){
    cat(paste("Incomplete object of class 'poweRbal_data' (element 'power'",
              "missing).\n"))
  } else {
    cat("Object of class 'poweRbal_data' with the following components:\n")
    print(names(x), ...)
    cat("TSS:\n")
    if(is.null(x$tss)){
      print(rownames(x$power), ...)
    } else {
      print(x$tss, ...)
    }
    cat("n:\n")
    if(is.null(x$n)){
      print(NULL, ...)
    } else {
      print(x$n, ...)
    }
    cat("Null model:\n")
    if(is.null(x$null_model)){
      print(NULL, ...)
    } else {
      print(toString(x$null_model), ...)
    }
    cat("Alternative models:\n")
    if(is.null(x$null_model)){
      print(colnames(x$power), ...)
    } else {
      print(sapply(1:length(x$alt_models), function(X) {
        paste(unlist(x$alt_models[[X]]), collapse=", ")}), ...)
    }
  }
}

#' Print and summary method for poweRbal_data objects
#'
#' This function provides a summary of an object of class \code{poweRbal_data}.
#' It offers a high-level overview of the contents and their structure.
#'
#' @param object An object of class \code{poweRbal_data} (see \code{x} for
#' more details).
#'
#' @return \code{summary.poweRbal_data} No return value, as the primary purpose
#' of this function is the side effect (printing summary).
#'
#' @rdname poweRbal_data-printsummary
#' @export
#'
#' @examples
#' # Summary of a 'poweRbal_data' object:
#' summary(pc1)
summary.poweRbal_data <- function(object, ...) {
  if(is.null(object$power)){
    cat(paste("Incomplete object of class 'poweRbal_data' (element 'power'",
              "missing).\n"))
  } else {
    cat("Summary of 'poweRbal_data' object with the following components:\n")
    print(toString(names(object)), ...)

    cat("\nElement types:\n")
    element_types <- sapply(object, class)
    print(toString(element_types),...)

    cat("\nTSS:\n")
    if(is.null(object$tss)){
      print(rownames(object$power), ...)
    } else {
      print(object$tss, ...)
    }
    cat("n:\n")
    if(is.null(object$n)){
      print(NULL, ...)
    } else {
      print(object$n, ...)
    }
    cat("Null model:\n")
    if(is.null(object$null_model)){
      print(NULL, ...)
    } else {
      print(paste(unlist(object$null_model), collapse=", "), ...)
    }
    cat("Alternative models:\n")
    if(is.null(object$alt_models)){
      print(colnames(object$power), ...)
    } else {
      print(sapply(1:length(object$alt_models), function(X) {
        paste(unlist(object$alt_models[[X]]), collapse=", ")}), ...)
    }

    cat("\nTest method: ")
    if(is.null(object$test_type)){
      cat("\n")
    } else {
      cat(object$test_type,"\n")
    }
    cat("Level of significance: ")
    if(is.null(object$sig_lvl)){
      cat("\n")
    } else {
      cat(object$sig_lvl,"\n")
    }
    cat("Correction method: ")
    if(is.null(object$correction)){
      cat("\n")
    } else {
      cat(object$correction,"\n")
    }
    cat("Distribution: ")
    if(is.null(object$distribs)){
      cat("\n")
    } else {
      if(object$distribs == "exact_if_possible"){
        if(is.null(object$n) && is.null(object$N_null) && is.null(object$null_model)){
          cat(object$distribs,"\n")
        } else {
          if(object$n <= 20 && (object$null_model[[1]] == "yule" ||
                           object$null_model[[1]] == "pda" ||
                           object$null_model[[1]] == "etm")){
            cat(object$distribs, " (in this case exact)\n", sep = "")
          } else {
            cat(object$distribs, " (in this case sampled with ",
                "'N_null'=", object$N_null, ")\n", sep = "")
          }
        }
      } else if(object$distribs == "sampled"){
        if(is.null(object$N_null)){
          cat(object$distribs,"\n")
        } else {
          cat(object$distribs," with 'N_null'=", object$N_null, "\n", sep = "")
        }
      } else {
        cat(object$distribs,"\n")
      }
    }

    cat("\n\nFirst few elements of main component(s):\n")
    for (name in (names(object)[names(object) %in% c("power","accept_regions",
                                           "CIradius",
                                           "actual_sample_sizes")])) {

      if(name == "actual_sample_sizes"){
        if(is.null(object$N_alt)){
          cat("Component '", name, "':\n", sep = "")
          print(utils::head(object[[name]]), ...)
          cat("\n")
        } else {
          cat("Component '", name, "' ('N_alt'=",object$N_alt,"):\n", sep = "")
          print(utils::head(object[[name]]), ...)
          cat("\n")
        }
      } else {
        cat("Component '", name, "':\n", sep = "")
        print(utils::head(object[[name]]), ...)
        cat("\n")
      }
    }


  }
}
#' Plot method for poweRbal_data objects
#'
#' This function generates a plot for an object of class \code{poweRbal_data}.
#' Creates a bar plot if \code{alt_model_params} and \code{x$alt_model_params}
#' = NULL and a line plot otherwise if this information is given.
#'
#' @param x An object of class \code{poweRbal_data}, which is a list containing
#' one mandatory element, \code{power}, and several optional elements:\cr
#' - \code{power}: A numeric matrix containing the power values (one row per
#' TSS and one column per alternative model).\cr
#' - \code{accept_regions}: A numeric matrix containing information on the
#' region of acceptance (one row per TSS and four columns).\cr
#' - \code{CIradius}: A numeric matrix containing the confidence interval
#' radii (one row per TSS and one column per alternative model).\cr
#' - \code{actual_sample_sizes}: A numeric vector containing the actual sample
#' sizes under each alternative model, as some models do not always
#' successfully generate trees.\cr
#' - \code{alt_model_params}: A numeric vector (one element per alternative
#' model) containing the values of a tree model parameter. This is only
#' suitable if the alternative models all belong to the same tree model
#' family and differ only in one parameter.\cr
#' - Other input data from the \code{powerComp()} function, such as \code{tss},
#' \code{null_model}, \code{alt_models}, \code{n}, \code{distribs},
#' \code{N_null}, \code{N_alt}, \code{test_type}, \code{correction}, and
#' \code{sig_lvl}.
#' @param tss_names Vector of characters/expression of the TSS names (default
#' = NULL). If none are provided, \code{x$tss} is used for the names
#' of the TSS if existent and otherwise the row names of \code{x$power}
#' are used.
#' @param tss_colors Vector of colors for the TSS (default = NULL).
#' @param sig_lvl Level of significance (default=0.05, must be >0 and <1)
#' depicted as a dashed horizontal line. Not depicted if set to NULL.
#' @param legend_pos Character specifying where the legend is displayed
#' (default = "topright"). No legend is displayed if set to NULL.
#'
#' @param alt_model_names Vector of characters/expression of the model names
#' (default = NULL). If none are provided, the column names of
#' \code{x$power} are used as names for the models.\cr
#' Only used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.
#'
#' @param alt_model_params Numeric vector containing the parameter
#' values of the representatives of the tree model (default = NULL). If none
#' are provided, \code{x$alt_model_params} is used if existent.
#' @param tss_ltys Vector of line types for the TSS (default = NULL).\cr
#' Not used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.
#' @param alt_model_family Vector of characters/expressions of the name of the
#' tree model family and of the parameter (default = NULL), e.g.
#' \code{c("Aldous\'", expression(beta))}.
#' If none is provided, the first column name of \code{x$power} is used.\cr
#' Not used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.
#' @param ... Additional arguments passed to the \code{plot} function.
#'
#' @return \code{plot.poweRbal_data} No return value, as the primary purpose
#' of this function is the side effect (plotting).
#'
#' @rdname poweRbal_data-plot
#' @export
#'
#' @examples
#' # Plotting a 'poweRbal_data' object:
#' pc1 <- powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                  alt_models = list(list("aldous",-1), "pda", "etm"),
#'                  n = 8L, N_null = 40L, N_alt = 20L)
#' plot(pc1)
#' # Plotting a power comparison with a tree model family
#' pc2 <- powerComp(tss = c("Sackin", "Colless", "B1I"),
#'                alt_models = list(list("aldous", -1.5),
#'                                  list("aldous", -1),list("aldous", -0.5),
#'                                  list("aldous", 0),list("aldous", 0.5)),
#'                n=20L, N_null = 20L, N_alt = 10L, distribs = "sampled")
#' # Create a bar plot or ...
#' plot(pc2)
#' # ... a line plot by specifying 'alt_model_params'.
#' plot(pc2, alt_model_params = c(-1.5,-1,-0.5,0,0.5),
#'      tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
#'      tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
#'      alt_model_family = c("Aldous\'", expression(beta)),
#'      ylim = c(0,1))
plot.poweRbal_data <- function(x, tss_names = NULL, tss_colors = NULL,
                               sig_lvl = NULL, legend_pos = "topright",
                               alt_model_names = NULL,
                               alt_model_params = NULL, tss_ltys = NULL,
                               alt_model_family = NULL, ...) {
  if(is.null(x$power)){
    stop("Power matrix not found.")
  }
  if(is.null(tss_names)){
    if(!is.null(x$tss)){
      tss_names <- x$tss
    } else {
      tss_names <- rownames(x$power)
    }
  }
  if(is.null(tss_colors)){
    if(nrow(x$power) <= length(.my_cols)){
      tss_colors <- .my_cols[1:nrow(x$power)]
    } else {
      tss_colors <- c(.my_cols, grDevices::rainbow(nrow(x$power)-
                                                     length(.my_cols)))
    }
  }
  if(is.null(alt_model_params) && is.null(x$alt_model_params)){
    # Create bar plot: ---------------------------------------------------------
    if(is.null(alt_model_names)){
      if(!is.null(x$alt_models)){
        alt_model_names <- sapply(1:length(x$alt_models), function(X) {
          paste(unlist(x$alt_models[[X]]), collapse=" ")})
      } else {
        alt_model_names <- colnames(x$power)
      }
    }
    barCenters <- graphics::barplot(height = x$power, beside = TRUE,
                                    col = tss_colors,
                                    names.arg = alt_model_names,...)
    barCenters
    if(!is.null(x$CIradius)){
      graphics::segments(barCenters, x$power - x$CIradius,
                         barCenters, x$power + x$CIradius, ...)
    }
    if(!is.null(sig_lvl)){
      graphics::abline(h = sig_lvl, col="black", lty=2)
    }
    if(!is.null(legend_pos)){
      graphics::legend(legend_pos, legend = tss_names, col = tss_colors,
                       pch = 16, bty = "n", pt.cex = 1.5, cex = 1,
                       text.col = "black", horiz = F, inset = c(0.05, 0))
    }
  } else {
    # Create line plot: --------------------------------------------------------
    if(is.null(alt_model_params)){
      alt_model_params <- x$alt_model_params
    }
    if(is.null(tss_ltys)){
      tss_ltys <- rep(1, nrow(x$power))
    }
    if(is.null(alt_model_family) && !is.null(x$alt_models)){
      alt_model_family <- c(x$alt_models[[1]][1], "parameter")
    }
    # Basic plot.
    if(is.null(alt_model_family)){
      graphics::plot(0, type="n", xlim=c(min(alt_model_params),
                                         max(alt_model_params)),
                     xaxt = "n", ...)
      graphics::axis(1, at = alt_model_params, labels=alt_model_params)
    } else {
      graphics::plot(0, type="n", main = alt_model_family[1],
                     xlab = alt_model_family[2],
                     xlim=c(min(alt_model_params),
                            max(alt_model_params)), xaxt = "n", ...)
      graphics::axis(1, at = alt_model_params, labels=alt_model_params)
    }
    if(!is.null(sig_lvl)){
      graphics::abline(h = sig_lvl, col="black", lty=2)
    }
    # Add points and lines to the plot for each TSS.
    for(i in 1:nrow(x$power)){
      if(!is.null(x$CIradius)){
        graphics::polygon(c(rev(alt_model_params), alt_model_params),
                          c(rev(x$power[i,]+x$CIradius[i,]),
                                x$power[i,]-x$CIradius[i,]),
                          col = scales::alpha(tss_colors[i], 0.1), border = NA)
      }
      graphics::lines(alt_model_params, x$power[i,], col = tss_colors[i],
                      lty=tss_ltys[i], lwd=2)
      graphics::points(alt_model_params, x$power[i,],
                       col = scales::alpha(tss_colors[i], 0.6),
                       pch=16,  cex=2)

    }
    if(!is.null(legend_pos)){
      graphics::legend(legend_pos, legend = tss_names, col = tss_colors,
                       pch = 16, bty = "n", pt.cex = 0.5, cex = 0.5,
                       text.col = "black", horiz = F, inset = c(0.05, 0))
    }
  }
}
.my_cols <- c("royalblue1", "lightpink3", "lightblue4", "darkblue",
              "deepskyblue4", "darkolivegreen", "mediumorchid", "brown",
              "darkolivegreen2", "tan3", "firebrick1", "darkgoldenrod1",
              "maroon4", "lemonchiffon3", "darkorange4", "chartreuse3",
              "mediumturquoise", "gray23", "deeppink3", "lightgoldenrod" )
