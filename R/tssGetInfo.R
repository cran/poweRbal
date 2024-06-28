#' Get information on included tree shape statistics
#'
#' \code{getTSSnames} - Returns the full names (character/expression) of the
#' TSS.
#'
#' @author Sophie Kersting
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSSnames} Vector of characters/expressions.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSnames(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSnames <- function(tss_shorts){
  namesTSS <- rep(NA, length(tss_shorts))
  for(i in 1:length(tss_shorts)){
    namesTSS[i] <- tssInfo[[tss_shorts[i]]]$name
  }
  return(namesTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSsimple} - Returns the simple names (character/expression) of the
#' TSS.
#'
#' @return \code{getTSSsimple} Vector of characters/expressions.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSsimple(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSsimple <- function(tss_shorts){
  namesTSS <- rep(NA, length(tss_shorts))
  for(i in 1:length(tss_shorts)){
    namesTSS[i] <- tssInfo[[tss_shorts[i]]]$simple
  }
  return(namesTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSScolors} - Returns the colors of the TSS.
#'
#' @return \code{getTSScolors} Vector of characters (color names).
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSScolors(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSScolors <- function(tss_shorts){
  colorsTSS <- rep(NA, length(tss_shorts))
  for(i in 1:length(tss_shorts)){
    colorsTSS[i] <- tssInfo[[tss_shorts[i]]]$col
  }
  return(colorsTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSsafe_n} - Returns the ranges of n that can be safely used.
#'
#' @return \code{getTSSsafe_n} Numeric matrix, one row per TSS and two columns
#' with lower and upper limit.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSsafe_n(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSsafe_n <- function(tss_shorts){
  safe_nTSS <- matrix(NA, nrow = length(tss_shorts), ncol = 2,
                     dimnames = list(tss_shorts, NULL))
  for(i in 1:length(tss_shorts)){
    safe_nTSS[i,] <- tssInfo[[tss_shorts[i]]]$safe_n
  }
  return(safe_nTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSStype} - Returns the ranges of n that can be safely used.
#'
#' @return \code{getTSStype} Vector of characters (types as factors).
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSStype(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSStype <- function(tss_shorts){
  typeTSS <- rep(NA, length(tss_shorts))
  for(i in 1:length(tss_shorts)){
    typeTSS[i] <- tssInfo[[tss_shorts[i]]]$type
  }
  typeTSS <- factor(typeTSS, levels = c("tss", "bali", "imbali"))
  return(typeTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSonly_bin} - Returns TRUE/FALSE vector: TRUE if TSS is only for
#' binary trees and FALSE otherwise.
#'
#' @return \code{getTSSonly_bin} Logical vector.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSonly_bin(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSonly_bin <- function(tss_shorts){
  only_binTSS <- rep(NA, length(tss_shorts))
  for(i in 1:length(tss_shorts)){
    only_binTSS[i] <- tssInfo[[tss_shorts[i]]]$only_binary
  }
  return(only_binTSS)
}
#' Get information on included tree shape statistics
#'
#' \code{getAllTSS} - Returns the short names of all TSS that are safe to
#' use for the specified \code{n}, have one of the specified \code{types} and
#' can be applied to (non-)binary trees (\code{not_only_bin}).
#'
#' @param n Integer value or vector of integer values, that
#' specifies the number(s) of leaves.
#' If NULL (default), then \code{getAllSafeTSS} returns the short names of all
#' TSS contained in \code{tssInfo}.
#' @param not_only_bin Select TRUE if you also want to analyze non-binary trees
#' and therefore want to filter out any TSS that only work on binary trees.
#' Otherwise, select FALSE (default) if all TSS are applicable.
#' @param types Character vector, that specifies all permissible TSS types. The
#' vector may contain a subset of \code{c("tss", "bali", "imbali")} to indicate
#' if balance indices, imbalance indices or mere TSS should be included. By
#' default all types are permissible.
#'
#' @return \code{getAllTSS} Character vector of short names of TSS
#' contained in \code{tssInfo}.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getAllTSS(n = c(3,30))
getAllTSS <- function(n = NULL, not_only_bin = FALSE,
                      types = c("tss", "bali", "imbali")){
  select_TSS <- NULL
  for(i in 1:length(tssInfo)){
    is_permissible <- TRUE
    # Check range of n.
    if(!is.null(n) &&
       (sum(n<tssInfo[[i]]$safe_n[1] |
           n>tssInfo[[i]]$safe_n[2] )>0)){
      is_permissible <- FALSE
    }
    # Check if applicable to (non-)binary trees
    if(not_only_bin &&
       tssInfo[[i]]$only_binary){
      is_permissible <- FALSE
    }
    # Check if type OK.
    if(!tssInfo[[i]]$type %in% types){
      is_permissible <- FALSE
    }
    # Add if permissible.
    if(is_permissible){
      select_TSS <- c(select_TSS, tssInfo[[i]]$short)
    }
  }
  return(select_TSS)
}
