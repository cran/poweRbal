# Internal helper: extract a single named field from tssInfo for a vector of shorts.
.getTSSfield <- function(tss_shorts, field) {
  sapply(tss_shorts, function(x) tssInfo[[x]][[field]], USE.NAMES = FALSE)
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSnames} - Returns the full names (character/expression) of the
#' TSS.
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
getTSSnames <- function(tss_shorts) {
  .getTSSfield(tss_shorts, "name")
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSsimple} - Returns the simple names (character/expression) of the
#' TSS.
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSSsimple} Vector of characters/expressions.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSsimple(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSsimple <- function(tss_shorts) {
  .getTSSfield(tss_shorts, "simple")
}
#' Get information on included tree shape statistics
#'
#' \code{getTSScolors} - Returns the colors of the TSS.
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSScolors} Vector of characters (color names).
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSScolors(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSScolors <- function(tss_shorts) {
  .getTSSfield(tss_shorts, "col")
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSsafe_n} - Returns the ranges of n that can be safely used.
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSSsafe_n} Numeric matrix, one row per TSS and two columns
#' with lower and upper limit.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSsafe_n(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSsafe_n <- function(tss_shorts) {
  matrix(t(sapply(tss_shorts, function(x) tssInfo[[x]]$safe_n)),
    ncol = 2, dimnames = list(tss_shorts, NULL)
  )
}
#' Get information on included tree shape statistics
#'
#' \code{getTSStype} - Returns the types of the TSS, i.e., whether they are
#' balance or imbalance indices, or simple tree shape statistics.
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSStype} Factor vector with levels
#' \code{c("tss", "bali", "imbali")}.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSStype(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSStype <- function(tss_shorts) {
  factor(.getTSSfield(tss_shorts, "type"), levels = c("tss", "bali", "imbali"))
}
#' Get information on included tree shape statistics
#'
#' \code{getTSSonly_bin} - Returns TRUE/FALSE vector: TRUE if TSS is only for
#' binary trees and FALSE otherwise.
#'
#' @param tss_shorts Vector of short names (characters) of TSS contained in
#' \code{tssInfo}.
#'
#' @return \code{getTSSonly_bin} Logical vector.
#'
#' @export
#' @rdname tssGetInfo
#'
#' @examples
#' getTSSonly_bin(tss_shorts = c("Sackin", "Colless", "B1I"))
getTSSonly_bin <- function(tss_shorts) {
  .getTSSfield(tss_shorts, "only_binary")
}
#' Get information on included tree shape statistics
#'
#' \code{getAllTSS} - Returns the short names of all TSS that are safe to
#' use for the specified \code{n}, have one of the specified \code{types} and
#' can be applied to (non-)binary trees (\code{not_only_bin}).
#'
#' @param n Integer value or vector of integer values, that
#' specifies the number(s) of leaves.
#' If a vector is provided, only TSS that are safe for \emph{all} values in the
#' vector are returned (i.e., the intersection of safe TSS across all specified
#' leaf counts).
#' If NULL (default), then \code{getAllTSS} returns the short names of all
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
#' getAllTSS(n = c(3, 30))
getAllTSS <- function(n = NULL, not_only_bin = FALSE,
                      types = c("tss", "bali", "imbali")) {
  select_TSS <- NULL
  for (i in seq_along(tssInfo)) {
    is_permissible <- TRUE
    # Check range of n.
    if (!is.null(n) &&
      (sum(n < tssInfo[[i]]$safe_n[1] |
        n > tssInfo[[i]]$safe_n[2]) > 0)) {
      is_permissible <- FALSE
    }
    # Check if applicable to (non-)binary trees
    if (not_only_bin &&
      tssInfo[[i]]$only_binary) {
      is_permissible <- FALSE
    }
    # Check if type OK.
    if (!tssInfo[[i]]$type %in% types) {
      is_permissible <- FALSE
    }
    # Add if permissible.
    if (is_permissible) {
      select_TSS <- c(select_TSS, tssInfo[[i]]$short)
    }
  }
  return(select_TSS)
}
