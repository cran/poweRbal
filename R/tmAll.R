# Internal helper: generate Ntrees trees using gen_func() and return the
# appropriate class. NULL results (from stochastic models that may fail) are
# silently dropped and a message reports the shortfall.
.wrap_trees <- function(Ntrees, gen_func) {
  if (Ntrees == 1L) {
    return(gen_func())
  }
  phylos <- lapply(seq_len(Ntrees), function(x) gen_func())
  nulls <- which(sapply(phylos, is.null))
  if (length(nulls) > 0) {
    phylos <- phylos[-nulls]
    message(paste0(
      "Only ", length(phylos), " trees generated instead of ",
      Ntrees, ".\n"
    ))
  }
  structure(phylos, class = "multiPhylo")
}
#' Generation of rooted binary trees under a given tree model
#'
#' \code{genTrees} - Is a wrapper function that generates
#' \code{Ntrees}-many rooted binary trees with the given number of \code{n}
#' leaves under any tree model \code{tm} contained in this package (more
#' details on the available models are given in the parameter information
#' for \code{tm}).
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0. \cr
#' Due to the restrictions of the \code{phylo} or \code{multiPhylo} format,
#' the number of leaves must be at least 2 since there
#' must be at least one edge.
#' @param Ntrees Integer value (default = 1) that specifies the desired number
#' of generated trees.
#' @param tm Character or list specifying the tree model under which the trees
#' should be generated as well as their parameters. Available are: \cr
#' - \code{"yule"} - Yule model. \cr
#' - \code{"pda"} - PDA model. \cr
#' - \code{"etm"} - ETM. \cr
#' - \code{list("aldous", BETA)} - Aldous' beta splitting model with parameter
#' \code{BETA} >= -2. \cr
#' - \code{list("ford", ALPHA)} - Ford's alpha model with parameter
#' \code{ALPHA} >= 0 and <= 1. \cr
#' - \code{list("alt-birth-death", BIRTHRATE, DEATHRATE)} or \cr
#' \code{list("alt-birth-death", BIRTHRATE, DEATHRATE, TRIES)} -
#' Alternative birth-death model with parameters \code{BIRTHRATE} >0 and
#' \code{DEATHRATE}>= 0. \cr
#' - \code{list("density", BIRTHRATE, EQUILIB)} or \cr
#' \code{list("density", BIRTHRATE, EQUILIB, TRIES, TIMEperTRY)} -
#' Density dependent model with parameters \code{BIRTHRATE} >0 and
#' \code{EQULIB}>= 1. \cr
#' - \code{list("BiSSE", BIRTHRATES, DEATHRATES, TRANSRATES)} or \cr
#' \code{list("BiSSE", BIRTHRATES, DEATHRATES, TRANSRATES, TRIES, TIMEperTRY)} -
#' BiSSE model with parameters \code{BIRTHRATES} (vector with 2 values >=0, one
#' value >0),
#' \code{DEATHRATES} (vector with 2 values >=0), and
#' \code{TRANSRATES} (vector with 2 values >=0, one
#' value >0). \cr
#' - \code{list("DCO_sym", ZETA)} or \cr
#' \code{list("DCO_sym", ZETA, STARTING_RATE)} - Symmetric direct-children-only
#' with parameter \code{ZETA} > 0 and optionally \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("DCO_asym", ZETA)} or \cr
#' \code{list("DCO_asym", ZETA, STARTING_RATE)} - Asymmetric
#' direct-children-only with parameter \code{ZETA} > 0 and optionally
#' \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("IF_sym", ZETA)} or \cr
#' \code{list("IF_sym", ZETA, STARTING_RATE)} - Symmetric inherited fertility
#' with parameter \code{ZETA} > 0 and optionally \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("IF_asym", ZETA)} or \cr
#' \code{list("IF_asym", ZETA, STARTING_RATE)} - Asymmetric inherited fertility
#' with parameter \code{ZETA} > 0 and optionally \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("IF-diff", ZETA)} or \cr
#' \code{list("IF-diff", ZETA, STARTING_RATE)} - Unequal fertility inheritance
#' with parameter \code{ZETA} >= 1 and optionally \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("biased", ZETA)} or \cr
#' \code{list("biased", ZETA, STARTING_RATE)} - Biased speciation
#' with parameter \code{ZETA} >=0 and <=1 and optionally
#' \code{STARTING_RATE} > 0 (default = 1). \cr
#' - \code{list("ASB", ZETA)} or \cr
#' \code{list("ASB", ZETA, STARTING_RATE)} - Age-step-based fertility
#' with parameter \code{ZETA} > 0 and optionally \code{STARTING_RATE} > 0
#' (default = 1). \cr
#' - \code{list("simpleBrown_sym", SIGMA)} or \cr
#' \code{list("simpleBrown_sym", SIGMA, STARTING_RATE)} - Symmetric simple
#' Brownian with parameter \code{SIGMA} >= 0 and optionally
#' \code{STARTING_RATE} > 0 (default = 1). \cr
#' - \code{list("simpleBrown_asym", SIGMA)} or \cr
#' \code{list("simpleBrown_asym", SIGMA, STARTING_RATE)} - Asymmetric simple
#' Brownian with parameter \code{SIGMA} >= 0 and optionally
#' \code{STARTING_RATE} > 0 (default = 1). \cr
#' - \code{list("lin-Brown_sym", SIGMA)} or \cr
#' \code{list("lin-Brown_sym", SIGMA, STARTING_RATE, STARTING_TRAIT)} - Sym.
#' punctuated(-intermittent) linear-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr
#' - \code{list("lin-Brown_asym", SIGMA)} or \cr
#' \code{list("lin-Brown_asym", SIGMA, STARTING_RATE, STARTING_TRAIT)} - Asym.
#' punctuated(-intermittent) linear-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr
#' - \code{list("lin-Brown-bounded_sym", SIGMA)} or \cr
#' \code{list("lin-Brown-bounded_sym", SIGMA, STARTING_RATE, STARTING_TRAIT)} -
#' Bounded sym. punctuated(-intermittent) linear-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr
#' - \code{list("lin-Brown-bounded_asym", SIGMA)} or \cr
#' \code{list("lin-Brown-bounded_asym", SIGMA, STARTING_RATE, STARTING_TRAIT)} -
#' Bounded asym. punctuated(-intermittent) linear-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr
#' - \code{list("log-Brown_sym", SIGMA)} or \cr
#' \code{list("log-Brown_sym", SIGMA, STARTING_RATE, STARTING_TRAIT)} - Sym.
#' punctuated(-intermittent) log-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr
#' - \code{list("log-Brown_asym", SIGMA)} or \cr
#' \code{list("log-Brown_asym", SIGMA, STARTING_RATE, STARTING_TRAIT)} - Asym.
#' punctuated(-intermittent) log-Brownian with parameter
#' \code{SIGMA} (vector with 2 values >=0) and optionally
#' \code{STARTING_RATE} > 0 (default = 1) and \code{STARTING_TRAIT}
#' (default = 10). \cr\cr
#' More information on each model and their parameters can be found in the
#' description of each model, accessible with \code{?genYuleTree},
#'  \code{?genPDATree}, \code{?genETMTree}, \code{?genAldousBetaTree},
#'  \code{?genFordsAlphaTree}, \cr \code{?genBirthDeathTree},
#'  \code{?genAltBirthDeathTree}, \code{?genGrowTree}.
#'
#' @return \code{genTrees}   If \code{Ntrees} is 1, then a single tree of
#' class \code{phylo} is returned.
#' If \code{Ntrees} is larger than 1, a list of class
#' \code{multiPhylo} containing the trees of class \code{phylo} is returned.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic
#'  models: Supplementary material. https://tinyurl.com/278cwdh8, 2024.
#'
#' @export
#' @rdname generateTrees
#'
#' @examples
#' genTrees(n = 5, Ntrees = 2, tm = list("aldous", 1))
#' genTrees(n = 5, tm = "pda")
genTrees <- function(n, Ntrees = 1L, tm) {
  if (n < 2 || n %% 1 != 0) {
    stop(paste(
      "A tree must have at least 2 leaves, i.e., n>=2 and n must be",
      "an integer"
    ))
  }
  if (Ntrees < 1 || Ntrees %% 1 != 0) {
    stop("Ntrees must be a positive integer")
  }

  # Generate trees under specified tree model.
  if (tm[[1]] == "yule") { # -------------------------------------------- Yule model
    return(.wrap_trees(Ntrees, function() genYuleTree(n = n)))
  } else if (tm[[1]] == "pda") { # -------------------------------------- PDA model
    return(.wrap_trees(Ntrees, function() genPDATree(n = n)))
  } else if (tm[[1]] == "etm") { # -------------------------------------- ETM model
    return(.wrap_trees(Ntrees, function() genETMTree(n = n)))
  } else if (tm[[1]] == "aldous") { # -------------------------------- Aldous model
    return(.wrap_trees(Ntrees, function() genAldousBetaTree(n = n, BETA = tm[[2]])))
  } else if (tm[[1]] == "ford") { # ---------------------------------- Ford model
    return(.wrap_trees(Ntrees, function() genFordsAlphaTree(n = n, ALPHA = tm[[2]])))
  } else if (tm[[1]] %in% c(
    "DCO_sym", "DCO_asym", "IF_sym", "IF_asym",
    "IF-diff", "biased", "ASB"
  )) { # ------ Zeta models
    if (length(tm) < 3) tm <- list(tm[[1]], tm[[2]], 1)
    return(.wrap_trees(Ntrees, function() {
      genGrowTree(
        n = n, ZETA = tm[[2]], STARTING_RATE = tm[[3]],
        use_built_in = tm[[1]]
      )
    }))
  } else if (tm[[1]] %in% c("simpleBrown_sym", "simpleBrown_asym")) { # SimBrown
    if (length(tm) < 3) tm <- list(tm[[1]], tm[[2]], 1)
    return(.wrap_trees(Ntrees, function() {
      genGrowTree(
        n = n, SIGMA = tm[[2]], STARTING_RATE = tm[[3]],
        use_built_in = tm[[1]]
      )
    }))
  } else if (tm[[1]] %in% c(
    "lin-Brown_sym", "lin-Brown_asym",
    "lin-Brown-bounded_sym", "lin-Brown-bounded_asym",
    "log-Brown_sym", "log-Brown_asym"
  )) { # -- Brownians
    if (length(tm) < 3) tm <- list(tm[[1]], tm[[2]], 1, 10)
    return(.wrap_trees(Ntrees, function() {
      genGrowTree(
        n = n, SIGMA = tm[[2]], STARTING_RATE = tm[[3]],
        STARTING_TRAIT = tm[[4]], use_built_in = tm[[1]]
      )
    }))
  } else if (tm[[1]] == "alt-birth-death") { # --- alternative birth-death model
    if (length(tm) < 4) tm <- list(tm[[1]], tm[[2]], tm[[3]], 5)
    return(.wrap_trees(Ntrees, function() {
      genAltBirthDeathTree(
        n = n, BIRTHRATE = tm[[2]], DEATHRATE = tm[[3]],
        TRIES = tm[[4]]
      )
    }))
  } else if (tm[[1]] == "density") { # --------------------------------- density
    if (length(tm) < 4) tm <- list(tm[[1]], tm[[2]], tm[[3]], 5, 0.01)
    return(.wrap_trees(Ntrees, function() {
      genDensityTree(
        n = n, BIRTHRATE = tm[[2]], EQUILIB = tm[[3]],
        TRIES = tm[[4]], TIMEperTRY = tm[[5]]
      )
    }))
  } else if (tm[[1]] == "BiSSE") { # ------------------------------------- BiSSE
    if (length(tm) < 5) tm <- list(tm[[1]], tm[[2]], tm[[3]], tm[[4]], 5, 0.01)
    return(.wrap_trees(Ntrees, function() {
      genBiSSETree(
        n = n, BIRTHRATES = tm[[2]], DEATHRATES = tm[[3]],
        TRANSRATES = tm[[4]], TRIES = tm[[5]], TIMEperTRY = tm[[6]]
      )
    }))
  } else {
    stop("Unknown tree model")
  }
}
