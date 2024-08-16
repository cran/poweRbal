#' Generation of rooted binary trees under tree growing models (no extinction)
#'
#' \code{genGrowTree} - Generates a rooted binary tree in \code{phylo}
#' format with the given number of \code{n} leaves under a specified
#' discrete-time tree growing model without extinction.
#' These tree growing models act at the leaves by varying their speciation
#' rates according to a parameter \code{ZETA} or variance \code{SIGMA}. They
#' may also depend on so-called trait values of the leaves (e.g., continuous or
#' discrete age, or another numeric trait that affects fitness).\cr
#' You may choose an already built-in model (see \code{use_built_in}) or
#' specify a (new) model by defining how the rates (and optionally traits)
#' change in every time step (see parameters \code{childRates} and
#' \code{otherRates} as well as \code{childTraits} and
#' \code{otherTraits}; see also Table 5 of the supplementary material of
#' the corresponding manuscript).
#'
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.\cr
#' Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
#' the number of leaves must be at least 2 since there must be at
#' least one edge.
#' @param STARTING_RATE Positive numeric value (default = 1) which specifies the
#' initial rate at which the speciation events occur (has only influence on the
#' edge length, not on the tree topology).
#' @param STARTING_TRAIT Numeric value (default = 10) which specifies
#' the initial state of a trait.
#' @param ZETA Constant non-negative numeric value (default = 1) which can
#' influence the speciation rates. Can also be a vector if used as such when
#' defining the functions \code{childRates}, \code{otherRates},
#' \code{childTraits}, and \code{otherTraits}.
#' @param SIGMA Constant positive numeric value (default = 0) which can influence
#' the speciation rates. Can also be a vector if used as such when defining the
#' functions \code{childRates}, \code{otherRates}, \code{childTraits}, and
#' \code{otherTraits}.
#' @param childRates A function that generates two speciation rates for the
#' children emerging from a speciation event based on various factors. \cr
#' Necessary if \code{use_built_in} is not specified.
#' \code{childTraits} works similarly but is executed before
#' \code{childRates}. \cr
#' All available parameters are: \cr
#' - the starting rate \code{sr},  \cr
#' - the starting trait value \code{st},  \cr
#' - the parent's rate \code{pr},  \cr
#' - the parent's trait value \code{pt},  \cr
#' - the children's trait values \code{ct} (vector \code{ct[1]} and
#' \code{ct[2]}), \cr
#' - the parameters zeta \code{ze}  \cr
#' - and sigma \code{si}. \cr\cr
#' All parameters have to appear in the function definition but not
#' necessarily in the body of the function. Trait values are NA, if
#' \code{childTraits} and \code{otherTraits} is not given.\cr
#' Example: \cr
#' \code{function (sr, st, pr, pt, ct, ze, si) return(c(pr*ze,
#' pr*(1-ze)))} for biased speciation.
#' @param childTraits An optional function (default = NULL) that generates two
#' trait values for the children emerging from a speciation event based on
#' various factors. \cr
#' See \code{childRates} for available parameters (except \code{ct}) and
#' explanations. Not necessary; is only applied if not NULL.\cr
#' Example: \cr
#' \code{function (sr, st, pr, pt, ze, si)
#' return(c(0, 0))} for age.
#' @param otherRates A function that generates a new speciation rate for all
#' leaves not affected by the speciation event (all but parent and children)
#' based on various factors. The function is applied after the speciation event,
#' i.e., after \code{childRates/Traits}.\cr
#' Necessary if \code{use_built_in} is not specified.
#' \code{otherTraits} works similarly. \cr
#' All available parameters are: \cr
#' - the starting rate \code{sr},  \cr
#' - the starting trait value \code{st},  \cr
#' - the leaf's old rate \code{or},  \cr
#' - the leaf's old trait value \code{ot}, \cr
#' - the parameters zeta \code{ze}  \cr
#' - and sigma \code{si}. \cr\cr
#' All parameters have to appear in the function definition but not
#' necessarily in the body of the function. Trait values are NA, if
#' \code{childTraits} and \code{otherTraits} is not given. \cr
#' Example: \cr
#' \code{function (sr, st, or, ot, ze, si) return(or*ze)} for
#' age-step-based fertility.
#' @param otherTraits  An optional function (default = NULL) that generates a new
#' trait value for all leaves not affected by the speciation event (all but
#' parent and children) based on various factors. \cr
#' See \code{otherRates} for available parameters and explanations. \cr
#' Not necessary; is only applied if not NULL.\cr
#' Example:\cr
#' \code{function (sr, st, or, ot, ze, si) return(ot+1)}
#' for discrete age (age in time steps).
#' @param use_built_in Optional (default = NULL): Character specifying which of
#' the already implemented models should be used. Overwrites \code{childRates},
#' \code{otherRates}, \code{childTraits}, and \code{otherTraits}.\cr
#' Here is a list of available models with their (abbreviated) underlying
#' functions given in parentheses  (in order \code{childRates},
#' \code{otherRates}; then \code{childTraits} and \code{otherTraits}
#' if necessary): \cr
#' - "DCO_sym": Symmetric direct-children-only, \code{ZETA}>0
#' (c(sr ze, sr ze), sr) \cr
#' - "DCO_asym": Asymmetric direct-children-only, \code{ZETA}>0
#' (c(sz, pr), sr) \cr
#' - "IF_sym": Symmetric inherited fertility, \code{ZETA}>0
#' (c(pr ze, pr ze), or) \cr
#' - "IF_asym": Asymmetric inherited fertility, \code{ZETA}>0
#' (c(pr ze, pr), or) \cr
#' - "IF-diff": Unequal fertility inheritance, \code{ZETA}>=1
#' (c(2 pr ze / (ze+1), 2 pr / (ze+1)), or)\cr
#' - "biased": Biased speciation, \code{ZETA} >=0 and <=1
#' (c(pr ze, pr (1-ze)), or)\cr
#' - "ASB": Age-step-based fertility, \code{ZETA}>0
#' (c(sr, sr), or ze)\cr
#' - "simpleBrown_sym": Symmetric simple Brownian, \code{SIGMA}> =0
#'    (c(max\{pr+ rnorm(1, mean=0, sd=si),1e-100\},
#'       max\{pr+ rnorm(1, mean=0, sd=si),1e-100\}), or) \cr
#' - "simpleBrown_asym": Asymmetric simple Brownian, \code{SIGMA}>=0
#'    (c(max\{pr+ rnorm(1, mean=0, sd=si),1e-100\}, pr), or) \cr
#' - "lin-Brown_sym": Sym. punctuated(-intermittent) linear-Brownian,
#' \code{SIGMA} vector with two values >=0 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])),
#'      10^(log(ct\[2\])+ rnorm(1, mean=0, sd=si\[1\]))), or; \cr
#'    c(max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\},
#'      max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\}), ot) \cr
#' - "lin-Brown_asym": Asym. punctuated(-intermittent) linear-Brownian,
#' \code{SIGMA} vector with two values >=0 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])), pr), or; \cr
#'    c(max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\}, pt, ot) \cr
#' - "lin-Brown-bounded_sym": Bounded sym. punctuated(-intermittent)
#' linear-Brownian,  \code{SIGMA} vector with two values >=0,
#' \code{STARTING_TRAIT} is automatically set to 10 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])),
#'      10^(log(ct\[2\])+ rnorm(1, mean=0, sd=si\[1\]))), or; \cr
#'    c(min\{max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\},20\},
#'      min\{max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\},20\}), ot) \cr
#' - "lin-Brown-bounded_asym": Bounded asym. punctuated(-intermittent)
#' linear-Brownian, \code{SIGMA} vector with two values >=0 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])), pr), or; \cr
#'    c(min\{max\{pt + rnorm(1, mean=0, sd=si\[2\]),1e-100\},20\}, pt), ot) \cr
#' - "log-Brown_sym": Sym. punctuated(-intermittent) log-Brownian,
#' \code{SIGMA} vector with two values >=0 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])),
#'      10^(log(ct\[2\])+ rnorm(1, mean=0, sd=si\[1\]))), or; \cr
#'    c(10^(log(pt)+ rnorm(1, mean=0, sd=si\[2\])),
#'      10^(log(pt)+ rnorm(1, mean=0, sd=si\[2\]))), ot) \cr
#' - "log-Brown_asym": Asym. punctuated(-intermittent) log-Brownian,
#' \code{SIGMA} vector with two values >=0 \cr
#'   (c(10^(log(ct\[1\])+ rnorm(1, mean=0, sd=si\[1\])), pr), or; \cr
#'    10^(c(log(pt)+ rnorm(1, mean=0, sd=si\[2\])), pt), ot)
#'
#' @return \code{genGrowTree} A single tree of class \code{phylo} is
#' returned.
#'
#' @references
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic models.
#'  arXiv:2406.05185, 2024.
#'  - S. J. Kersting, K. Wicke, and M. Fischer. Tree balance in phylogenetic
#'  models: Supplementary material. https://tinyurl.com/278cwdh8, 2024.
#'  - M. G. B. Blum and O. Francois. On statistical tests of phylogenetic tree
#'  imbalance: the Sackin and other indices revisited. Mathematical Biosciences,
#'  195(2):141–153, 2005.
#'  - S. B. Heard. Patterns in phylogenetic tree balance with variable and
#'  evolving speciation rates. Evolution, 50(6):2141–2148, 1996.
#'  - S. J. Kersting. Genetic programming as a means for generating improved
#'  tree balance indices (Master’s thesis, University of Greifswald), 2020.
#'  - M. Kirkpatrick and M. Slatkin. Searching for evolutionary patterns in the
#'  shape of a phylogenetic tree. Evolution, 47(4):1171–1181, 1993.
#'
#'
#' @export
#' @rdname tmGrowNoExtDT
#'
#' @examples
#' genGrowTree(n = 5, use_built_in = "IF_sym", ZETA = 2)
genGrowTree <- function(n,
                        STARTING_RATE = 1, STARTING_TRAIT = 10,
                        ZETA = 1, SIGMA = 0,
                        childRates, otherRates,
                        childTraits = NULL, otherTraits = NULL,
                        use_built_in = NULL){
  if(n < 2 || n%%1!=0){
    stop(paste("A tree must have at least 2 leaves, i.e., n>=2 and n must be",
               "an integer."))
  }
  if(STARTING_RATE<=0){
    stop(paste("The speciation starting rate must be >0."))
  }
  if(!is.null(use_built_in)){
    if (use_built_in == "DCO_sym") {
      if(ZETA<=0){
        stop("The parameter ZETA must be >0 for DCO (sym).")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si) return(c(sr*ze, sr*ze))
      otherRates <- function(sr, st, or, ot, ze, si) return(sr)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "DCO_asym") {
      if(ZETA<=0){
        stop("The parameter ZETA must be >0 for DCO (asym).")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si) return(c(sr*ze, pr))
      otherRates <- function(sr, st, or, ot, ze, si) return(sr)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "IF_sym") {
      if(ZETA<=0){
        stop("The parameter ZETA must be >0 for IF (sym).")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si) return(c(pr*ze, pr*ze))
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "IF_asym") {
      if(ZETA<=0){
        stop("The parameter ZETA must be >0 for IF (asym).")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si) return(c(pr*ze, pr))
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "IF-diff") {
      if(ZETA<1){
        stop("The parameter ZETA must be >=1 for IF-diff.")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(2*pr*ze/(ze+1), 2*pr/(ze+1)))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "biased") {
      if(ZETA<0 || ZETA >1){
        stop("The parameter ZETA must be in [0,1] for biased speciation.")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(pr*ze, pr*(1-ze)))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "ASB") {
      if(ZETA<=0){
        stop("The parameter ZETA must be >0 for age-step-based fertility.")
      }
      SIGMA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si) return(c(sr, sr))
      otherRates <- function(sr, st, or, ot, ze, si) return(or * ze)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "simpleBrown_sym") {
      if(SIGMA<0){
        stop("The parameter SIGMA must be >0 for sym. simple Brownian.")
      }
      ZETA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(max(c(1e-100, pr + stats::rnorm(1, mean = 0, sd = si))),
                 max(c(1e-100, pr + stats::rnorm(1, mean = 0, sd = si)))))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "simpleBrown_asym") {
      if(SIGMA<0){
        stop("The parameter SIGMA must be >=0 for asym. simple Brownian.")
      }
      ZETA <- NULL; STARTING_TRAIT <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(max(c(1e-100, pr + stats::rnorm(1, mean = 0, sd = si))), pr))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- NULL
      otherTraits <- NULL
    } else if (use_built_in == "lin-Brown_sym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for sym. linear-Brownian."))
      }
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean = 0, sd = si[1])),
                 10^(log10(ct[2]) + stats::rnorm(1, mean = 0, sd = si[1]))))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(max(c(pt + stats::rnorm(1, mean = 0, sd = si[2]),1e-100)),
                 max(c(pt + stats::rnorm(1, mean = 0, sd = si[2]),1e-100))))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else if (use_built_in == "lin-Brown_asym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for sym. linear-Brownian."))
      }
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean=0, sd = si[1])), pr))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(max(c(pt + stats::rnorm(1, mean = 0, sd = si[2]),1e-100)), pt))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else if (use_built_in == "lin-Brown-bounded_sym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for bounded sym. linear-Brownian."))
      }
      STARTING_TRAIT <- 10
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean = 0, sd = si[1])),
                 10^(log10(ct[2]) + stats::rnorm(1, mean = 0, sd = si[1]))))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(min(c(max(c(pt + stats::rnorm(1, mean = 0, sd = si[2]),0)),20)),
                 min(c(max(c(pt + stats::rnorm(1, mean = 0, sd = si[2]),0)),20))))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else if (use_built_in == "lin-Brown-bounded_asym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for bounded sym. linear-Brownian."))
      }
      STARTING_TRAIT <- 10
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean = 0, sd = si[1])),
                 pr))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(min(c(max(c(pt + stats::rnorm(1, mean = 0, sd =si[2]),0)),20)),
                 pt))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else if (use_built_in == "log-Brown_sym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for sym. log-Brownian."))
      }
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean = 0, sd = si[1])),
                 10^(log10(ct[2]) + stats::rnorm(1, mean = 0, sd = si[1]))))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(10^(log10(pt) + stats::rnorm(1, mean = 0, sd = si[2])),
                 10^(log10(pt) + stats::rnorm(1, mean = 0, sd = si[2]))))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else if (use_built_in == "log-Brown_asym") {
      if(length(SIGMA)!= 2 || SIGMA[1]<0 || SIGMA[2]<0){
        stop(paste("The parameter SIGMA must be a vector of two values >=0",
                   "for sym. log-Brownian."))
      }
      ZETA <- NULL
      childRates <- function(sr, st, pr, pt, ct, ze, si){
        return(c(10^(log10(ct[1]) + stats::rnorm(1, mean = 0, sd = si[1])),
                 pr))}
      otherRates <- function(sr, st, or, ot, ze, si) return(or)
      childTraits <- function(sr, st, pr, pt, ze, si){
        return(c(10^(log10(pt) + stats::rnorm(1, mean = 0, sd = si[2])), pt))}
      otherTraits <-  function(sr, st, or, ot, ze, si) return(ot)
    } else {
      stop("Unknown tree growing model.")
    }
  }
  # Create the edge matrix -----------------------------------------------------
  m <- matrix(rep(NA,(2*n-2)*2), nrow = 2*n-2, ncol = 2)
  # Initialize vector for current leaves and their rates -----------------------
  curr_leaves <- 1 # Vector of the current leaves (start: only one node)
  curr_rates <- STARTING_RATE
  # Initialize vector for current trait values ---------------------------------
  curr_traits <- rep(NA, 2*n-1) # set to NA if traits are not used
  TRAITS_ARE_USED <- FALSE
  if(!is.null(childTraits) && !is.null(otherTraits)){
    TRAITS_ARE_USED <- TRUE
    curr_traits <- STARTING_TRAIT
  }
  # Do (n-1) speciation steps --------------------------------------------------
  for(i in 1:(n-1)) {
    leaf_index <- sample(1:length(curr_leaves), size = 1, replace = F,
                         prob = curr_rates) # Choose leaf for speciation
    # New numbers for children
    c_1 <- 2*i
    c_2 <- 1+2*i
    # Fill out matrix row by row (edges: parent->child)
    m[2*i-1,] <- c(curr_leaves[leaf_index], c_1)
    m[2*i,] <- c(curr_leaves[leaf_index], c_2)
    # Remove parent from current leaves and insert children.
    curr_leaves <- c(curr_leaves[-leaf_index], c_1, c_2)
    # Update the traits if necessary, remove parent value and insert children's
    # trait values.
    if(TRAITS_ARE_USED){
      new_child_traits <- childTraits(sr = STARTING_RATE, st = STARTING_TRAIT,
                                      pr = curr_rates[leaf_index],
                                      pt = curr_traits[leaf_index],
                                      ze = ZETA, si = SIGMA)
      if(i>1){
        new_other_traits <- sapply((1:i)[-leaf_index],
                                   function(x) {
                                     otherTraits(sr = STARTING_RATE,
                                                 st = STARTING_TRAIT,
                                                 or = curr_rates[x],
                                                 ot = curr_traits[x],
                                                 ze = ZETA, si = SIGMA)})
      } else {
        new_other_traits <- NULL
      }
      curr_traits <- c(new_other_traits, new_child_traits)
    }
    # Update the rates, remove parent rate and insert children's rates.
    new_child_rates <- childRates(sr = STARTING_RATE, st = STARTING_TRAIT,
                                  pr = curr_rates[leaf_index],
                                  pt = curr_traits[leaf_index],
                                  ct = curr_traits[i:(i+1)],
                                  ze = ZETA, si = SIGMA)
    if(i>1){
      new_other_rates <- sapply((1:i)[-leaf_index],
                                 function(x) {
                                   otherRates(sr = STARTING_RATE,
                                               st = STARTING_TRAIT,
                                               or = curr_rates[x],
                                               ot = curr_traits[x],
                                               ze = ZETA, si = SIGMA)})
    } else {
      new_other_rates <- NULL
    }
    curr_rates <- c(new_other_rates, new_child_rates)
  }
  # Create the phylo object and enumerate cladewise ----------------------------
  phy <- list(edge = m, tip.label = paste("t", sample.int(n,n), sep = ""),
              Nnode = as.integer(n-1))
  attr(phy, "class") <- "phylo"
  return(enum2cladewise(phy, root = 1))
}
