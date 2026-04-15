#' Functions for computing the region of acceptance
#'
#' \code{getAccRegion} - Computes the region of acceptance based on quantiles
#' for a specified level of significance and method.
#'
#' @param tss Vector containing the names (as character) of the tree shape
#' statistics that should be compared. You may either use the short names
#' provided in \code{tssInfo} to use the already included TSS, or use the
#' name of a list object containing similar information as the entries in
#' \code{tssInfo}. Example:\cr
#' Use \code{"new_tss"} as the name for the list object
#' \code{new_tss} containing at least the function
#' \code{new_tss$func = function(tree){...}},
#' and optionally also the information \code{new_tss$short},
#' \code{new_tss$simple}, \code{new_tss$name}, \code{new_tss$type},
#' \code{new_tss$only_binary}, and \code{new_tss$safe_n}.
#' @param null_model The null model that is to be used to determine the power
#' of the tree shape statistics. In general, it must be a function that
#' produces rooted binary trees in \code{phylo} format. \cr
#' If the respective model is included in this
#' package, then specify the model and its parameters by using a character
#' or list. Available are all options listed under parameter \code{tm} in
#' the documentation of function \code{genTrees} (type \code{?genTrees}).\cr
#' If you want to include your own tree model, then use the
#' name of a list object containing the function (with the two input parameters
#' \code{n} and \code{Ntrees}). Example: \cr
#' Use \code{"new_tm"} for the list object
#' \code{new_tm <- list(func = function(n, Ntrees){...})}.
#' @param n Integer value that specifies the desired number of leaves, i.e.,
#' vertices with in-degree 1 and out-degree 0.
#' @param distribs Determines how the distributions (and with that the
#' bounds of the critical region) are computed. Available are: \cr
#' - "exact_if_possible" (default): Tries to compute the exact distribution
#' under the null model if possible. Currently, this is only implemented for
#' \code{null_model = "yule"}, \code{"pda"}, or \code{"etm"}, and
#' \code{n}<=20. In all other cases the distribution is approximated
#' by sampling \code{N_null} many trees under the null model as in the
#' option "sampled" below. \cr
#' - "sampled": \code{N_null} many trees are sampled under the
#' null model to approximate the distribution.
#' @param N_null Sample size (integer >=10) if distributions are sampled
#' (default = 10000L).
#' @param N_alt Sample size (integer >=10) for the alternative models to
#' estimate the power (default = 1000L). Only needed here if the
#' \code{test_type} is "two-tailed-unbiased".
#' @param N_intervals Number (integer >=3, default = 1000L) of different
#' quantile/cut-off pairs investigated as potential bounds of the region of
#' acceptance. This parameter is only necessary if the \code{test_type} is
#' "two-tailed-unbiased".
#' @param test_type Determines the method. Available are: \cr
#' - "two-tailed" (default): The lower and upper bound of the region of
#' acceptance are determined based on the (empirical) distribution function
#' such that P(TSS < lower bound) <= \code{sig_lvl}/2 and
#' P(TSS > upper bound) <= \code{sig_lvl}/2. See parameter \code{correction}
#' for specifying how conservative the test should be: the null
#' hypothesis can either be rejected only if the values are strictly outside of
#' this region of acceptance (can be too conservative) or it can also be
#' rejected (with certain probabilities) if the value equals the lower or
#' upper bound.\cr
#' - "two-tailed-unbiased": Experimental - Use with caution!\cr
#' The region of acceptance is optimized to yield an unbiased test, i.e., a test
#' that identifies non-null models with a probability of at least
#' \code{sig_lvl}.
#' The region of acceptance is determined similar to the default method.
#' However, it need not be symmetrical, i.e., not necessarily
#' cutting off \code{sig_lvl}/2 on both sides. Also see parameter
#' \code{correction} for specifying how conservative the test should be.
#' @param correction Specifies the desired correction method.
#' Available are: \cr
#' - "small-sample" (default): This method tries to ensure that the critical
#' region, i.e., the range of values for which the null hypothesis is rejected,
#' is as close to \code{sig_lvl} as possible (compared with "none" below, which
#' can be too conservative). The idea is that the null hypothesis is also
#' rejected with certain probabilities if the value matches a bound of the
#' region of acceptance. \cr
#' - "none": No correction method is applied. With that the test might be
#' slightly too conservative as the null hypothesis is maintained if the values
#' are >= the lower and <= the upper bound.
#' @param sig_lvl Level of significance (default=0.05, must be >0 and <1).
#'
#' @return \code{getAccRegion} Numeric matrix (one row per TSS) with four
#' columns: The first two columns contain the interval limits of the region
#' of acceptance, i.e., we reject the null hypothesis for values strictly
#' outside of this interval. The third and fourth columns contain the
#' probabilities to reject the null hypothesis if values equal the lower or
#' upper bound, respectively.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' getAccRegion(tss = c("Sackin", "Colless", "B1I"), n = 6L)
#' getAccRegion(
#'   tss = c("Sackin", "Colless", "B1I"), n = 6L, null_model = "etm",
#'   N_null = 20L, correction = "none", distribs = "sampled"
#' )
#' getAccRegion(
#'   tss = c("Sackin", "Colless", "B1I"), n = 6L, N_null = 20L,
#'   test_type = "two-tailed-unbiased", N_intervals = 5L,
#'   N_alt = 10L
#' )
getAccRegion <- function(tss,
                         null_model = "yule", n, distribs = "exact_if_possible",
                         N_null = 10000L, N_alt = 1000L, N_intervals = 1000L,
                         test_type = "two-tailed", correction = "small-sample",
                         sig_lvl = 0.05) {
  if (!is.list(null_model) && !is.character(null_model)) {
    stop("null_model must be a character string or a list")
  }
  if (!is.numeric(n) || n < 2 || n %% 1 != 0) {
    stop("n must be an integer >= 2")
  }
  if (distribs == "exact_if_possible") {
    if (n <= 20 && (null_model[[1]] == "yule" || null_model[[1]] == "pda" ||
      null_model[[1]] == "etm")) {
      return(getAccRegion_exact(
        tss = tss, null_model = null_model, n = n,
        N_alt = N_alt, N_intervals = N_intervals,
        test_type = test_type,
        correction = correction, sig_lvl = sig_lvl
      ))
    } else {
      message(paste(
        "Exact computation of the distribution not available",
        "for the given choice of null model and n.",
        "Parameter distribs = 'sampled' used instead."
      ))
      return(getAccRegion_sampled(
        tss = tss, null_model = null_model, n = n,
        N_null = N_null, N_alt = N_alt,
        N_intervals = N_intervals,
        test_type = test_type,
        correction = correction, sig_lvl = sig_lvl
      ))
    }
  } else if (distribs == "sampled") {
    return(getAccRegion_sampled(
      tss = tss, null_model = null_model, n = n,
      N_null = N_null, N_alt = N_alt,
      N_intervals = N_intervals,
      test_type = test_type,
      correction = correction, sig_lvl = sig_lvl
    ))
  } else {
    stop("Unknown distribs method")
  }
}
#' Functions for computing the region of acceptance
#'
#' \code{getAccRegion_sampled} - Computes a sampling-based region of acceptance
#' for the given null model based on quantiles for a specified level of
#' significance and method.
#'
#' @return \code{getAccRegion_sampled} Numeric matrix (one row per TSS) with
#' four columns - similar as \code{getAccRegion}.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' getAccRegion_sampled(
#'   tss = c("Sackin", "Colless", "B1I"), n = 6L,
#'   N_null = 20L, correction = "none"
#' )
getAccRegion_sampled <- function(tss, null_model = "yule", n, N_null,
                                 N_alt = 1000L, N_intervals = 1000L,
                                 test_type = "two-tailed",
                                 correction = "small-sample", sig_lvl = 0.05) {
  # --- Compute the TSS values of trees sampled under the null model.
  null_vals <- getTSSdata(
    tss = tss, n = n, Ntrees = N_null,
    tm = null_model
  )
  if (is.null(dim(null_vals))) {
    null_vals <- matrix(null_vals,
      nrow = 1, byrow = TRUE,
      dimnames = list(tss[[1]], NULL)
    )
  }
  #  --------- Compute the unique values and their probabilities.
  all_unique_null_vals <- list()
  all_unique_null_probs <- list()
  for (i_tss in seq_len(nrow(null_vals))) {
    tbl_i_tss <- table(null_vals[i_tss, ])
    all_unique_null_vals[[i_tss]] <- as.numeric(names(tbl_i_tss))
    all_unique_null_probs[[i_tss]] <- unname(tbl_i_tss) / sum(tbl_i_tss)
  }
  #  --------- Calculate Acceptance Intervals.
  acc_regions <- matrix(NA,
    nrow = length(tss), ncol = 4,
    dimnames = list(
      rownames(null_vals),
      c(
        "lower_limit", "upper_limit",
        "lower_rej_prob", "upper_rej_prob"
      )
    )
  )
  if (test_type == "two-tailed") {
    for (i_tss in seq_along(tss)) {
      acc_regions[i_tss, ] <- computeAccRegion(
        unique_null_vals = all_unique_null_vals[[i_tss]],
        unique_null_probs = all_unique_null_probs[[i_tss]],
        correction = correction,
        cutoff_left = sig_lvl / 2,
        cutoff_right = sig_lvl / 2
      )
    }
  } else if (test_type == "two-tailed-unbiased") {
    acc_regions <- .getAccRegion_unbiased(
      tss = tss,
      all_unique_null_vals = all_unique_null_vals,
      all_unique_null_probs = all_unique_null_probs,
      n, N_alt = N_alt, N_intervals = N_intervals,
      correction = correction, sig_lvl = sig_lvl
    )
  } else {
    stop("Unknown test type")
  }
  return(acc_regions)
}
#' Functions for computing the region of acceptance
#'
#' \code{getAccRegion_exact} - Computes the exact region of acceptance for the
#' given null model based on quantiles for a specified level of significance
#' and method. Currently, this is only implemented for
#' \code{null_model = "yule"}, \code{"pda"}, or \code{"etm"}, and \code{n}<=20.
#'
#' The restriction to \code{n <= 20} arises from the rapid growth of the
#' Wedderburn-Etherington numbers (the count of distinct unlabelled rooted
#' binary tree topologies with \code{n} leaves), which makes exhaustive
#' enumeration infeasible beyond this limit. The implementation iterates over
#' all topologies using the Furnas rank rather than storing them
#' simultaneously, which keeps memory usage tractable for \code{n <= 20}.
#' For larger \code{n}, use \code{\link{getAccRegion}} with
#' \code{distribs = "sampled"}.
#'
#' @return \code{getAccRegion_exact} Numeric matrix (one row per TSS) with
#' four columns - similar as \code{getAccRegion}.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' getAccRegion_exact(
#'   tss = c("Sackin", "Colless", "B1I"),
#'   null_model = "etm", n = 8L
#' )
getAccRegion_exact <- function(tss, null_model = "yule", n,
                               N_alt = 1000L, N_intervals = 1000L,
                               test_type = "two-tailed",
                               correction = "small-sample", sig_lvl = 0.05) {
  if (n < 2 || n > 20) {
    stop("Not computing the exact region of acceptance for n > 20")
  }
  if (!null_model %in% c("yule", "pda", "etm")) {
    stop("Cannot compute exact region of acceptance for this null model")
  }
  # --- Compute the exact null distribution:
  n_tops <- treebalance::we_eth(n)
  # Process each topology individually to avoid holding all n_tops phylo objects
  # in memory simultaneously. The memo cache is cleared on exit so that cached
  # topology objects do not persist beyond this call.
  on.exit(memoise::forget(.furnasI_inv), add = TRUE)
  # Pre-allocate storage for TSS values and topology probabilities.
  null_vals <- matrix(NA_real_,
    nrow = length(tss), ncol = n_tops,
    dimnames = list(tss, NULL)
  )
  prob_top_n <- rep(NA_real_, n_tops)
  for (k in seq_len(n_tops)) {
    top_k <- .furnasI_inv(rank = k, n = n)
    # TSS values for this topology
    for (i_tss in seq_along(tss)) {
      if (!is.null(tssInfo[[tss[i_tss]]])) {
        null_vals[i_tss, k] <- tssInfo[[tss[i_tss]]]$func(top_k)
      } else if (exists(tss[i_tss])) {
        tss_obj <- get(tss[i_tss])
        if (!is.null(tss_obj$func)) {
          null_vals[i_tss, k] <- tss_obj$func(top_k)
          if (k == 1L && !is.null(tss_obj$short)) {
            rownames(null_vals)[i_tss] <- tss_obj$short
          }
        } else {
          stop(paste0("No function provided for TSS ", tss[i_tss]))
        }
      } else {
        stop(paste("Unknown TSS:", tss[i_tss]))
      }
    }
    # Topology probability under the null model
    if (null_model == "etm") {
      prob_top_n[k] <- 1 / n_tops
    } else {
      # Number of phylogenetic trees with topology top_k
      phyl_k <- factorial(n) / 2^(n - 1 - symNodesI(top_k))
      if (null_model == "yule") {
        # Clade sizes minus 1
        clades_m1 <- treebalance::get.subtreesize(top_k) - 1
        prob_top_n[k] <- (2^(n - 1) / factorial(n)) *
          prod(1 / clades_m1[clades_m1 != 0]) * phyl_k
      } else { # pda: accumulate unnormalized, normalize below
        prob_top_n[k] <- phyl_k
      }
    }
  }
  if (null_model == "pda") {
    prob_top_n <- prob_top_n / sum(prob_top_n) # normalize
  }
  # --------- Get the unique values and their probabilities:
  acc_regions <- matrix(NA,
    nrow = length(tss), ncol = 4,
    dimnames = list(
      rownames(null_vals),
      c(
        "lower_limit", "upper_limit",
        "lower_rej_prob", "upper_rej_prob"
      )
    )
  )
  all_unique_null_vals <- list()
  all_unique_null_probs <- list()
  for (i_tss in seq_along(tss)) {
    unique_ranks <- rank(null_vals[i_tss, ], ties.method = "min")
    unique_null_vals <- rep(NA_real_, max(unique_ranks))
    unique_null_probs <- rep(0, max(unique_ranks))
    for (i in 1:ncol(null_vals)) {
      unique_null_vals[unique_ranks[i]] <- null_vals[i_tss, i]
      unique_null_probs[unique_ranks[i]] <- unique_null_probs[unique_ranks[i]] +
        prob_top_n[i]
    }
    filled <- !is.na(unique_null_vals)
    all_unique_null_probs[[i_tss]] <- unique_null_probs[filled]
    all_unique_null_vals[[i_tss]] <- unique_null_vals[filled]
  }
  if (test_type == "two-tailed") {
    for (i_tss in seq_along(tss)) {
      acc_regions[i_tss, ] <- computeAccRegion(
        unique_null_vals = all_unique_null_vals[[i_tss]],
        unique_null_probs = all_unique_null_probs[[i_tss]],
        correction = correction,
        cutoff_left = sig_lvl / 2,
        cutoff_right = sig_lvl / 2
      )
    }
  } else if (test_type == "two-tailed-unbiased") {
    acc_regions <- .getAccRegion_unbiased(
      tss = tss,
      all_unique_null_vals = all_unique_null_vals,
      all_unique_null_probs = all_unique_null_probs,
      n, N_alt = N_alt, N_intervals = N_intervals,
      correction = correction, sig_lvl = sig_lvl
    )
  } else {
    stop("Unknown test type")
  }
  return(acc_regions)
}
.furnasI_inv <- memoise::memoise(function(rank, n) {
  if (rank < 1 || rank %% 1 != 0) {
    stop("Tree cannot be calculated, because rank is not valid")
  }
  if (n < 1 || n %% 1 != 0) {
    stop(paste(
      "Tree cannot be calculated, because number of leaves",
      "is no positive integer"
    ))
  }
  if (rank > treebalance::we_eth(n)) {
    stop(paste(
      "Tree cannot be calculated, because rank", rank,
      "is larger than the available number of trees",
      treebalance::we_eth(n), "for n =", n
    ))
  }
  if (n == 1) {
    return(ape::read.tree(text = "();"))
  }
  we_nums_mult <- sapply(1:ceiling(n / 2), function(x) treebalance::we_eth(x)) *
    rev(sapply(floor(n / 2):(n - 1), function(x) treebalance::we_eth(x)))
  rsums <- cumsum(we_nums_mult)
  alpha <- min(which(rsums >= rank))
  rsums_alpha1 <- ifelse(alpha > 1, rsums[alpha - 1], 0)
  beta <- n - alpha
  if (alpha < beta) {
    b_temp <- (rank - rsums_alpha1) %% treebalance::we_eth(beta)
    a_temp <- (rank - rsums_alpha1 - b_temp) / treebalance::we_eth(beta)
    if (b_temp > 0) {
      r_alpha <- a_temp + 1
      r_beta <- b_temp
    } else if (b_temp == 0) {
      r_alpha <- a_temp
      r_beta <- treebalance::we_eth(beta)
    }
  } else if (alpha == beta) {
    temp_val <- 2 * (treebalance::we_eth(beta) - rank + rsums_alpha1) +
      (1 - 2 * treebalance::we_eth(beta))^2 / 4
    m <- max(c(0, ceiling(treebalance::we_eth(beta) - 1 / 2 - sqrt(temp_val))))
    r_alpha <- m + 1
    r_beta <- rank - rsums_alpha1 - (r_alpha - 1) * treebalance::we_eth(beta) +
      (r_alpha - 2) * (r_alpha - 1) / 2 + r_alpha - 1
  }
  tL <- .furnasI_inv(rank = r_alpha, n = alpha)
  tR <- .furnasI_inv(rank = r_beta, n = beta)
  return(treebalance::tree_merge(tL, tR))
})
#' Functions for computing the region of acceptance
#'
#' \code{computeAccRegion} - Computes the bounds of the region of acceptance
#' given the empirical distribution function (specified by the unique values
#' and their probabilities under the null model) for specified cut-offs
#' (e.g., 0.025 on both sides for a symmetric two-tailed test).
#' For values strictly outside of the interval the null hypothesis is
#' rejected. \cr
#' This function also computes the probabilities to
#' reject the null hypothesis if the value equals the lower or upper bound of
#' the region of acceptance. This probability is 0 for correction method
#' "none" and for "small-sample" it ensures that the probability of rejection
#' exactly corresponds with the specified cut-offs.
#'
#' @param unique_null_vals Numeric vector containing all the unique values under
#' the null model.
#' @param unique_null_probs Numeric vector containing the corresponding
#' probabilities of the unique values under the null model.
#' @param cutoff_left Numeric value (>=0, <1) specifying the cut-off of the
#' distribution for the lower bound of the region of acceptance. The sum of
#' the two cut-offs must be <1.
#' @param cutoff_right Numeric value (>=0, <1) specifying the cut-off of the
#' distribution for the upper bound of the region of acceptance. The sum of
#' the two cut-offs must be <1.
#'
#' @return \code{computeAccRegion} Numeric vector with
#' four columns - similar as \code{getAccRegion}.
#'
#' @export
#' @rdname powerRegAcc
#'
#' @examples
#' computeAccRegion(
#'   unique_null_vals = c(1, 2, 3, 4, 5),
#'   unique_null_probs = c(0.1, 0.4, 0.1, 0.2, 0.2),
#'   correction = "small-sample",
#'   cutoff_left = 0.15, cutoff_right = 0.15
#' )
computeAccRegion <- function(unique_null_vals, unique_null_probs, correction,
                             cutoff_left, cutoff_right) {
  if ((cutoff_left < 0 || cutoff_left >= 1) ||
    (cutoff_right < 0 || cutoff_right >= 1) ||
    (cutoff_left + cutoff_right >= 1)) {
    stop("Cannot compute the region of acceptance for the given cut-offs")
  }
  acc_region <- rep(NA, 4)
  unique_null_cumprob <- cumsum(unique_null_probs)
  # ------------ Lower and uppper bound:
  index_low_bound <- max(c(which(unique_null_cumprob <= cutoff_left), 0)) + 1
  index_upp_bound <- min(c(
    which(unique_null_cumprob >= (1 - cutoff_right)),
    length(unique_null_cumprob)
  ))
  acc_region[c(1, 2)] <- unique_null_vals[c(
    index_low_bound,
    index_upp_bound
  )]
  if (acc_region[1] > acc_region[2]) {
    warning(
      "Lower bound > upper bound for the region of acceptance ",
      "(all null values may be identical). Bounds are set equal (lower bound)."
    )
    acc_region[2] <- acc_region[1]
  }
  # ------------ The probabilities of rejecting H0 if value == bound:
  if (correction == "small-sample") {
    acc_region[c(3, 4)] <- c(
      (cutoff_left - (if (index_low_bound == 1) {
        0
      } else {
        unique_null_cumprob[index_low_bound - 1]
      })) /
        unique_null_probs[index_low_bound],
      (cutoff_right - (1 - unique_null_cumprob[index_upp_bound])) /
        unique_null_probs[index_upp_bound]
    )
    # For cutoffs close to 0, there may be small computation errors that
    # result in probabilities < 0.
    neg_probs <- which(acc_region[c(3, 4)] < 0)
    if (length(neg_probs) > 0) {
      warning(paste(
        "Rejection probabilities were negative due to floating-point ",
        "rounding near the significance level boundary; (",
        paste(acc_region[c(3, 4)], collapse = ", "), ") clamped to 0."
      ))
      acc_region[neg_probs + 2] <- 0
    }
  } else if (correction == "none") {
    acc_region[c(3, 4)] <- c(0, 0)
  } else {
    stop("Unknown correction method")
  }
  return(acc_region)
}
#-------------------------------------------------------------------------------
# Sub-functions for .getAccRegion_unbiased()

# Compute the N_intervals candidate acceptance regions for one TSS.
.evaluate_interval_candidates <- function(unique_null_vals, unique_null_probs,
                                          quantile_pairs, N_intervals,
                                          correction) {
  intervals <- matrix(NA,
    nrow = N_intervals, ncol = 4,
    dimnames = list(
      seq(0, quantile_pairs[N_intervals, 1],
        length.out = N_intervals
      ),
      c(
        "lower_limit", "upper_limit",
        "lower_rej_prob", "upper_rej_prob"
      )
    )
  )
  for (q_pair in seq_len(N_intervals)) {
    intervals[q_pair, ] <- computeAccRegion(
      unique_null_vals  = unique_null_vals,
      unique_null_probs = unique_null_probs,
      correction        = correction,
      cutoff_left       = quantile_pairs[q_pair, 1],
      cutoff_right      = 1 - quantile_pairs[q_pair, 2]
    )
  }
  intervals
}

# Compute acceptance probabilities (N_intervals x n_alt_models) for one TSS.
.compute_accept_p <- function(intervals, alt_vals_tss_row,
                              N_intervals, N_alt) {
  n_alt <- length(alt_vals_tss_row)
  accept_p <- matrix(NA, nrow = N_intervals, ncol = n_alt)
  for (q_pair in seq_len(N_intervals)) {
    for (col_atm in seq_len(n_alt)) {
      vals <- alt_vals_tss_row[[col_atm]]
      accept_p[q_pair, col_atm] <-
        (sum(vals > intervals[q_pair, 1] & vals < intervals[q_pair, 2]) +
          sum(vals == intervals[q_pair, 1]) * (1 - intervals[q_pair, 3]) +
          sum(vals == intervals[q_pair, 2]) * (1 - intervals[q_pair, 4])) / N_alt
    }
  }
  accept_p
}

# Select the best interval index from candidates and log messages.
.select_best_interval <- function(intervals, accept_p, unbiased_threshold,
                                  quantile_pairs, tss_name) {
  n_intervals <- nrow(intervals)
  # Perfect candidates: all alt models meet the threshold.
  candidates <- rowSums(accept_p <= unbiased_threshold) == ncol(accept_p)
  # Least-bad fallback if no perfect candidate exists.
  if (!any(candidates)) {
    deviations <- pmax(accept_p - unbiased_threshold, 0)
    candidate_dev <- rowSums(deviations)
    candidates[candidate_dev == min(candidate_dev)] <- TRUE
  }
  cand_idx <- which(candidates)
  if (length(cand_idx) > 1) {
    best_candidate <- cand_idx[ceiling(length(cand_idx) / 2)]
    message(paste0(
      "There are several suitable candidates for the ",
      "region of acceptance for ", tss_name, ":"
    ))
    message("Lower bound,   upper bound,   lower_rej_prob,   upper_rej_prob")
    message(paste0(
      intervals[cand_idx[1], 1], ",  ", intervals[cand_idx[1], 2], ",  ",
      intervals[cand_idx[1], 3], ",  ", intervals[cand_idx[1], 4],
      "<- leftmost"
    ))
    message(paste0(
      intervals[best_candidate, 1], ",  ", intervals[best_candidate, 2], ",  ",
      intervals[best_candidate, 3], ",  ", intervals[best_candidate, 4],
      "<- chosen"
    ))
    message(paste0(
      intervals[cand_idx[length(cand_idx)], 1], ",  ",
      intervals[cand_idx[length(cand_idx)], 2], ",  ",
      intervals[cand_idx[length(cand_idx)], 3], ",  ",
      intervals[cand_idx[length(cand_idx)], 4], "<- rightmost"
    ))
  } else {
    best_candidate <- cand_idx
    message(paste0(
      "There is a single suitable candidate for the ",
      "region of acceptance for ", tss_name, ":"
    ))
    message("Lower bound,   upper bound,   lower_rej_prob,   upper_rej_prob")
    message(paste(intervals[best_candidate, 1], intervals[best_candidate, 2],
      intervals[best_candidate, 3], intervals[best_candidate, 4],
      sep = ", "
    ))
  }
  message(paste0(
    "The chosen region of acceptance corresponds to the quantiles ",
    quantile_pairs[best_candidate, 1], " and ",
    quantile_pairs[best_candidate, 2], ".\n"
  ))
  intervals[best_candidate, ]
}

# Evaluate multiple candidates for the region of acceptance of an unbiased test.
.getAccRegion_unbiased <- function(tss,
                                   all_unique_null_vals, all_unique_null_probs,
                                   n, N_alt, N_intervals,
                                   correction, sig_lvl) {
  unbiased_threshold <- 1 - sig_lvl
  message(paste(
    "It may take a while to compute the region of acceptance",
    "using the method 'two-tailed-unbiased'."
  ))
  message(paste("Number of observed possible intervals:", N_intervals, "\n"))
  # --------- The alternative tree models:
  alt_models <- c(
    "yule", "pda", "etm",
    lapply(
      X = c(-0.5, 0.5),
      FUN = function(X) {
        list("aldous", X)
      }
    ),
    lapply(
      X = c(0.1, 0.5, 0.9),
      FUN = function(X) {
        list("ford", X)
      }
    ),
    lapply(
      X = c(0.1, 0.2, 0.3),
      FUN = function(X) {
        list("alt-birth-death", 1, X)
      }
    ),
    lapply(
      X = c(0.01, 1.1),
      FUN = function(X) {
        list("DCO_sym", X)
      }
    ),
    lapply(
      X = c(0.6, 1.1),
      FUN = function(X) {
        list("IF_asym", X)
      }
    ),
    lapply(
      X = c(1.1),
      FUN = function(X) {
        list("IF-diff", X)
      }
    ),
    lapply(
      X = c(0.05, 0.15, 0.25, 0.35),
      FUN = function(X) {
        list("biased", X)
      }
    ),
    lapply(
      X = c(0.99, 1.01),
      FUN = function(X) {
        list("ASB", X)
      }
    ),
    lapply(
      X = c(0),
      FUN = function(X) {
        list("lin-Brown_asym", c(X, 1))
      }
    ),
    lapply(
      X = c(0, 0.5),
      FUN = function(X) {
        list("log-Brown_sym", c(X, 0.1))
      }
    ),
    lapply(
      X = c(0.25),
      FUN = function(X) {
        list("log-Brown_asym", c(X, 0.1))
      }
    ),
    lapply(
      X = c(0.01, 0.4),
      FUN = function(X) {
        list("alt-birth-death", 1, X)
      }
    )
  )
  # --------- Compute the corresponding trees and TSS values.
  # Deduplicate alt_models to avoid generating trees twice for identical models.
  alt_models_keys <- sapply(alt_models, function(m) paste(unlist(m), collapse = "|||"))
  unique_mask <- !duplicated(alt_models_keys)
  alt_vals_unique <- lapply(
    X = alt_models[unique_mask],
    FUN = function(X) {
      getTSSdata(
        tss = tss, n = n,
        Ntrees = N_alt, tm = X
      )
    }
  )
  alt_vals <- alt_vals_unique[match(alt_models_keys, alt_models_keys[unique_mask])]
  # --------- Optimize the regions of acceptance per TSS.
  acc_regions <- matrix(NA,
    nrow = length(tss), ncol = 4,
    dimnames = list(
      rownames(all_unique_null_vals[[1]]),
      c(
        "lower_limit", "upper_limit",
        "lower_rej_prob", "upper_rej_prob"
      )
    )
  )
  quantile_pairs <- matrix(
    c(
      seq(0, sig_lvl, length.out = N_intervals),
      seq(1 - sig_lvl, 1, length.out = N_intervals)
    ),
    nrow = N_intervals, ncol = 2
  )
  for (i_tss in seq_along(tss)) {
    intervals <- .evaluate_interval_candidates(
      all_unique_null_vals[[i_tss]], all_unique_null_probs[[i_tss]],
      quantile_pairs, N_intervals, correction
    )
    # Extract the i_tss row from each alt model's TSS data matrix.
    alt_vals_row <- lapply(alt_vals, function(m) m[i_tss, ])
    accept_p <- .compute_accept_p(intervals, alt_vals_row, N_intervals, N_alt)
    acc_regions[i_tss, ] <- .select_best_interval(
      intervals, accept_p, unbiased_threshold, quantile_pairs, tss[i_tss]
    )
  }
  return(acc_regions)
}
