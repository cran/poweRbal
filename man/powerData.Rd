% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerData.R
\name{getTSSdata}
\alias{getTSSdata}
\alias{getTSSdata_trees}
\title{Functions for generating the TSS data under a tree model}
\usage{
getTSSdata(tss, n, Ntrees = 1L, tm)

getTSSdata_trees(tss, treeList)
}
\arguments{
\item{tss}{Vector containing the names (as character) of the tree shape
statistics that should be compared. You may either use the short names
provided in \code{tssInfo} to use the already included TSS, or use the
name of a list object containing similar information as the entries in
\code{tssInfo}. Example:\cr
Use \code{"new_tss"} as the name for the list object
\code{new_tss} containing at least the function
\code{new_tss$func = function(tree){...}},
and optionally also the information \code{new_tss$short},
\code{new_tss$simple}, \code{new_tss$name}, \code{new_tss$type},
\code{new_tss$only_binary}, and \code{new_tss$safe_n}.}

\item{n}{Integer value that specifies the desired number of leaves, i.e.,
vertices with in-degree 1 and out-degree 0.}

\item{Ntrees}{Integer value (default = 1) that specifies the desired number
of generated trees.}

\item{tm}{If the respective model is included in this package, then specify
the model and its parameters by using a character or list. Available are all
options listed under parameter \code{tm} in
the documentation of function \code{genTrees} (type \code{?genTrees}).\cr
If you want to include your own tree model, then use the
name of a list object containing the function (with the two input parameters
\code{n} and \code{Ntrees}). Example: \cr
Use \code{"new_tm"} for the list object
\code{new_tm <- list(func = function(n, Ntrees){...})}.}

\item{treeList}{List of trees of class \code{multiphylo}.}
}
\value{
\code{getTSSdata} Numeric matrix of TSS values (one row per TSS).

\code{getTSSdata_trees} Numeric matrix of TSS values
(one row per TSS).
}
\description{
\code{getTSSdata} - Compute the tree shape statistics of trees generated
under a tree model for each given TSS.

\code{getTSSdata_trees} - Compute the tree shape statistics for each given
TSS and all given trees.
}
\examples{
# Example using tree models and TSS included in this package:
getTSSdata(tss = c("Colless", "Sackin"), n = 5L, Ntrees = 3L,
           tm = list("aldous", -1))
# Example using a "new" tree model and a "new" TSS provided by the user:
my_aldous <- list(func = function(n, Ntrees){
                        trees <- lapply(1:Ntrees,
                           function(x){genAldousBetaTree(n = n, BETA =5L)})
                        attr(trees, "class") <- "multiPhylo"
                        return(trees)})
my_avd <- list(func = treebalance::avgVertDep, short = "My AVD")
getTSSdata(tss = c("Colless", "my_avd"), n = 5L, Ntrees = 3L,
           tm = "my_aldous")
# Example using TSS provided in tssInfo.
getTSSdata_trees(tss = c("Colless", "Sackin"),
           treeList = genTrees(n = 5L, Ntrees = 3L, tm = "yule"))
# Example using a "new" TSS provided by the user.
my_avd <- list(func = treebalance::avgVertDep, short = "My AVD")
getTSSdata_trees(tss = c("Colless", "my_avd"),
           treeList = genTrees(n = 5L, Ntrees = 3L,
                               tm = list("IF_sym", 2)))
}
