% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmYule.R
\name{genYuleTree}
\alias{genYuleTree}
\title{Generation of rooted binary trees under the Yule model}
\usage{
genYuleTree(n)
}
\arguments{
\item{n}{Integer value that specifies the desired number of leaves, i.e.,
vertices with in-degree 1 and out-degree 0.\cr
Due to the restrictions of the \code{phylo} or \code{multiphylo} format,
the number of leaves must be at least 2 since there must be at
least one edge.}
}
\value{
\code{genYuleTree} A single tree of class \code{phylo} is returned.
}
\description{
\code{genYuleTree} - Generates a rooted binary tree in \code{phylo}
format with the given number of \code{n} leaves under the Yule model.
The Yule process is a simple birth-process in which all species have the same
speciation rate.
}
\examples{
genYuleTree(n = 5)
}
\references{
\itemize{
\item This function uses the \code{rtree} function of the \code{ape} package
(E. Paradis, K. Schliep. “ape 5.0: an environment for modern
phylogenetics and evolutionary analyses in R.”
Bioinformatics, 35, 526-528, 2019). \cr
\item G. U. Yule. A mathematical theory of evolution, based on the conclusions of
Dr. J. C. Willis, F. R. S. Philosophical Transactions of the Royal Society
of London. Series B, Containing Papers of a Biological
Character, 213(402-410):21–87, 1925. \cr
\item E. F. Harding. The probabilities of rooted tree-shapes generated by random
bifurcation. Advances in Applied Probability, 3(1):44–77, 1971.
}
}
