#' Tree shape statistics
#'
#' \code{tssInfo} - List that provides information on available tree shape
#' statistics (TSS) from the package 'treebalance'.
#' Most of them are either balance or imbalance indices. The indices are grouped
#' by their families and otherwise sorted alphabetically by their full names.\cr
#' The following information is provided: \cr
#' - short: Abbreviation of the name (plain characters). \cr
#' - simple: Simplified full name (plain characters). \cr
#' - name: Full name (partly expressions as some names use special symbols).\cr
#' - func: Function of the TSS. \cr
#' - type: Either "tss", "bali", or "imbali" expressing what type of tree shape
#'         statistic it is. \cr
#' - only_binary: TRUE if TSS is suitable only for binary trees, FALSE if also
#'                applicable to arbitrary rooted trees. \cr
#' - safe_n : Integer vector with two entries specifying the range of leaf
#'            numbers \code{n} for which the TSS can be (safely) used, without
#'            warnings for too few leaves or values reaching Inf for too many
#'            leaves.\cr
#'            c(4,800), for example means that this TSS should only be applied
#'            on trees with 4 to 800 leaves. 'Inf' as the second entry means
#'            that there is no specific upper limit, but that the size of the
#'            tree itself and the computation time are the limiting factors.\cr
#' - col: Color for the TSS (related TSS have similar colors).
#'
#' @references
#'  - M. Fischer, L.Herbst, S. J. Kersting, L. Kühn, and K. Wicke,
#'  Tree Balance Indices - A Comprehensive Survey. Springer, 2023.
#'  ISBN: 978-3-031-39799-8
#'
#' @import treebalance
#' @export
#' @rdname tssInfo
#'
#' @examples
#' tssInfo$ALD$name
#' tssInfo$ALD$func(genYuleTree(6))
tssInfo <- list(
  #-----------------------------------------------------------------------------
  ALD = list(short="ALD", simple="Average leaf depth",
             name="Average leaf depth",
             func=treebalance::avgLeafDepI,
             type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "brown3"),
  #-----------------------------------------------------------------------------
  AVD = list(short="AVD", simple="Average vertex depth",
             name="Average vertex depth",
             func=treebalance::avgVertDep,
             type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "brown1"),
  #-----------------------------------------------------------------------------
  B1I = list(short="B1I", simple="B1 index",
             name=expression(paste(italic("B")[1], " index")),
             func=treebalance::B1I,
             type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "deeppink4"),
  #-----------------------------------------------------------------------------
  B2I = list(short="B2I", simple="B2 index",
             name=expression(paste(italic("B")[2], " index")),
             func=treebalance::B2I,
             type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "deeppink2"),
  #-----------------------------------------------------------------------------
  CherI = list(short="CherI", simple="Cherry index", name="Cherry index",
               func=treebalance::cherryI,
               type = "tss", only_binary = FALSE, safe_n = c(2,Inf),
               col = "dimgray"),
  #-----------------------------------------------------------------------------
  lnCPr = list(short="lnCPr", simple="ln(Colijn-Plazotta rank)",
               name="ln(Colijn-Plazotta rank)",
               func=function(tree){as.numeric(log(treebalance::colPlaLab(tree)))},
               type = "imbali", only_binary = TRUE, safe_n = c(2,13),
               col = "lemonchiffon1"),
  #-----------------------------------------------------------------------------
  Colless = list(short="Colless", simple="Colless index",
                 name="Colless index",
                 func=treebalance::collessI,
                 type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                 col = "darkblue"),
  #-----------------------------------------------------------------------------
  CollLike = list(short="CollLike",
                  simple="Colless-like index (f(n)=e^n, D=variance)",
                  name=expression(paste("Colless-like index (", italic("f"),"(",
                             italic("n"),")=",italic("e")^italic("n"),
                             ", ",italic("D"),"=variance)")),
                  func=function(tree){treebalance::collesslikeI(tree,
                                                        f.size = "exp",
                                                        dissim = "var")},
                  type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                  col = "deepskyblue4"),
  #-----------------------------------------------------------------------------
  corrColl = list(short="corrColl", simple="Corrected Colless index",
                  name="Corrected Colless index",
                  func=function(tree){treebalance::collessI(tree,
                                                        method = "corrected")},
                  type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                  col = "royalblue1"),
  #-----------------------------------------------------------------------------
  quadColl = list(short="quadColl", simple="Quadratic Colless index",
                  name="Quadratic Colless index",
                  func=function(tree){treebalance::collessI(tree,
                                                         method = "quadratic")},
                  type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                  col = "royalblue3"),
  #-----------------------------------------------------------------------------
  ewColl = list(short="ewColl", simple="Equal weights Colless index",
                name="Equal weights Colless index",
                func=treebalance::ewCollessI,
                type = "imbali", only_binary = TRUE, safe_n = c(3,Inf),
                col = "royalblue4"),
  #-----------------------------------------------------------------------------
  Furnas = list(short="Furnas", simple="Furnas rank",
                name="Furnas rank",
                func=function(tree){as.numeric(treebalance::furnasI(tree))},
                type = "bali", only_binary = TRUE, safe_n = c(2,792),
                col = "lavenderblush3"),
  #-----------------------------------------------------------------------------
  maxDep = list(short="maxDep", simple="Maximum depth",
                name="Maximum depth",
                func=treebalance::maxDepth,
                type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
                col = "turquoise1"),
  #-----------------------------------------------------------------------------
  maxWid = list(short="maxWid", simple="Maximum width",
                name="Maximum width",
                func=treebalance::maxWidth,
                type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
                col = "aquamarine1"),
  #-----------------------------------------------------------------------------
  mWomD = list(short="mWomD", simple="Maximum width over maximum depth",
               name="Maximum width over maximum depth",
               func=treebalance::mWovermD,
               type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
               col = "aquamarine3"),
  #-----------------------------------------------------------------------------
  modMaxDelW = list(short="modMaxDelW",
                    simple="Modified maximum difference in widths",
                    name="Modified maximum difference in widths",
                    func=treebalance::maxDelW,
                    type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
                    col = "aquamarine4"),
  #-----------------------------------------------------------------------------
  MeanIp = list(short="MeanIp", simple="Mean I' index",
                name=expression(paste("Mean ",italic("I"),
                                      "' index")),
                func=function(tree){treebalance::IbasedI(tree, method = "mean",
                                    correction = "prime", logs = F)},
                type = "imbali", only_binary = FALSE, safe_n = c(4,Inf),
                col = "olivedrab4"),
  #-----------------------------------------------------------------------------
  RogersJ = list(short="RogersJ", simple="Rogers J",
                 name=expression(paste("Rogers ",italic("J"))),
                 func=treebalance::rogersI,
                 type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                 col = "tan1"),
  #-----------------------------------------------------------------------------
  RQI = list(short="RQI", simple="Rooted quartet index",
             name="Rooted quartet index",
             func=treebalance::rQuartetI,
             type = "bali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "lightpink3"),
  #-----------------------------------------------------------------------------
  Sackin = list(short="Sackin", simple="Sackin index",
                name="Sackin index",
                func=treebalance::sackinI,
                type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
                col = "darkred"),
  #-----------------------------------------------------------------------------
  sShape = list(short="sShape", simple="s-shape statistic",
                name=expression(paste(hat(italic("s")),
                                      "-shape statistic")),
                func=treebalance::sShapeI,
                type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
                col = "darkgoldenrod3"),
  #-----------------------------------------------------------------------------
  SNI = list(short="SNI", simple="Symmetry nodes index",
             name="Symmetry nodes index",
             func=treebalance::symNodesI,
             type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
             col = "tan3"),
  #-----------------------------------------------------------------------------
  stairs1 = list(short="stairs1", simple="stairs1",  name="stairs1",
                 func=treebalance::stairs1,
                 type = "imbali", only_binary = TRUE, safe_n = c(2,Inf),
                 col = "lightblue4"),
  #-----------------------------------------------------------------------------
  stairs2 = list(short="stairs2", simple="stairs2", name="stairs2",
                 func=treebalance::stairs2,
                 type = "bali", only_binary = TRUE, safe_n = c(2,Inf),
                 col = "lightblue1"),
  #-----------------------------------------------------------------------------
  TIP = list(short="TIP", simple="Total internal path length",
             name="Total internal path length",
             func=treebalance::totIntPathLen,
             type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "firebrick3"),
  #-----------------------------------------------------------------------------
  TPL = list(short="TPL", simple="Total path length",
             name="Total path length",
             func=treebalance::totPathLen,
             type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "firebrick1"),
  #-----------------------------------------------------------------------------
  TotCoph = list(short="TotCoph", simple="Total cophenetic index",
                 name="Total cophenetic index",
                 func=treebalance::totCophI,
                 type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
                 col = "darkgoldenrod1"),
  #-----------------------------------------------------------------------------
  TotIp = list(short="TotIp", simple="Total I' index",
               name=expression(paste("Total ",italic("I"),
                                     "' index")),
               func=function(tree){treebalance::IbasedI(tree, method = "total",
                                   correction = "prime", logs = F)},
               type = "imbali", only_binary = FALSE, safe_n = c(4,Inf),
               col = "olivedrab2"),
  #-----------------------------------------------------------------------------
  VLD = list(short="VLD", simple="Variance of leaf depth",
             name="Variance of leaf depth",
             func=treebalance::varLeafDepI,
             type = "imbali", only_binary = FALSE, safe_n = c(2,Inf),
             col = "tomato1")
  #-----------------------------------------------------------------------------
)
