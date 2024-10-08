% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poweRbal_data-methods.R
\name{plot.poweRbal_data}
\alias{plot.poweRbal_data}
\title{Plot method for poweRbal_data objects}
\usage{
\method{plot}{poweRbal_data}(
  x,
  tss_names = NULL,
  tss_colors = NULL,
  sig_lvl = NULL,
  legend_pos = "topright",
  alt_model_names = NULL,
  alt_model_params = NULL,
  tss_ltys = NULL,
  alt_model_family = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{poweRbal_data}, which is a list containing
one mandatory element, \code{power}, and several optional elements:\cr
\itemize{
\item \code{power}: A numeric matrix containing the power values (one row per
TSS and one column per alternative model).\cr
\item \code{accept_regions}: A numeric matrix containing information on the
region of acceptance (one row per TSS and four columns).\cr
\item \code{CIradius}: A numeric matrix containing the confidence interval
radii (one row per TSS and one column per alternative model).\cr
\item \code{actual_sample_sizes}: A numeric vector containing the actual sample
sizes under each alternative model, as some models do not always
successfully generate trees.\cr
\item \code{alt_model_params}: A numeric vector (one element per alternative
model) containing the values of a tree model parameter. This is only
suitable if the alternative models all belong to the same tree model
family and differ only in one parameter.\cr
\item Other input data from the \code{powerComp()} function, such as \code{tss},
\code{null_model}, \code{alt_models}, \code{n}, \code{distribs},
\code{N_null}, \code{N_alt}, \code{test_type}, \code{correction}, and
\code{sig_lvl}.
}}

\item{tss_names}{Vector of characters/expression of the TSS names (default
= NULL). If none are provided, \code{x$tss} is used for the names
of the TSS if existent and otherwise the row names of \code{x$power}
are used.}

\item{tss_colors}{Vector of colors for the TSS (default = NULL).}

\item{sig_lvl}{Level of significance (default=0.05, must be >0 and <1)
depicted as a dashed horizontal line. Not depicted if set to NULL.}

\item{legend_pos}{Character specifying where the legend is displayed
(default = "topright"). No legend is displayed if set to NULL.}

\item{alt_model_names}{Vector of characters/expression of the model names
(default = NULL). If none are provided, the column names of
\code{x$power} are used as names for the models.\cr
Only used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.}

\item{alt_model_params}{Numeric vector containing the parameter
values of the representatives of the tree model (default = NULL). If none
are provided, \code{x$alt_model_params} is used if existent.}

\item{tss_ltys}{Vector of line types for the TSS (default = NULL).\cr
Not used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.}

\item{alt_model_family}{Vector of characters/expressions of the name of the
tree model family and of the parameter (default = NULL), e.g.
\code{c("Aldous\'", expression(beta))}.
If none is provided, the first column name of \code{x$power} is used.\cr
Not used if \code{alt_model_params} and \code{x$alt_model_params} = NULL.}

\item{...}{Additional arguments passed to the \code{plot} function.}
}
\value{
\code{plot.poweRbal_data} No return value, as the primary purpose
of this function is the side effect (plotting).
}
\description{
This function generates a plot for an object of class \code{poweRbal_data}.
Creates a bar plot if \code{alt_model_params} and \code{x$alt_model_params}
= NULL and a line plot otherwise if this information is given.
}
\examples{
# Plotting a 'poweRbal_data' object:
pc1 <- powerComp(tss = c("Sackin", "Colless", "B1I"),
                 alt_models = list(list("aldous",-1), "pda", "etm"),
                 n = 8L, N_null = 40L, N_alt = 20L)
plot(pc1)
# Plotting a power comparison with a tree model family
pc2 <- powerComp(tss = c("Sackin", "Colless", "B1I"),
               alt_models = list(list("aldous", -1.5),
                                 list("aldous", -1),list("aldous", -0.5),
                                 list("aldous", 0),list("aldous", 0.5)),
               n=20L, N_null = 20L, N_alt = 10L, distribs = "sampled")
# Create a bar plot or ...
plot(pc2)
# ... a line plot by specifying 'alt_model_params'.
plot(pc2, alt_model_params = c(-1.5,-1,-0.5,0,0.5),
     tss_names = getTSSnames(c("Sackin", "Colless", "B1I")),
     tss_colors = getTSScolors(c("Sackin", "Colless", "B1I")),
     alt_model_family = c("Aldous\'", expression(beta)),
     ylim = c(0,1))
}
