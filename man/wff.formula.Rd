% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wff.R
\name{wff.formula}
\alias{wff.formula}
\title{WGCNA based fuzzy forest algorithm}
\usage{
\method{wff}{formula}(formula, data = NULL, ...)
}
\arguments{
\item{formula}{Formula object.}

\item{data}{data used in the analysis.}

\item{...}{Additional arguments}
}
\value{
An object of type \code{\link[fuzzyforest]{fuzzy_forest}}.  This
object is a list containing useful output of fuzzy forests.
In particular it contains a data.frame with list of selected features.
It also includes the random forest fit using the selected features.
}
\description{
Implements formula interface for \code{\link[fuzzyforest]{wff}}.
}
\note{
See \code{\link[fuzzyforest]{ff}} for additional arguments.
Note that the matrix, \code{Z}, of features that do not go through
the screening step must specified separately from the formula.
\code{test_features} and \code{test_y} are not supported in formula
interface.  As in the \code{randomForest} package, for large data sets
the formula interface may be substantially slower.

This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
}
\examples{
data(ctg)
y <- ctg$NSP
X <- ctg[, 2:22]
dat <- as.data.frame(cbind(y, X))
WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
screen_params <- screen_control(drop_fraction = drop_fraction,
                                keep_fraction = .25, min_ntree = min_ntree,
                                ntree_factor = ntree_factor,
                                mtry_factor = mtry_factor)
select_params <- select_control(drop_fraction = drop_fraction,
                                number_selected = 5,
                                min_ntree = min_ntree,
                                ntree_factor = ntree_factor,
                                mtry_factor = mtry_factor)
\donttest{
wff_fit <- wff(y ~ ., data=dat,
               WGCNA_params = WGCNA_params,
               screen_params = screen_params,
               select_params = select_params,
               final_ntree = 500)

#extract variable importance rankings
vims <- wff_fit$feature_list

#plot results
modplot(wff_fit)
}
data(ctg)
y <- ctg$NSP
X <- ctg[, 2:22]
dat <- as.data.frame(cbind(y, X))
WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
screen_params <- screen_control(drop_fraction = drop_fraction,
                                keep_fraction = .25, min_ntree = min_ntree,
                                ntree_factor = ntree_factor,
                                mtry_factor = mtry_factor)
select_params <- select_control(drop_fraction = drop_fraction,
                                number_selected = 5,
                                min_ntree = min_ntree,
                                ntree_factor = ntree_factor,
                                mtry_factor = mtry_factor)
\donttest{
wff_fit <- wff(y ~ ., data=dat,
               WGCNA_params = WGCNA_params,
               screen_params = screen_params,
               select_params = select_params,
               final_ntree = 500)

#extract variable importance rankings
vims <- wff_fit$feature_list

#plot results
modplot(wff_fit)
}
}
\seealso{
\code{\link[fuzzyforest]{wff}},
         \code{\link[fuzzyforest]{print.fuzzy_forest}},
         \code{\link[fuzzyforest]{predict.fuzzy_forest}},
         \code{\link[fuzzyforest]{modplot}}
}
