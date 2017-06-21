#' Fits fuzzy forests using WGCNA to cluster features into
#' distinct modules.  Note that a formula interface for
#' WGCNA based fuzzy forests also exists: \code{\link[fuzzyforest]{wff.formula}}.
#'
#' @title WGCNA based fuzzy forest algorithm
#' @name wff
#' @export
#' @param X                 A data.frame. Each column corresponds to a feature
#'                          vector.  WGCNA will be used to cluster the
#'                          features in X.  As a result, the features should be
#'                          all be numeric.  Non-numeric features may be input
#'                          via Z.
#' @param y                 Response vector.  For classification, y should be a
#'                          factor.  For regression, y should be
#'                          numeric.
#' @param Z                 Additional features that are not to be screened out
#'                          at the screening step.  WGCNA is not carried out on
#'                          features in Z.
#' @param WGCNA_params      Parameters for WGCNA.
#'                          See \code{\link[WGCNA]{blockwiseModules}} and
#'                          \code{\link[fuzzyforest]{WGCNA_control}} for details.
#'                          \code{WGCNA_params} is an object of type
#'                          \code{WGCNA_control}.
#' @param screen_params     Parameters for screening step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{screen_control}} for details.
#'                          \code{screen_params} is an object of type
#'                          \code{screen_control}.
#' @param select_params     Parameters for selection step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{select_control}} for details.
#'                          \code{select_params} is an object of type
#'                          \code{select_control}.
#' @param final_ntree       Number of trees grown in the final random forest.
#'                          This random forest contains all selected features.
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum terminal nodesize. 1 if classification.
#'                          5 if regression.  If the sample size is very large,
#'                          the trees will be grown extremely deep.
#'                          This may lead to issues with memory usage and may
#'                          lead to significant increases in the time it takes
#'                          the algorithm to run.  In this case,
#'                          it may be useful to increase \code{nodesize}.
#' @param test_features     A data.frame containing features from a test set.
#'                          The data.frame should contain the features in both
#'                          X and Z.
#' @param test_y            The responses for the test set.
#' @param ...               Additional arguments currently not used.
#' @return An object of type \code{\link[fuzzyforest]{fuzzy_forest}}.  This
#' object is a list containing useful output of fuzzy forests.
#' In particular it contains a data.frame with list of selected features.
#' It also includes the random forest fit using the selected features.
#' @references
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5-32.
#'
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' @examples
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#' WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' \donttest{
#' wff_fit <- wff(X, y, WGCNA_params = WGCNA_params,
#'                 screen_params = screen_params,
#'                 select_params = select_params,
#'                 final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- wff_fit$feature_list
#'
#' #plot results
#' modplot(wff_fit)
#' }
#' @seealso \code{\link[fuzzyforest]{wff.formula}},
#'          \code{\link[fuzzyforest]{print.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{predict.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{modplot}}
#' @note
#' This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
#> NULL

#' @export
#' @rdname wff
wff.default <- function(X, y, Z=NULL, WGCNA_params=WGCNA_control(power=6),
                        screen_params=screen_control(min_ntree=500),
                        select_params=select_control(min_ntree=500),
                        final_ntree=5000, num_processors=1, nodesize,
                        test_features=NULL, test_y=NULL, ...) {
  if ( !("package:WGCNA" %in% search()) ) {
    stop("WGCNA must be loaded and attached. Type library(WGCNA) to do so.")
  }
  if(class(X) != "data.frame"){
    stop("X must be a data.frame")
  }
  if((!is.null(Z)) && (class(Z) != "data.frame")){
    stop("Z must be a data.frame")
  }
  numeric_test <- sapply(X, is.numeric)
  if (sum(numeric_test) != dim(X)[2]) {
    stop("To carry out WGCNA, all columns of X must be numeric.")
  }
  CLASSIFICATION <- is.factor(y)
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  WGCNA_control <- WGCNA_params
  screen_control <- screen_params
  select_control <-  select_params
  WGCNA_args <- list(X,WGCNA_control$power)
  WGCNA_args <- c(WGCNA_args, WGCNA_control$extra_args)
  names(WGCNA_args) <- c("datExpr", "power", names(WGCNA_control$extra_args))
  bwise <- do.call("blockwiseModules", WGCNA_args)
  module_membership <- bwise$colors
  screen_drop_fraction <- screen_control$drop_fraction
  screen_keep_fraction <- screen_control$keep_fraction
  screen_mtry_factor <- screen_control$mtry_factor
  screen_ntree_factor <- screen_control$ntree_factor
  screen_min_ntree <- screen_control$min_ntree
  out <- ff(X, y, Z, module_membership,
                    screen_control, select_control, final_ntree,
                    num_processors, nodesize=nodesize,
                    test_features=test_features, test_y=test_y)
  out$WGCNA_object <- bwise
  return(out)
}

#' @export
#' @rdname wff
wff <- function(X, ...) {
  UseMethod("wff", X)
}

#' Implements formula interface for \code{\link[fuzzyforest]{wff}}.
#'
#' @title WGCNA based fuzzy forest algorithm
#' @export
#' @param formula           Formula object.
#' @param data              data used in the analysis.
#' @param ...               Additional arguments
#' @return An object of type \code{\link[fuzzyforest]{fuzzy_forest}}.  This
#' object is a list containing useful output of fuzzy forests.
#' In particular it contains a data.frame with list of selected features.
#' It also includes the random forest fit using the selected features.
#' @examples
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#' dat <- as.data.frame(cbind(y, X))
#' WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' \donttest{
#' wff_fit <- wff(y ~ ., data=dat,
#'                WGCNA_params = WGCNA_params,
#'                screen_params = screen_params,
#'                select_params = select_params,
#'                final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- wff_fit$feature_list
#'
#' #plot results
#' modplot(wff_fit)
#' }
#' @note See \code{\link[fuzzyforest]{ff}} for additional arguments.
#' Note that the matrix, \code{Z}, of features that do not go through
#' the screening step must specified separately from the formula.
#' \code{test_features} and \code{test_y} are not supported in formula
#' interface.  As in the \code{randomForest} package, for large data sets
#' the formula interface may be substantially slower.
#'
#' This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
#' @seealso \code{\link[fuzzyforest]{wff}},
#'          \code{\link[fuzzyforest]{print.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{predict.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{modplot}}
#' @examples
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#' dat <- as.data.frame(cbind(y, X))
#' WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' \donttest{
#' wff_fit <- wff(y ~ ., data=dat,
#'                WGCNA_params = WGCNA_params,
#'                screen_params = screen_params,
#'                select_params = select_params,
#'                final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- wff_fit$feature_list
#'
#' #plot results
#' modplot(wff_fit)
#' }
wff.formula <- function(formula, data=NULL, ...){
  #code is stolen from randomForest by way of e1071
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  ## Catch test_features and test_y in arguments.
  if (any(c("test_features", "test_y") %in% names(m)))
    stop("test_features/test_y not supported through the formula interface")
  names(m)[2] <- "formula"
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  #rn <- 1:nrow(m)

  y <- model.response(m)
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  ## Drop any "negative" terms in the formula.
  ## test with:
  ## randomForest(Fertility~.-Catholic+I(Catholic<50),data=swiss,mtry=2)
  m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)),
                   data.frame(m))
  ## if (!is.null(y)) m <- m[, -1, drop=FALSE]
  for (i in seq(along=m)) {
    if (is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
  }
  ret <- wff.default(m, y, ...)
  cl <- match.call()
  cl[[1]] <- as.name("fuzzy_forest")
  ret$call <- cl
  ret$terms <- Terms
  class(ret) <- c("fuzzy_forest.formula", "fuzzy_forest")
  return(ret)
}




