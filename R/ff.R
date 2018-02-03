#' Fits the fuzzy forests algorithm. Note that a formula interface for
#' fuzzy forests also exists: \code{\link[fuzzyforest]{ff.formula}}.
#'
#' @title Fuzzy forests algorithm
#' @name ff
#' @param X                 A data.frame.
#'                          Each column corresponds to a feature vectors.
#' @param y                 Response vector.  For classification, y should be a
#'                          factor.  For regression, y should be
#'                          numeric.
#' @param Z                 A data.frame. Additional features that are not to be
#'                          screened out at the screening step.
#' @param module_membership A character vector giving the module membership of
#'                          each feature.
#' @param screen_params     Parameters for screening step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{screen_control}} for
#'                          details. \code{screen_params} is an object of type
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
#' In particular it contains a data.frame with a list of selected the features.
#' It also includes a random forest fit using the selected features.
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
#' @seealso \code{\link[fuzzyforest]{ff.formula}},
#'          \code{\link[fuzzyforest]{print.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{predict.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{modplot}}
#' @examples
#' #ff requires that the partition of the covariates be previously determined.
#' #ff is also handy if the user wants to test out multiple settings of WGCNA
#' #prior to running fuzzy forests.
#'
#' library(mvtnorm)
#' gen_mod <- function(n, p, corr) {
#'   sigma <- matrix(corr, nrow=p, ncol=p)
#'   diag(sigma) <- 1
#'   X <- rmvnorm(n, sigma=sigma)
#'   return(X)
#' }
#'
#' gen_X <- function(n, mod_sizes, corr){
#'   m <- length(mod_sizes)
#'   X_list <- vector("list", length = m)
#'   for(i in 1:m){
#'     X_list[[i]] <- gen_mod(n, mod_sizes[i], corr[i])
#'   }
#'   X <- do.call("cbind", X_list)
#'   return(X)
#' }
#'
#' err_sd <- .5
#' n <- 500
#' mod_sizes <- rep(25, 4)
#' corr <- rep(.8, 4)
#' X <- gen_X(n, mod_sizes, corr)
#' beta <- rep(0, 100)
#' beta[c(1:4, 76:79)] <- 5
#' y <- X%*%beta + rnorm(n, sd=err_sd)
#' X <- as.data.frame(X)
#'
#' Xtest <- gen_X(n, mod_sizes, corr)
#' ytest <- Xtest%*%beta + rnorm(n, sd=err_sd)
#' Xtest <- as.data.frame(Xtest)
#'
#' cdist <- as.dist(1 - cor(X))
#' hclust_fit <- hclust(cdist, method="ward.D")
#' groups <- cutree(hclust_fit, k=4)
#'
#' screen_c <- screen_control(keep_fraction = .25,
#'                            ntree_factor = 1,
#'                            min_ntree = 250)
#' select_c <- select_control(number_selected = 10,
#'                            ntree_factor = 1,
#'                            min_ntree = 250)
#' \donttest{
#' ff_fit <- ff(X, y, module_membership = groups,
#'              screen_params = screen_c,
#'              select_params = select_c,
#'              final_ntree = 250)
#' #extract variable importance rankings
#' vims <- ff_fit$feature_list
#'
#' #plot results
#' modplot(ff_fit)
#'
#' #obtain predicted values for a new test set
#' preds <- predict(ff_fit, new_data=Xtest)
#'
#' #estimate test set error
#' test_err <- sqrt(sum((ytest - preds)^2)/n)
#' }
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
#> NULL

#' @export
#' @rdname ff
ff.default <- function(X, y, Z=NULL, module_membership,
                       screen_params = screen_control(min_ntree=500),
                       select_params = select_control(min_ntree=500),
                       final_ntree = 5000,
                       num_processors=1, nodesize, test_features=NULL,
                       test_y=NULL, ...) {
  CLASSIFICATION <- is.factor(y)
  if ( !((mode(y)=="numeric") || is.factor(y)) ) {
    stop("y must be a numeric vector or factor")
  }
  if( (!CLASSIFICATION) && (length(unique(y)) < 5) ) {
    warning("y has 5 or fewer unique values.  In this case, we recommend
            classification instead of regression.  For classification,
            y must be a factor.")
  }
  if(!is.data.frame(X)) {
    stop("X must be a data.frame.")
  }
  if(!is.null(Z)) {
    if (!is.data.frame(Z)) {
      stop("Z must be a data.frame.")
    }
  }
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
  screen_control <- screen_params
  select_control <-  select_params
  module_list <- unique(module_membership)
  if(num_processors > 1) {
    #set up parallel backend
    cl <- parallel::makeCluster(num_processors)
    parallel::clusterCall(cl, library, package = "randomForest", character.only = TRUE)
    doParallel::registerDoParallel(cl)
    #close parallel backend on exit
    on.exit(try(parallel::stopCluster(cl), silent=TRUE))
  }
  survivors <- vector('list', length(module_list))
  drop_fraction <- screen_control$drop_fraction
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  if(ncol(X)*keep_fraction < select_control$number_selected){
    warning(c("ncol(X)*keep_fraction < number_selected", "\n",
              "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }

  for (i in 1:length(module_list)) {
    module <- X[, which(module_membership == module_list[i]), drop=FALSE]
    num_features <- ncol(module)
    #TUNING PARAMETER mtry_factor
    if(CLASSIFICATION == TRUE) {
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
      if(missing(nodesize)){
        nodesize <- 1
      }
    }
    if(CLASSIFICATION == FALSE) {
      mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
      if(missing(nodesize)){
        nodesize <- 5
      }
    }
    #TUNING PARAMETER ntree_factor
    ntree <- max(num_features*ntree_factor, min_ntree)
    #TUNING PARAMETER keep_fraction
    target <- ceiling(num_features * keep_fraction)
    while (num_features >= target){
      if(num_processors > 1) {
        rf <- foreach(ntree = rep(ntree/num_processors, num_processors),
                      .combine = combine, .packages = 'randomForest') %dorng% {
                        randomForest(module, y, ntree = ntree, mtry = mtry,
                                     importance = TRUE, scale = FALSE, nodesize=nodesize) }
      }
      if(num_processors == 1) {
        rf <- randomForest::randomForest(module, y, ntree = ntree, mtry = mtry,
                                         importance = TRUE, scale = FALSE,
                                         nodesize = nodesize)
      }
      var_importance <- randomForest::importance(rf, type=1, scale=FALSE)[, 1]
      var_importance <- var_importance[order(var_importance,
                                             decreasing=TRUE)]
      reduction <- ceiling(num_features*drop_fraction)
      if(num_features - reduction > target) {
        trimmed_varlist <- var_importance[1:(num_features - reduction)]
        features <- names(trimmed_varlist)
        module <- module[, which(names(module) %in% features)]
        num_features <- length(features)
        if(CLASSIFICATION == TRUE) {
          mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
        }
        if(CLASSIFICATION == FALSE) {
          mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
        }
        ntree <- max(num_features*ntree_factor, min_ntree)
      }
      else {
        num_features <- target - 1
        mod_varlist <- var_importance[1:target]
        features <- names(var_importance)[1:target]
        survivors[[i]] <- cbind(features, mod_varlist)
        row.names(survivors[[i]]) <- NULL
        survivors[[i]] <- as.data.frame(survivors[[i]])
        survivors[[i]][, 1] <- as.character(survivors[[i]][, 1])
        survivors[[i]][, 2] <- as.numeric(as.character(survivors[[i]][, 2]))
      }
    }
  }
  survivor_list <- survivors
  names(survivor_list) <- module_list
  survivors <- do.call('rbind', survivors)
  survivors <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2] <- as.numeric(survivors[, 2])
  names(survivors) <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors[, 1]]
  if(!is.null(Z)) {
    X_surv <- cbind(X_surv, Z, stringsAsFactors=FALSE)
  }
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  select_results <- do.call("select_RF", select_args)
  final_list <- select_results[[1]][, 1, drop=F]
  selection_list <- select_results[[2]]
  row.names(final_list) <- NULL
  colnames(final_list) <- c("feature_name")
  final_list <- as.data.frame(final_list, stringsAsFactors=FALSE)
  #VIMs from last tree in recursive feature elimination should be
  #replaced.
  final_list <- cbind(final_list,
                      matrix(rep(".", 2*dim(final_list)[1]), ncol=2),
                      stringsAsFactors=F)
  final_X <- X[, names(X) %in% final_list[, 1], drop=FALSE]
  #Some selected features may be from Z
  if(!is.null(Z)) {
    final_X <- cbind(final_X, Z[, names(Z) %in% final_list[, 1], drop=FALSE],
                     stringsAsFactors=FALSE)
  }
  current_p <- dim(final_X)[2]
  if(CLASSIFICATION == TRUE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*sqrt(current_p)),
                      current_p)
  }
  if(CLASSIFICATION == FALSE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p/3),
                      current_p)
  }
  if(!is.null(test_features)) {
    test_features <- test_features[, which(names(test_features) %in%
                                             names(final_X))]
  }
  final_rf <- randomForest::randomForest(x=final_X, y=y, mtry=final_mtry, ntree=final_ntree,
                            importance=TRUE, nodesize=nodesize,
                            xtest=test_features, ytest=test_y)
  final_importance <- randomForest::importance(final_rf, type=1, scale = F)
  final_list[, 1] <- row.names(final_importance)
  final_list[, 2] <- final_importance[, 1]
  #Now it's very important to associate the right module to the right
  #feature.  The ordering must be correct.  This is made trickier by
  #by the fact that when Z is not null, there exist elements in the
  #the VIM list that aren't in X.

  #select_X is a vector with selected features in order of X.
  select_X <- names(X)[which(names(X) %in% final_list[, 1])]
  #select_mods is a vector with associated module memberships in order of X.
  select_mods <- module_membership[which(names(X) %in% final_list[, 1])]
  #select_order is a vector with selected features given according to
  #the order returned by randomForest.
  select_order <- final_list[, 1][which(final_list[,1] %in% names(X))]
  #select_mods is a vector with module memberships reordered according
  #to the order returned by randomForest
  select_mods <- select_mods[match(select_order, select_X)]
  #Here is where the module membership is entered into the table.
  #Note that for elements of Z, the entry will be "."
  final_list[, 3][final_list[, 1] %in% names(X)] <- select_mods
  names(final_list)[2:3] <- c("variable_importance", "module_membership")
  #Reorder vims so that they are in decreasing order.
  final_list <- final_list[order(final_list[, 2], decreasing=T), ]
  module_membership <- as.data.frame(cbind(names(X), module_membership),
                                     stringsAsFactors=FALSE)
  names(module_membership) <- c("feature_name", "module")
  out <- fuzzy_forest(final_list, final_rf, module_membership,
                      survivor_list=survivor_list,
                      selection_list=selection_list)
  return(out)
}

#' Fits fuzzy forest algorithm.
#'
#' @rdname ff
#' @export
ff <- function(X, ...) {
  UseMethod("ff", X)
}

#' Implements formula interface for \code{\link[fuzzyforest]{ff}}.
#'
#' @title Fuzzy forests algorithm
#' @export
#' @param formula           Formula object.
#' @param data              data used in the analysis.
#' @param module_membership A character vector giving the module membership
#'                          of each feature.
#' @param ...               Additional arguments
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
#' @seealso \code{\link[fuzzyforest]{ff}},
#'          \code{\link[fuzzyforest]{print.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{predict.fuzzy_forest}},
#'          \code{\link[fuzzyforest]{modplot}}
#' @examples
#' #ff requires that the partition of the covariates be previously determined.
#' #ff is also handy if the user wants to test out multiple settings of WGCNA
#' #prior to running fuzzy forests.
#' library(mvtnorm)
#' gen_mod <- function(n, p, corr) {
#'   sigma <- matrix(corr, nrow=p, ncol=p)
#'   diag(sigma) <- 1
#'   X <- rmvnorm(n, sigma=sigma)
#'   return(X)
#' }
#'
#' gen_X <- function(n, mod_sizes, corr){
#'   m <- length(mod_sizes)
#'   X_list <- vector("list", length = m)
#'   for(i in 1:m){
#'     X_list[[i]] <- gen_mod(n, mod_sizes[i], corr[i])
#'   }
#'   X <- do.call("cbind", X_list)
#'   return(X)
#' }
#'
#' err_sd <- .5
#' n <- 500
#' mod_sizes <- rep(25, 4)
#' corr <- rep(.8, 4)
#' X <- gen_X(n, mod_sizes, corr)
#' beta <- rep(0, 100)
#' beta[c(1:4, 76:79)] <- 5
#' y <- X%*%beta + rnorm(n, sd=err_sd)
#' X <- as.data.frame(X)
#' dat <- as.data.frame(cbind(y, X))
#'
#' Xtest <- gen_X(n, mod_sizes, corr)
#' ytest <- Xtest%*%beta + rnorm(n, sd=err_sd)
#' Xtest <- as.data.frame(Xtest)
#'
#' cdist <- as.dist(1 - cor(X))
#' hclust_fit <- hclust(cdist, method="ward.D")
#' groups <- cutree(hclust_fit, k=4)
#'
#' screen_c <- screen_control(keep_fraction = .25,
#'                            ntree_factor = 1,
#'                            min_ntree = 250)
#' select_c <- select_control(number_selected = 10,
#'                            ntree_factor = 1,
#'                            min_ntree = 250)
#' \donttest{
#' ff_fit <- ff(y ~ ., data=dat,
#'              module_membership = groups,
#'              screen_params = screen_c,
#'              select_params = select_c,
#'              final_ntree = 250)
#' #extract variable importance rankings
#' vims <- ff_fit$feature_list
#'
#' #plot results
#' modplot(ff_fit)
#'
#' #obtain predicted values for a new test set
#' preds <- predict(ff_fit, new_data=Xtest)
#'
#' #estimate test set error
#' test_err <- sqrt(sum((ytest - preds)^2)/n)
#' }
#' @note See \code{\link[fuzzyforest]{ff}} for additional arguments.
#' Note that the matrix, \code{Z}, of features that do not go through
#' the screening step must specified separately from the formula.
#' \code{test_features} and \code{test_y} are not supported in formula
#' interface.  As in the \code{randomForest} package, for large data sets
#' the formula interface may be substantially slower.
#'
#' This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
ff.formula <- function(formula, data=NULL, module_membership, ...){
  #code is stolen from randomForest by way of e1071
  if (!inherits(formula, "formula"))
    stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  ## Catch test_features and test_y in arguments.
  if (any(c("test_features", "test_y") %in% names(m)))
    stop("xtest/ytest not supported through the formula interface")
  names(m)[2] <- "formula"
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$module_membership <- NULL
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
  ret <- ff.default(m, y, module_membership=module_membership, ...)
  cl <- match.call()
  cl[[1]] <- as.name("fuzzy_forest")
  ret$call <- cl
  ret$terms <- Terms
  class(ret) <- c("fuzzy_forest.formula", "fuzzy_forest")
  return(ret)
}
