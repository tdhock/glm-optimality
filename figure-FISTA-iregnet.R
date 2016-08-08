devtools::install_packages("anujkhare/iregnet")
devtools::install_packages("tdhock/PeakSegJoint")

library(iregnet)
data(H3K4me3.PGP.immune.4608, package="PeakSegJoint")
features.list <- list()
targets.list <- list()
for(problem.name in names(H3K4me3.PGP.immune.4608)){
  problem <- H3K4me3.PGP.immune.4608[[problem.name]]
  features.list[[problem.name]] <- colSums(problem$features)
  targets.list[[problem.name]] <- problem$target
}
X.nan <- do.call(rbind, features.list)
feature.is.finite <- apply(is.finite(X.nan), 2, all)
penalty.learning <- list(
  X.mat=X.nan[, feature.is.finite],
  y.mat=do.call(rbind, targets.list))

## TODO: add this data set to the iregnet package?
save(penalty.learning, file="penalty.learning.RData")
prompt(penalty.learning, file="penalty.learning.Rd")

## Fit iregnet model normal AFT model.
with(penalty.learning, {
  fit <- iregnet(X.mat, y.mat)
})

### TODO: compare iregnet solution with the following FISTA solver
### code which was copied and adapted from
### https://github.com/tdhock/PeakSegJoint/blob/master/R/fista.R

IntervalRegressionMatrix <- function
### Solve the squared hinge loss interval regression problem for one
### regularization parameter: w* = argmin_w L(w) + regularization *
### ||w||_1 where L(w) is the average squared hinge loss with respect
### to the targets, and ||w||_1 is the L1-norm of the weight vector
### (excluding the first element, which is the un-regularized
### intercept or bias term).
(features,
### Scaled numeric feature matrix (problems x features). The first
### column/feature should be all ones and will not be regularized.
 targets,
### Numeric target matrix (problems x 2).
 initial.param.vec,
### initial guess for weight vector (features).
 regularization,
### Degree of L1-regularization.
 threshold=1e-2,
### When the stopping criterion gets below this threshold, the
### algorithm stops and declares the solution as optimal.
 max.iterations=1e5,
### Error if the algorithm has not found an optimal solution after
### this many iterations.
 weight.vec=NULL,
### A numeric vector of weights for each training example.
 Lipschitz=NULL,
### A numeric scalar or NULL, which means to compute Lipschitz as the
### mean of the squared L2-norms of the rows of the feature matrix.
 verbose=2
### Cat messages: for restarts and at the end if >= 1, and for every
### iteration if >= 2.
 ){
  stopifnot(is.matrix(features))
  stopifnot(is.numeric(features))
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(is.matrix(targets))
  stopifnot(nrow(targets) == n.problems)
  stopifnot(ncol(targets) == 2)

  if(is.null(weight.vec)){
    weight.vec <- rep(1, n.problems)
  }
  stopifnot(is.numeric(weight.vec))
  stopifnot(length(weight.vec) == n.problems)

  if(is.null(Lipschitz)){
    Lipschitz <- mean(rowSums(features * features) * weight.vec)
  }
  stopifnot(is.numeric(Lipschitz))
  stopifnot(length(Lipschitz) == 1)

  stopifnot(is.numeric(max.iterations))
  stopifnot(length(max.iterations) == 1)

  stopifnot(is.numeric(threshold))
  stopifnot(length(threshold) == 1)

  stopifnot(is.numeric(initial.param.vec))
  stopifnot(length(initial.param.vec) == n.features)

  ## Return 0 for a negative number and the same value otherwise.
  positive.part <- function(x){
    ifelse(x<0, 0, x)
  }
  squared.hinge <- function(x){
    ifelse(x<1,(x-1)^2,0)
  }
  squared.hinge.deriv <- function(x){
    ifelse(x<1,2*(x-1),0)
  }  
  calc.loss <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge(linear.predictor-targets[,1])
    right.term <- squared.hinge(targets[,2]-linear.predictor)
    both.terms <- left.term+right.term
    weighted.loss.vec <- both.terms * weight.vec
    mean(weighted.loss.vec)
  }
  calc.grad <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge.deriv(linear.predictor-targets[,1])
    right.term <- squared.hinge.deriv(targets[,2]-linear.predictor)
    full.grad <- features * (left.term-right.term) * weight.vec
    colSums(full.grad)/nrow(full.grad)
  }    
  calc.penalty <- function(x){
    regularization * sum(abs(x[-1]))
  }
  calc.cost <- function(x){
    calc.loss(x) + calc.penalty(x)
  }
  soft.threshold <- function(x,thresh){
    ifelse(abs(x) < thresh, 0, x-thresh*sign(x))
  }
  ## do not threshold the intercept.
  prox <- function(x,thresh){
    x[-1] <- soft.threshold(x[-1],thresh)
    x
  }
  ## p_L from the fista paper.
  pL <- function(x,L){
    grad <- calc.grad(x)
    prox(x - grad/L, regularization/L)
  }
  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  iterate.count <- 1
  stopping.crit <- threshold
  last.iterate <- this.iterate <- y <- initial.param.vec
  this.t <- 1
  while({
    ##browser(expr=is.na(stopping.crit))
    ##str(stopping.crit)
    stopping.crit >= threshold
  }){
    ## here we implement the FISTA method with constant step size, as
    ## described by in the Beck and Tebolle paper.
    last.iterate <- this.iterate
    this.iterate <- pL(y, Lipschitz)
    last.t <- this.t
    this.t <- (1+sqrt(1+4*last.t^2))/2
    y <- this.iterate + (last.t - 1)/this.t*(this.iterate-last.iterate)
    ## here we calculate the subgradient optimality condition, which
    ## requires 1 more gradient evaluation per iteration.
    after.grad <- calc.grad(this.iterate)
    w.dist <- dist2subdiff.opt(this.iterate[-1],after.grad[-1])
    zero.at.optimum <- c(abs(after.grad[1]),w.dist)
    stopping.crit <- max(zero.at.optimum)

    if(verbose >= 2){
      cost <- calc.cost(this.iterate)
      cat(sprintf("%10d cost %10f crit %10.7f\n",
                  iterate.count,
                  cost,
                  stopping.crit))
    }
    iterate.count <- iterate.count + 1
    if(iterate.count > max.iterations){
      Lipschitz <- Lipschitz * 1.5
      iterate.count <- 1
      if(verbose >= 1){
        cat(max.iterations, "iterations, increasing Lipschitz.",
            "crit =", stopping.crit, "\n")
      }
    }
    if(any(!is.finite(this.iterate)) || 1e100 < stopping.crit){
      if(verbose >= 1){
        cat("restarting with bigger Lipschitz.\n")
      }
      iterate.count <- 1
      stopping.crit <- threshold
      last.iterate <- this.iterate <- y <- initial.param.vec
      this.t <- 1
      Lipschitz <- Lipschitz * 1.5
    }
  }
  if(verbose >= 1){
    cat("solution with crit =", stopping.crit, "\n")
  }
  this.iterate
### Numeric vector of scaled weights w of the affine function f_w(X) =
### X %*% w for a scaled feature matrix X with the first row entirely
### ones.
}

