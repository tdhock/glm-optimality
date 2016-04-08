## This script defines functions for computing an optimality criterion
## based on the subdifferential, 0 \in \partial C(\beta, w). If there
## are d input variables, w is the d dimensional weight vector, and
## \beta is the the un-regularized intercept term.

## The cost C(\beta, w) = \lambda P_\alpha(w) + L(\beta, w), where
## \lambda is the non-negative degree of regularization, P_\alpha(w) =
## \alpha ||w||_1 + 0.5 (1-\alpha)||w||_2^2 is the elastic net
## penalty, and L(\beta, w) is the average loss with respect to the
## training examples.

## There are three criteria, and we will compute the entire vector of
## d+1 criteria, then take the norm of that vector.

## 1. For the un-regularized intercept, we compute | \nabla_\beta
## L(\beta, w) |.

## 2. For the weights that are zero we compute (|\nabla_w L(\beta, w)
## + \lambda(1-\alpha)w|-\lambda\alpha)_+ where ()_+ is the positive
## part function.

## 3. For the weights that are NOT zero we compute
## |\lambda[(1-\alpha)w + \alpha\sign(w)]|.

loss.list <- list(
  gaussian=function(r, y)0.5 * (r-y)^2,
  poisson=function(r, y)exp(r) - y*r,
  logistic=function(r, y)log(1+exp(-y*r)))
derivative.list <- list(
  gaussian=function(r, y)r-y,
  poisson=function(r, y)exp(r)-y,
  logistic=function(r, y)-y/(1+exp(y*r)))
curve(derivative.list$logistic(x, 1), -2, 2)
curve(derivative.list$logistic(x, -1), -2, 2)

subdifferentialCriteria <- function
### Compute subdifferential optimality criteria for an elastic net
### regularized GLM.
(y,
### Numeric vector of output labels (n observations)
 X,
### Numeric matrix of input features (n observations x p features).
 w,
### Numeric vector of weights (p features).
 lambda,
### Numeric lambda regularization parameter (scalar).
 alpha=1,
### Numeric elastic net parameter (between 0 and 1).
 derivative.fun=derivative.list$gaussian
### observation-specific derivative function(prediction.vec,
### label.vec) from derivative.list.
 ){
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(X))
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(w))
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1) #TODO support lambda vector w matrix.
  stopifnot(nrow(X) == length(y))
  stopifnot(ncol(X) == length(w))
  stopifnot(is.numeric(alpha))
  stopifnot(length(alpha)==1)
  stopifnot(0 <= alpha && alpha <= 1)
  pred.vec <- X %*% w
  observation.gradient.vec <- derivative.fun(pred.vec, y)
  weight.gradient.vec <- t(X) %*% observation.gradient.vec / nrow(X)
  positive.part <- function(x)ifelse(x<0, 0, x)
  common.term <- gradient.vec + lambda * (1-alpha) * weight.vec
  ifelse(weight.vec==0,
         positive.part(abs(common.term) - lambda*alpha),
         abs(common.term + lambda*alpha*sign(weight.vec)))
### p-vector of subdifferential optimality criteria. Lower values mean
### closer to the global optimum of the elastic net problem.
}

