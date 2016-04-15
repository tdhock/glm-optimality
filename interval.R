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

derivative.list <- list(
  gaussian=function(r, y.lo, y.hi){
    res.lo <- r - y.lo
    res.hi <- r - y.hi
    numerator <- dnorm(res.hi) - dnorm(res.lo)
    denominator <- pnorm(res.lo) - pnorm(res.hi)
    ifelse(res.lo == res.hi, res.lo, numerator/denominator)
  })
curve(derivative.list$gaussian(x, 0, 0), -10, 10, asp=1)
curve(derivative.list$gaussian(x, -3, 3), -10, 10, asp=1, add=TRUE, col="blue")
points(-10, -7, col="blue") #
curve(derivative.list$gaussian(x, -5, Inf), -10, 10, asp=1, add=TRUE, col="red")
curve(derivative.list$gaussian(x, -Inf, 5), -10, 10, asp=1, add=TRUE, col="violet")

loss.list <- list(
  gaussian=function(r, y.lo, y.hi){
    res.lo <- r - y.lo
    res.hi <- r - y.hi
    ifelse(res.lo == res.hi,
           0.5 * res.lo * res.lo,
           -log(pnorm(res.lo) - pnorm(res.hi)))
  })
curve(loss.list$gaussian(x, 0, 0), -10, 10, asp=1)
curve(loss.list$gaussian(x, -0.3, 0.3), -10, 10, asp=1, add=TRUE, col="blue")
curve(loss.list$gaussian(x, -5, Inf), -10, 10, asp=1, add=TRUE, col="red")
curve(loss.list$gaussian(x, -Inf, 5), -10, 10, asp=1, add=TRUE, col="violet")

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

