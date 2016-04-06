library(penalized)
library(glmnet)
library(spams)
library(lars)
data(prostate,package="ElemStatLearn")
pros <- subset(prostate,select=-train,train==TRUE)
ycol <- which(names(pros)=="lpsa")
X.unscaled <- as.matrix(pros[-ycol])
y.unscaled <- pros[[ycol]]
M <- matrix(
  colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
X.centered <- X.unscaled - M
sd.vec <- apply(X.unscaled, 2, sd)
S <- diag(1/sd.vec)
X.scaled <- X.centered %*% S
dimnames(X.scaled) <- dimnames(X.unscaled)
m <- mean(y.unscaled)
sigma <- sd(y.unscaled)
y.scaled <- (y.unscaled - m)/sigma

## X and y will be used in the various solvers.
X <- X.scaled
y <- y.scaled

fit.lars <- lars(X, y, type="lasso", normalize=FALSE)
beta.scaled <- coef(fit.lars)
pred.mat <- predict(fit.lars, X.scaled)$fit
pred.manual <- with(fit.lars, {
  ## Even when normalize=FALSE, the mean of the variables is
  ## subtracted away (but the variance stays the same).
  scale(X.scaled, meanx, FALSE) %*% t(beta.scaled) + mu
})
rbind(pred.mat[1,], pred.manual[1,])
stopifnot(pred.mat == pred.manual)
## If the prediction function for an input vector x is f(x) = b +
## beta'x, then the lars intercept b = meanx'beta + mu.
intercept.scaled <- with(fit.lars, -beta.scaled %*% meanx + mu)

lars.path.list <- list()
for(step.i in 1:nrow(beta.scaled)){
  coef.vec <- fit.lars$beta[step.i,]
  ## The loss in lars is half of the total squared error (0.5 * RSS)
  ## but the cost function in glmnet is half of the mean squared error
  ## (0.5 * RSS / N), so we divide the lars lambda value by N so that
  ## it is comparable to the glmnet lambda values.
  lambda <- fit.lars$lambda[step.i] / nrow(X.scaled)
  if(is.na(lambda))lambda <- 0
  lars.path.list[[paste(step.i)]] <- data.table(
    coef=c(intercept.scaled[[step.i]], coef.vec),
    arclength=sum(abs(coef.vec)),
    lambda, variable=c("(Intercept)", names(coef.vec)))
}
lars.path <- do.call(rbind, lars.path.list)

fit.glmnet <- glmnet(X, y, standardize=FALSE)
glmnet.path.list <- list()
for(lambda.i in 1:ncol(fit.glmnet$beta)){
  coef.vec <- coef(fit.glmnet)[, lambda.i]
  glmnet.path.list[[paste(lambda.i)]] <- data.table(
    lambda=fit.glmnet$lambda[[lambda.i]],
    variable=names(coef.vec),
    coef=coef.vec,
    arclength=sum(abs(coef.vec[-1]))
    )
}
glmnet.path <- do.call(rbind, glmnet.path.list)

lars.at.glmnet.lambda <- predict(
  fit.lars,
  s=fit.glmnet$lambda * nrow(X.scaled),
  type="coefficients", mode="lambda")
coef.mat.list <- list(
  glmnet=t(fit.glmnet$beta),
  lars=lars.at.glmnet.lambda$coefficients)
lapply(coef.mat.list, head)
coef.mat.lambda <- fit.glmnet$lambda

## penalized pacakge.
fit.penalized.list <- with(fit.glmnet, penalized(
  y, X, lambda1=min(lambda), steps=length(lambda),
  standardize=FALSE))
fit.penalized <- penalized(
  y, X,
  lambda1=1,
  model="linear",
  standardize=FALSE)
ppred <- predict(fit.penalized, X)
residuals <- ppred[,1] - y
n <- length(residuals)
ss <- sum(residuals * residuals)
## The equation for the loglik below "is essentially the likelihood,
## but with sigma=1 plugged in" and comes from
## https://github.com/cran/penalized/blob/master/R/linear.R#L18
(loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin)))
fit.penalized@loglik
stopifnot(fit.penalized@loglik == loglik)
## without adding inside the log, still equal.
stopifnot(fit.penalized@loglik == (-n/2) * (log(2*pi/n) + 1 + log(ss)))

## some checks using un-regularized glm solver.
fit.unregularized <- glm(y ~ X)
upred <- predict(fit.unregularized)
sum(dnorm(y, ppred[,1], ppred[,2], log=TRUE)) #not the loglik!
sum(dnorm(y, ppred[,1], 1, log=TRUE)) #not the loglik!

W0 <- cbind(rep(0, ncol(X)))
y.spams <- cbind(y)
spams.path.list <- list()
for(lambda in coef.mat.lambda){
  lambda1 <- lambda * nrow(X)
  W0 <- spams.fistaFlat(y.spams, X, W0, loss="square", regul="l1", lambda1=lambda1)
  ## From the docs, this solves argmin_w 0.5 || y - Xw ||_2^2 + lambda1 || w ||_1.
  pred.vec <- X %*% W0
  dual.norm <- function(x)max(abs(x))
  primal.norm <- function(x)sum(abs(x))
  gradient <- (pred.vec - y.spams)/nrow(X)
  inside.dual.norm <- t(X) %*% gradient
  dual.norm.value <- dual.norm(inside.dual.norm)
  ## Julien's thesis does the derivation for the square loss
  ## \tilde f(z) = ||y-z||_2^2 / 2n
  ## which has the Fenchel/convex conjugate function
  ## \tilde f^*(k) = n ||k||_2^2 / 2 + k'y.
  ## so the total cost function is
  ## g(w) = ||y - Xw||_2^2 / 2n + lambda * ||w||_1
  ## Note that is lambda not lambda1!
  grad.coef <- min(1, lambda/dual.norm.value)
  dual.vec <- grad.coef * gradient
  residual.vec <- y-pred.vec
  primal.cost <- sum(residual.vec * residual.vec) / 2 / nrow(X) + lambda * primal.norm(W0)
  loss.conjugate <- nrow(X)/2 * sum(dual.vec * dual.vec) + sum(y * dual.vec)
  penalty.conjugate.arg <- - t(X) %*% dual.vec / lambda
  penalty.conjugate <- dual.norm(penalty.conjugate.arg)
  ##stopifnot(penalty.conjugate <= 1)
  duality.gap <- primal.cost + loss.conjugate 
  cat(sprintf("lambda=%f, penalty conjugate=%f, duality gap=%e\n",
              lambda, penalty.conjugate, duality.gap))
  spams.path.list[[paste(lambda)]] <- data.table(
    lambda, variable=colnames(X), coef=as.numeric(W0), arclength=sum(abs(W0)))
}
spams.path <- do.call(rbind, spams.path.list)

library(ggplot2)
addColumn <- function(dt, pkg){
  data.table(dt, pkg=factor(pkg, c("spams", "glmnet", "penalized")))
}
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pkg ~ ., scales="free")+
  geom_point(
    aes(lambda, coef, colour=variable),
    addColumn(glmnet.path, "glmnet"),
    shape=1
    )+
  geom_point(
    aes(lambda, coef, colour=variable),
    addColumn(spams.path, "spams"),
    shape=1
    )+
  geom_line(
    aes(lambda, coef, colour=variable),
    lars.path)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pkg ~ ., scales="free")+
  geom_point(
    aes(arclength, coef, colour=variable),
    addColumn(glmnet.path, "glmnet"),
    shape=1
    )+
  geom_point(
    aes(arclength, coef, colour=variable),
    addColumn(spams.path, "spams"),
    shape=1
    )+
  geom_line(
    aes(arclength, coef, colour=variable),
    lars.path)

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
  gaussian=function(r, y)r-y,
  poisson=function(r, y)exp(r)-y)

getCriteria <- function
### Compute optimality criteria for param.vec.
(param.vec,
### Numeric vector of regression model parameters. The first element
### is the intercept which is not regularized, and the other elements
### are the elastic net regularized weights.
 derivative,
### Function which takes a vector of 
 partial.weights, lambda, alpha){
  intercept <- param.vec[1]
  
  weight.vec <- param.vec[-1]
  weight.vec
