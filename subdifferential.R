library(penalized)
library(glmnet)
data(prostate,package="ElemStatLearn")
pros <- subset(prostate,select=-train,train==TRUE)
ycol <- which(names(pros)=="lpsa")
X.unscaled <- as.matrix(pros[-ycol])
y.unscaled <- pros[[ycol]]
library(lars)
fit.unscaled <- lars(X.unscaled, y.unscaled, type="lasso", normalize=TRUE)
## the returned beta matrix is scaled to be in original X units, so to
## recover the scaled coefficients for plotting we need to use normx.
beta.unscaled <- scale(coef(fit.unscaled), FALSE, 1/fit.unscaled$normx)
pred.mat <- predict(fit.unscaled, X.unscaled)$fit
pred.manual <- with(fit.unscaled, {
  ## For prediction we do not need to use normx since the returned
  ## beta matrix is already scaled to be in units of the original X
  ## variables.
  scale(X.unscaled, meanx, FALSE) %*% t(beta) + mu
})
rbind(pred.mat[1,], pred.manual[1,])
stopifnot(pred.mat == pred.manual)

M <- matrix(
  colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
X.centered <- X.unscaled - M
l2norm.vec <- apply(X.centered, 2, function(x)sqrt(sum(x*x)))
stopifnot(sum(abs(
  l2norm.vec - fit.unscaled$normx
  )) < 1e6)
sd.vec <- apply(X.unscaled, 2, sd)
S <- diag(1/sd.vec)
X.scaled <- X.centered %*% S
dimnames(X.scaled) <- dimnames(X.unscaled)
m <- mean(y.unscaled)
sigma <- sd(y.unscaled)
y.scaled <- (y.unscaled - m)/sigma
fit.scaled <- lars(X.unscaled, y.unscaled, type="lasso", normalize=FALSE)
beta.scaled <- coef(fit.scaled)
pred.mat <- predict(fit.scaled, X.scaled)$fit
pred.manual <- with(fit.scaled, {
  ## Even when normalize=FALSE, the mean of the variables is
  ## subtracted away (but the variance stays the same).
  scale(X.scaled, meanx, FALSE) %*% t(beta.scaled) + mu
})
rbind(pred.mat[1,], pred.manual[1,])
stopifnot(pred.mat == pred.manual)
## If the prediction function for an input vector x is f(x) = b +
## beta'x, then the lars intercept b = meanx'beta + mu.
intercept.scaled <- with(fit.scaled, -beta.scaled %*% meanx + mu)

lars.path.list <- list()
for(step.i in 1:nrow(beta.scaled)){
  coef.vec <- fit.scaled$beta[step.i,]
  ## The loss in lars is half of the total squared error (0.5 * RSS)
  ## but the cost function in glmnet is half of the mean squared error
  ## (0.5 * RSS / N), so we divide the lars lambda value by N so that
  ## it is comparable to the glmnet lambda values.
  lambda <- fit.scaled$lambda[step.i] / nrow(X.scaled)
  if(is.na(lambda))lambda <- 0
  lars.path.list[[paste(step.i)]] <- data.table(
    coef=c(intercept.scaled[[step.i]], coef.vec),
    arclength=sum(abs(coef.vec)),
    lambda, variable=c("(Intercept)", names(coef.vec)))
}
lars.path <- do.call(rbind, lars.path.list)

gfit.scaled <- glmnet(X.unscaled, y.unscaled, standardize=FALSE)
glmnet.path.list <- list()
for(lambda.i in 1:ncol(gfit.scaled$beta)){
  coef.vec <- coef(gfit.scaled)[, lambda.i]
  glmnet.path.list[[paste(lambda.i)]] <- data.table(
    lambda=gfit.scaled$lambda[[lambda.i]],
    variable=names(coef.vec),
    coef=coef.vec,
    arclength=sum(abs(coef.vec[-1]))
    )
}
glmnet.path <- do.call(rbind, glmnet.path.list)

pfit.list <- with(gfit.scaled, penalized(
  y.unscaled, X.unscaled, lambda1=min(lambda), steps=length(lambda),
  standardize=FALSE))
pfit <- penalized(
  y.unscaled, X.unscaled,
  lambda1=1,
  model="linear",
  standardize=FALSE)
pfit10 <- penalized(
  y.unscaled, X.unscaled,
  lambda1=10,
  model="linear",
  standardize=FALSE)
ufit <- glm(y.unscaled ~ X.unscaled)
upred <- predict(ufit)
ppred <- predict(pfit, X.unscaled)
stopifnot(max(abs(upred - ppred[,1])) < 1e-10)
sum(dnorm(y.unscaled, ppred[,1], ppred[,2], log=TRUE)) #not the lik!
sum(dnorm(y.unscaled, ppred[,1], 1, log=TRUE)) #not the lik!
residuals <- ppred[,1] - y.unscaled
n <- length(residuals)
ss <- sum(residuals * residuals)
(loglik <- (-n/2) * (log(2*pi/n) + 1 + log(ss + .Machine$double.xmin)))
(-n/2) * (log(2*pi/n) + 1 + log(ss))

coef(pfit)
gfit.scaled$beta[, 10]

library(ggplot2)
addColumn <- function(dt, pkg){
  data.table(dt, pkg=factor(pkg, c("glmnet", "penalized")))
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
