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
lars.path.list <- list()
for(step.i in 1:nrow(beta.unscaled)){
  standardized.coef <- beta.unscaled[step.i,]
  arclength <- sum(standardized.coef)
  lambda <- fit.unscaled$lambda[step.i]
  if(is.na(lambda))lambda <- 0
  lars.path.list[[paste(step.i)]] <- data.table(
    coef=fit.unscaled$beta[step.i,],
    standardized.coef, arclength, lambda, variable=names(standardized.coef))
}
lars.path <- do.call(rbind, lars.path.list)

gfit.unscaled <- glmnet(X.unscaled, y.unscaled)
glmnet.path.list <- list()
for(lambda.i in 1:nrow(gfit.unscaled$beta)){
  glmnet.path.list[[paste(lambda.i)]] <- data.table(
    lambda=gfit.unscaled$lambda[[lambda.i]],
    variable=rownames(gfit.unscaled$beta),
    coef=gfit.unscaled$beta[, lambda.i]
    )
}
glmnet.path <- do.call(rbind, glmnet.path.list)

library(ggplot2)
addColumn <- function(dt, facet){
  data.table(dt, facet=factor(facet, c("coef", "standardized")))
}
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(facet ~ ., scales="free")+
  geom_point(
    aes(lambda, coef, colour=variable),
    glmnet.path,
    shape=1
    )+
  geom_line(
    aes(lambda, standardized.coef, colour=variable),
    addColumn(lars.path, "standardized"))+
  geom_line(
    aes(lambda, coef, colour=variable),
    addColumn(lars.path, "coef"))

M <- matrix(
  colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
sd.vec <- apply(X.unscaled, 2, sd)
S <- diag(1/sd.vec)
X.scaled <- (X.unscaled - M) %*% S
m <- mean(y.unscaled)
sigma <- sd(y.unscaled)
y.scaled <- (y.unscaled - m)/sigma
fit.scaled <- lars(X.scaled, y.scaled, type="lasso", normalize=FALSE)
beta.scaled <- coef(fit.scaled)
pred.mat <- predict(fit.scaled, X.scaled)$fit
pred.manual <- X.scaled %*% t(beta.scaled) + fit.scaled$mu
rbind(pred.mat[1,], pred.manual[1,])
pred.mat == pred.manual

rbind(beta.scaled[,1], beta.unscaled[,1])
beta.scaled[,1]/beta.unscaled[,1]
stopifnot(beta.scaled == beta.unscaled)

arclength <- rowSums(abs(beta.unscaled))
path <- data.frame(melt(beta.unscaled), arclength)
names(path)[1:3] <- c("step","variable","standardized.coef")
library(ggplot2)
ggplot(path,aes(arclength,standardized.coef,colour=variable))+
  geom_line(aes(group=variable))+
  ggtitle("LASSO path for prostate cancer data calculated using the LARS")+
  xlim(0,20)

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
