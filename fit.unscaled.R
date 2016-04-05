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
l2norm.vec <- apply(X.centered, 2, function(x)sqrt(sum(x*x)))
stopifnot(sum(abs(
  l2norm.vec - fit.unscaled$normx
  )) < 1e6)
