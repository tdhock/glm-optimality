works_with_R("3.2.3",
             data.table="1.9.7",
             ggplot2="1.0.1",
             directlabels="2015.6.17",
             microbenchmark="1.4.2.1",
             glmnet="1.9.5",
             spams="2.5")

## The spams package can be downloaded from
## http://spams-devel.gforge.inria.fr/hitcounter2.php?file=33815/spams-R-v2.5-svn2014-07-04.tar.gz
data(spam, package="ElemStatLearn")
is.y <- names(spam) == "spam"
X.unscaled <- as.matrix(spam[,!is.y])
y.unscaled <- spam[,is.y]

M <- matrix(
  colMeans(X.unscaled), nrow(X.unscaled), ncol(X.unscaled), byrow=TRUE)
X.centered <- X.unscaled - M
sd.vec <- apply(X.unscaled, 2, sd)
S <- diag(1/sd.vec)
X.scaled <- X.centered %*% S
dimnames(X.scaled) <- dimnames(X.unscaled)

## How to get a model equivalent to this one, using spams?
fit.glmnet.cv <- cv.glmnet(X.unscaled, y.unscaled, family="binomial")

## X and y will be used in the various solvers.
X <- X.scaled
y <- y.unscaled

fit.glmnet <- glmnet(X, y, standardize=FALSE, family="binomial", intercept=FALSE)
## glmnet : -loglik/nobs + lambda*penalty.
y.spams <- cbind(ifelse(y=="email", -1, 1))
spams.path.list <- list()
W0 <- cbind(rep(0, ncol(X)))
coef.dt.list <- list()
for(lambda.i in seq_along(fit.glmnet$lambda)){
  lambda <- fit.glmnet$lambda[[lambda.i]]
  cat(sprintf("%4d / %4d lambda=%f\n", lambda.i, length(fit.glmnet$lambda), lambda))
  ## spams :  argmin (1/m)sum_{j=1}^m log(1+e^(-y_j x^j' w)) + lambda1 psi(w),
  W0 <- spams.fistaFlat(y.spams, X, W0, loss="logistic", regul="l1", lambda1=lambda)
  coef.mat <- rbind(glmnet=fit.glmnet$beta[, lambda.i], spams=as.vector(W0))
  ## TODO: compute optimality criteria.
  for(pkg in rownames(coef.mat)){
    coef.vec <- coef.mat[pkg, ]
    coef.dt.list[[paste(pkg, lambda)]] <- data.table(
      pkg, lambda, variable=colnames(X),
      coef=as.numeric(coef.vec),
      arclength=sum(abs(coef.vec)))
  }
}
coef.dt <- do.call(rbind, coef.dt.list)

scatter.dt <- dcast(coef.dt, lambda + variable ~ pkg)
scatter.dt[, spams.0 := spams == 0]
scatter.dt[, glmnet.0 := glmnet == 0]
scatter.dt[glmnet.0 != spams.0,]
fig.scatter <- ggplot()+
  geom_abline(aes(slope=1, intercept=0),
              color="grey")+
  coord_equal()+
  geom_point(aes(glmnet, spams),
             shape=1,
             data=scatter.dt)

fig.path <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pkg ~ .)+
  geom_point(aes(arclength, coef, color=variable),
             data=coef.dt)
pdf("figure-logistic-criteria.pdf")
print(fig.scatter)
dev.off()
