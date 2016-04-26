works_with_R("3.2.3", ggplot2="1.0.1", data.table="1.9.7")

## There is definitely a bug in the pnorm function, since pnorm(x)
## should be the same as 1-pnorm(-x), but they are not the same! e.g.
rbind(pnorm(9.3)-pnorm(8.7),
      pnorm(-8.7)-pnorm(-9.3))
rbind(pnorm(4.3)-pnorm(3.7),
      pnorm(-3.7)-pnorm(-4.3))
fun.list <- list(
  pnorm=pnorm,
  pnorm2=function(x)1-pnorm(-x))
## on the probability scale you can't see any difference.
curve(pnorm(x), -10, 10, lwd=2, n=1001)
curve(1-pnorm(-x), -10, 10, col="red", n=1001, add=TRUE)
## on the log scale there is a difference on the left.
curve(pnorm(x, log.p=TRUE), -10, 10, lwd=2, n=1001)
curve(log(1-pnorm(-x)), -10, 10, col="red", n=1001, add=TRUE)
fun.values.list <- list()
input.vec <- seq(-10, 10, by=0.1)
for(fun.name in names(fun.list)){
  fun <- fun.list[[fun.name]]
  fun.values.list[[fun.name]] <- data.table(
    input.vec,
    log.probability=log(fun(input.vec)),
    fun.name)
}
fun.values <- do.call(rbind, fun.values.list)
ggplot()+
  geom_line(aes(input.vec, log.probability, color=fun.name),
            data=fun.values)

loss.list <- list(
  gaussian.diff=function(r, y.lo, y.hi){
    res.lo <- r - y.lo
    res.hi <- r - y.hi
    is.equal <- res.lo == res.hi
    numerator <- dnorm(res.hi) - dnorm(res.lo)
    ## There is some numerical instability! Sometimes the denominator can be zero.
    denominator <- pnorm(res.lo) - pnorm(res.hi)
    loss <- ifelse(
      is.equal,
      0.5 * res.lo * res.lo,
      -log(denominator))
    derivative <- ifelse(is.equal, res.lo, numerator/denominator)
    y.diff <- (y.hi - y.lo)/2
    is.interval <- -Inf < y.lo & y.lo < y.hi & y.hi < Inf
    constant <- ifelse(
      is.interval,
      log(pnorm(y.diff)-pnorm(-y.diff)),
      0)
    cbind(loss, derivative, constant)
  },
  gaussian.mean=function(r, y.lo, y.hi){
    res.lo <- r - y.lo
    res.hi <- r - y.hi
    is.equal <- res.lo == res.hi
    numerator <- dnorm(res.hi) - dnorm(res.lo)
    denominator <- pnorm(r, y.lo) - pnorm(r, y.hi)
    loss <- ifelse(
      is.equal,
      0.5 * res.lo * res.lo,
      -log(denominator))
    derivative <- ifelse(is.equal, res.lo, numerator/denominator)
    y.diff <- (y.hi - y.lo)/2
    is.interval <- -Inf < y.lo & y.lo < y.hi & y.hi < Inf
    constant <- ifelse(
      is.interval,
      log(pnorm(y.diff)-pnorm(-y.diff)),
      0)
    cbind(loss, derivative, constant)
  },
  logistic=function(r, y.lo, y.hi){
    res.lo <- y.lo - r
    res.hi <- r - y.hi
    is.equal <- res.lo == -res.hi
    exp.lo <- exp(res.lo)
    exp.hi <- exp(res.hi)
    log.lo <- log(1+exp.lo)
    log.hi <- log(1+exp.hi)
    loss <- ifelse(
      is.equal,
      2*log.lo - res.lo,
      log.lo + log.hi)
    derivative <- ifelse(
      is.equal,
      (1-exp.lo)/(1+exp.lo),
      exp.hi/(1+exp.hi)-exp.lo/(1+exp.lo))
    is.interval <- -Inf < y.lo & y.lo < y.hi & y.hi < Inf
    y.diff <- (y.lo-y.hi)/2
    constant <- ifelse(
      is.equal,
      -2*log(2),
      ifelse(
        is.interval,
        -2*log(1+exp(y.diff)),
        0))
    cbind(loss, derivative, constant)
  })
prediction.vec <- seq(-10, 10, l=101)
y.mat <- rbind(
  c(-Inf, 3),
  c(-3, Inf),
  c(-3, 3),
  c(-0.3, 0.3),
  c(3, 3))
loss.lines.list <- list()
for(loss.name in names(loss.list)){
  loss.fun <- loss.list[[loss.name]]
  for(y.i in 1:nrow(y.mat)){
    y.row <- y.mat[y.i,]
    y.name <- paste0("(", y.row[1], ",", y.row[2], ")")
    stopifnot(y.row[1] <= y.row[2])
    loss.mat <- loss.fun(prediction.vec, y.row[1], y.row[2])
    loss.lines.list[[paste(loss.name, y.i)]] <- data.table(
      loss.name, y.i, y.name,
      prediction=prediction.vec,
      loss=loss.mat[, "loss"] + loss.mat[, "constant"],
      derivative=loss.mat[, "derivative"])
  }
}
loss.lines <- do.call(rbind, loss.lines.list)
tangent.lines <- loss.lines[prediction==2,]
linetypes <- c(
  logistic="dashed",
  gaussian.diff="solid",
  gaussian.mean="dotted")
fig.interval.loss <- ggplot()+
  ggtitle("interval regression surrogate loss functions")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ y.name, scales="free")+
  scale_linetype_manual(values=linetypes)+
  geom_line(aes(prediction, loss, color=loss.name, linetype=loss.name),
            data=loss.lines)+
  scale_x_continuous(breaks=c(-5, 0, 5))
pdf("figure-interval-loss.pdf", w=8, h=3)
print(fig.interval.loss)
dev.off()

fig.loss.derivative <- ggplot()+
  ggtitle("tangent line verifies derivative computation")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(loss.name ~ y.name, scales="free")+
  scale_linetype_manual(values=linetypes)+
  geom_abline(aes(slope=derivative, intercept=loss-derivative*prediction,
                  color=line),
              data=data.table(tangent.lines, line="tangent"))+
  geom_point(aes(prediction, loss), data=tangent.lines, shape=1)+
  geom_line(aes(prediction, loss, color=line),
            data=data.table(loss.lines, line="loss"))+
  scale_color_manual(values=c(loss="black", tangent="grey50"))+
  scale_x_continuous(breaks=c(-5, 0, 5))
pdf("figure-interval-loss-derivative.pdf", w=8, h=5)
print(fig.loss.derivative)
dev.off()
