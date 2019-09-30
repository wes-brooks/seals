s <- function(t) exp(t) / (1 + exp(t))
s2 <- function(t) log(t)
threshold <- 6980
n <- 35
B <- 1000
xx <- runif(n, 1, 4)

yy <- 10 + 100 * exp(xx) + rnorm(n, 0, 100)


plot(xx, yy)

lm1 <- lm(yy ~ xx)
sigma2 <- sum(residuals(lm1)**2) / lm1$df.residual
sum(dnorm(lm1$residuals, mean=0, sd = sqrt(sigma2), log = TRUE))


X <- as.matrix(cbind(1, xx))
a <- t(sigma2 * chol2inv(qr.R(qr(X))))
a.eig <- eigen(a)
a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)


boots <- vector()
# simulations:
for (i in 1:B) {
  b <- as.matrix(rnorm(2))
  b <- t(b) %*% a.sqrt
  
  y_star <- lm1$coefficients[1] +  b[1] + xx * (lm1$coefficients[2] + b[2]) + rnorm(n, mean=0, sd=sqrt(sigma2))
  lm_star <- lm(y_star ~ xx)
  boots <- c(boots, sum(dnorm(lm_star$residuals, mean=0, sd = sqrt(sigma2), log = TRUE)))
}