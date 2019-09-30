require("lme4")
require("lattice")

units1 <- as.matrix(rnorm(5, 0, 1))
coefs <- as.matrix(-2:2)
n <- 200
Z1 <- matrix(0, nrow=n, ncol=length(units1))
Z1[1:40, 1] <- Z1[41:80, 2] <- Z1[81:120, 3] <- Z1[121:160, 4] <- Z1[161:200, 5] <- 1


X <- matrix(rnorm(n * length(units1), 0, 1), n, length(units1))
err <- rnorm(n)

dat <- as.vector(X %*% coefs + Z1 %*% units1) + err

mdf <- data.frame(dat, X)

exp_unit <- as.factor(rowSums(t(apply(Z1, 1, function(row) row * 1:5))))

mylmm <- lmer(dat ~ X1 + X2 + X3 + X4 + X5 + (1|exp_unit), data=mdf)



############
units2 <- as.matrix(rnorm(3, 0, 1))
Z2 <- matrix(0, nrow=n, ncol=length(units2))
Z2[1:80, 1] <- Z2[81:160, 2] <- Z2[161:200, 3] <- 1

dat2 <- as.vector(X %*% coefs + Z1 %*% units1 + Z2 %*% units2) + err

mdf2 <- data.frame(dat2, X)

exp_unit2 <- as.factor(rowSums(t(apply(Z2, 1, function(row) row * 1:3))))

mylmm2 <- lmer(dat ~ X1 + X2 + X3 + X4 + X5 + (1|exp_unit) + (1|exp_unit2), data=mdf2)




############
units3 <- as.matrix(rnorm(3, 0, 1))
Z3 <- matrix(0, nrow=n, ncol=length(units3))
Z3[1:30, 1] <- Z3[31:111, 2] <- Z3[112:200, 3] <- 1

dat3 <- as.vector(X %*% coefs + Z1 %*% units1 + Z3 %*% units3) + err

mdf3 <- data.frame(dat3, X)

exp_unit3 <- as.factor(rowSums(t(apply(Z3, 1, function(row) row * 1:3))))

mylmm3 <- lmer(dat ~ X1 + X2 + X3 + X4 + X5 + (1|exp_unit) + (1|exp_unit3), data=mdf3)
# mylmm4 <- lmer(dat ~ X1 + X2 + X3 + X4 + X5 + (1|exp_unit) + (1|exp_unit3), data=mdf3)
