units <- as.factor(rep(c("A", "B", "C", "D"), each=20))
units2 <- as.factor(rep(c("A", "B", "C", "D"), 20))
n = length(units)

#tmp <- model.matrix(~units*units2, contrasts.arg=list(units=contr.sum(4), units2=contr.sum(4)))
tmp <- cbind(1, model.matrix(~units-1), model.matrix(~units2-1), model.matrix(~units-units2-1))
#bb <- sample(4:9, size=ncol(tmp), replace=TRUE)
bb <- rnorm(ncol(tmp), mean=5, sd=5)
ref <- tmp %*% as.matrix(bb)

rho <- -0.6
xx <- rnorm(n, mean=0, sd=1)
xx2 <- rnorm(n, mean=0, sd=1)
xxx <- as.matrix(cbind(xx, xx2)) %*% matrix(c(sqrt(2), rho*sqrt(2), rho*sqrt(2), 1), 2, 2)

eps <- rnorm(n, mean=10, sd=2)
beta <- 2
fixed <- xxx[,1]*beta + eps

Y <- fixed + ref


m1 <- lmer(Y ~ xxx[,2] + (1 | units))
stored <- m1@devcomp$cmp[['REML']]

m1_vcov <- vcov.merMod(m1)
m1_vcov.eig <- eigen(m1_vcov)
m1_vcov.sqrt <- m1_vcov.eig$vectors %*% diag(sqrt(m1_vcov.eig$values)) %*% solve(m1_vcov.eig$vectors)

boots <- vector()
for (i in 1:100) {
  b_star <- as.vector(m1_vcov.sqrt %*% as.matrix(rnorm(2))) + fixef(m1)
  
  ranef_var <- as.vector(attr(ranef(m1, condVar=TRUE)[[1]], 'postVar') )
  z_star <- as.vector(ranef(m1)$units) + rnorm(length(ranef_var), 0, sd=sqrt(ranef_var))
  e_star <- rnorm(n, mean=0, sd=lme4:::sigma.merMod(m1))
  y_star <- b_star[1] + b_star[2] * xxx[,2] + as.vector(model.matrix(~units-1) %*% as.matrix(z_star)) + e_star
  
  m_star <- lmer(y_star ~ xxx[,2] + (1 | units))
  

    
  boots <- c(boots, m_star@devcomp$cmp[['REML']])
}
