## ------------------------------------------------------------------------
source("R/auxiliary-functions.R")
require(ggplot2)


## ------------------------------------------------------------------------
seal_locs <- read.csv("data/dbay_20040622_seal_locs.csv")
pred <- read.csv("data/20040622_nonpup_predict.csv")
photo_sites <- read.csv("data/photo_sites.csv")


## ------------------------------------------------------------------------
plot(pred$x, pred$y, col='grey70', bty='n', xlab="Easting (m)", ylab="Northing (m)", pch=15)
points(seal_locs[,c('ET_X', 'ET_Y')])
points(photo_sites, col='red', pch=22)


## ------------------------------------------------------------------------
# Each photo is a bin. We will put each seal into the nearest bin, and keep count of how many seals are in each.
seal_counts <- rep(0, nrow(photo_sites))
for (i in 1:nrow(seal_locs)) {
  photo_dist <- apply(photo_sites, 1, function(s) (s[1] - seal_locs$ET_X[i])**2 +
                        (s[2] - seal_locs$ET_Y[i])**2 )
  seal_counts[which.min(photo_dist)] <- seal_counts[which.min(photo_dist)] + 1
}


## ------------------------------------------------------------------------
grid_coarse <- read.csv("data/grid_coarse.csv")
grid_med <- read.csv("data/grid_med.csv")


## ------------------------------------------------------------------------
plot(pred$x, pred$y, col='grey70', bty='n', xlab="Easting (m)", ylab="Northing (m)", pch=15)
points(grid_coarse, col='green', pch=15)
points(grid_med, col='blue', pch=15)


## ------------------------------------------------------------------------
# set the bandwidths (here called 'r' for range)
r_coarse <- 2000 
r_med <- 800

# Set up spatial basis function model matrices:
coarse_mm <- t(apply_basis(photo_sites, grid_coarse, r_coarse, gaussian))
med_mm <- t(apply_basis(photo_sites, grid_med, r_med, gaussian))

# create the total model matrix
basis <- as.data.frame(coarse_mm)
basis <- cbind(basis, med_mm)
names(basis) <- paste0("X", 1:ncol(basis))


## ------------------------------------------------------------------------
basis$cnt <- seal_counts
seal_density_model <- glm(cnt ~ ., family='poisson', data=basis)


## ------------------------------------------------------------------------
# evaluate the basis functions at each location
coarse_pm <- t(apply_basis(pred[,1:2], grid_coarse, r_coarse, gaussian))
med_pm <- t(apply_basis(pred[,1:2], grid_med, r_med, gaussian))

# create the predictive model matrix
X_pred <- data.frame(coarse_pm)
X_pred <- cbind(X_pred, med_pm)
names(X_pred) <- paste0("X", 1:ncol(X_pred))

# make predictions and sum them up
link_pred <- predict(seal_density_model, X_pred)
popest <- sum(0.053 * 0.053 / (0.12 * 0.08) * exp(link_pred))


## ------------------------------------------------------------------------
ggdf <- pred
ggdf$fit <- link_pred
ggplot(ggdf) + aes(x=x,y=y, color=exp(fit)) + geom_point()


## ------------------------------------------------------------------------
# estimate the overdispersion parameter
over_dispersion <- mean(residuals(seal_density_model, type='pearson')[predict(seal_density_model) >
                                                                          quantile(predict(seal_density_model), 0.75)]**2)


## ------------------------------------------------------------------------
# Approximate the covariance matrix of the coefficients
smm <- model.matrix(seal_density_model)
precision <- matrix(0, ncol(smm), ncol(smm))
for (i in 1:nrow(smm)) {
  precision <- precision + smm[i,] %*% t(smm[i,]) * fitted.values(seal_density_model)[i] 
}
precision <- precision * 0.08 * 0.12 #scale for the size of the cells
A <- solve(precision)


## ------------------------------------------------------------------------
# compute the prediction standard error that comes from uncertainty in the coefficient estimates
cc <- colSums(apply(cbind(1, X_pred), 2, function(c) c * exp(link_pred))) * 
  (nrow(pred) * 0.053 * 0.053 - nrow(basis * 0.080 * 0.120)) / nrow(pred)
pred_var <- (t(as.matrix(cc)) %*% A %*% as.matrix(cc))


## ------------------------------------------------------------------------
# Compute the standard error coming from the Poisson assumption.
# Note that w have to adjust for the different sizes of the photos from the raster cells.
pois_var <- (sum(0.053*0.053 / 0.12/0.08*exp(link_pred)) - sum(basis$cnt))

# Compute the overall standard error of our estimate:
total_std_err <- over_dispersion * sqrt(pred_var + pois_var)

