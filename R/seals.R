seal_locs <- read.csv("data/dbay_20040622_seal_locs.csv")
pred <- read.csv("data/20040622_nonpup_predict.csv")
source('R/auxiliary-functions.R')



# filter to just the points near a prediction point
rad2 <- diff(range(pred$x)/100)**2 + diff(range(pred$y)/100)**2


# close <- filter_close(grid, pred, rad2)

####################
# put down a photo grid and rotate it to align with the observations
xpt2 <- seq(min(seal_locs[,'ET_X']), max(seal_locs[,'ET_X']), by=360)
ypt2 = seq(min(seal_locs[,'ET_Y']), max(seal_locs[,'ET_Y']), by=100)
g2 <- cbind(rep(xpt2, length(ypt2)), rep(ypt2, each=length(xpt2)))

#####################3
# Rotate the photo sites to roughly match the observed data
apt <- c(min(xpt2) - 1, min(ypt2)-1)
xnew <- ynew <- vector()
for (r in 1:nrow(g2)) {
  x = g2[r,1]
  y=g2[r,2]
  xnew <- c(xnew, sqrt((x - apt[1])**2 + (y - apt[2])**2) * sin(atan((x-apt[1]) /
                  (y-apt[2])) + 5*pi/180))
  ynew <- c(ynew, sqrt((x - apt[1])**2 + (y - apt[2])**2) * cos(atan((x-apt[1]) /
                  (y-apt[2])) + 5*pi/180))
}

xnew <- xnew + apt[1] - 600
ynew <- ynew + apt[2]
photo_sites <- cbind(x=xnew, y=ynew)

# valid_photos <- filter_close(photo_sites, pred[,1:2], 4 * sqrt(rad2))

#####################################
# for each observation, add it to the nearest photo
seal_counts <- rep(0, length(xnew))
for (i in 1:nrow(seal_locs)) {
  photo_dist <- apply(cbind(xnew, ynew), 1, function(s) (s[1] - seal_locs$ET_X[i])**2 +
                        (s[2] - seal_locs$ET_Y[i])**2 )
  seal_counts[which.min(photo_dist)] <- seal_counts[which.min(photo_dist)] + 1
}

#####################################
# Filter to accept only photos over the raster region
photo_filter <- filter_close(photo_sites, pred,  4*sqrt(rad2))  | seal_counts > 0
photo_sites <- photo_sites[photo_filter,]
seal_counts <- seal_counts[photo_filter]


# plot the observed points and the photo sites
plot(seal_locs[,c('ET_X','ET_Y')], xlim=range(seal_locs$ET_X), ylim=range(seal_locs$ET_Y))
par(new=TRUE)
plot(photo_sites, col='red', xlim=range(seal_locs$ET_X), ylim=range(seal_locs$ET_Y))


# Let's set up some spatial basis functions. First, put down knot locations.
# xlim <- range(pred$x)
# ylim <- range(pred$y)

steps_coarse <- 3
steps_med <- 8
# steps_coarse <- 5
# steps_med <- 9
#steps_fine <- 12

# Establish the coarse grid of knots
# x_coarse <- min(pred$x) + diff(range(pred$x)) * seq(0.1, 0.9, length.out=steps_coarse )
# y_coarse <- min(pred$y) + diff(range(pred$y)) * seq(0.1, 0.9, length.out=steps_coarse )
# grid_coarse <- cbind(x=rep(x_coarse, steps_coarse), y=rep(y_coarse, each=steps_coarse))
# 
# # filter to just the points near a prediction point
# rad2 <- diff(range(pred[,1])/100)**2 + diff(range(pred[,2])/100)**2
# close <- filter_close(grid_coarse, pred, rad2)
# grid_coarse <- grid_coarse[close,]

grid_coarse <- data.frame(x=vector(), y=vector())
grid_coarse[1,] <- c(796000, 1196000)
grid_coarse[2,] <- c(795000, 1197000)
grid_coarse[3,] <- c(797000, 1196000)
grid_coarse[4,] <- c(796600, 1200000)
grid_coarse[5,] <- c(797000, 1190000)
grid_coarse[6,] <- c(798500, 1193000)
grid_coarse[7,] <- c(798100, 1198000)







#Add a knot to fill in a blank space:
# grid_coarse <- rbind(grid_coarse, c(795200, 1197000))


# Establish the medium grid of knots
x_med <- min(pred$x) + diff(range(pred$x)) * seq(0.15, 0.85, length.out=steps_med )
y_med <- min(pred$y) + diff(range(pred$y)) * seq(0.15, 0.85, length.out=steps_med )
# x_med <- seq(xlim[1], xlim[2], length.out=steps_med )
# y_med <- seq(ylim[1], ylim[2], length.out=steps_med)
grid_med <- cbind(x=rep(x_med, steps_med), y=rep(y_med, each=steps_med))

# filter to just the points within the valid region
close <- filter_close(grid_med, seal_locs[,c('ET_X', 'ET_Y')], 4 * rad2)
grid_med <- grid_med[close,]



# Establish the fine grid of knots
# x_fine <- seq(xlim[1], xlim[2], length.out=steps_fine )
# y_fine <- seq(ylim[1], ylim[2], length.out=steps_fine)
# grid_fine <- cbind(x=rep(x_fine, steps_fine), y=rep(y_fine, each=steps_fine))
# 
# # filter to just the points within the valid region
# close <- filter_close(grid_fine, pred, rad2)
# grid_fine <- grid_fine[close,]



# Now apply the basis functions to the photo sites

#initialize the ranges
r_coarse <- 2000 
r_med <- 800
# r_fine <- 400


# Set up spatial basis function model matrices:
coarse_mm <- t(apply_basis(photo_sites, grid_coarse, r_coarse, gaussian))
med_mm <- t(apply_basis(photo_sites, grid_med, r_med, gaussian))
# fine_mm <- t(apply_basis(photo_sites, grid_fine, r_fine, gaussian))

# Filter out the columns with zero influence
coarse_filter <- colSums(coarse_mm) > 0.001
med_filter <- colSums(med_mm) > 0.001
# fine_filter <- colSums(fine_mm) > 0.001

coarse_mm <- coarse_mm[, coarse_filter]
med_mm <- med_mm[, med_filter]
# fine_mm <- fine_mm[, fine_filter]

# X <- cbind(coarse_mm, med_mm, fine_mm)
basis <- as.data.frame(coarse_mm)
basis <- cbind(basis, med_mm)
names(basis) <- paste0("X", 1:ncol(basis))
basis$cnt <- seal_counts

# Estimate a model:
wt <- rep(0.12*0.08, nrow(basis))
seal_density_model <- glm(cnt ~ ., family='poisson', data=basis)

over_dispersion <- mean(residuals(seal_density_model, type='pearson')[predict(seal_density_model) > quantile(predict(seal_density_model), 0.75)]**2)
# [1] 5.676487

smm <- model.matrix(seal_density_model)
covar <- matrix(0, ncol(smm), ncol(smm))
for (i in 1:nrow(smm)) {
  covar <- covar + smm[i,] %*% t(smm[i,]) * fitted.values(seal_density_model)[i] 
}
covar <- covar * 0.08 * 0.12 #scale for the size of the cells
A <- solve(covar)


nll <- function(bws, y, wt, apply_basis, photo_sites, grid_coarse, grid_med, gaussian) {
  ifelse(bws[1] < bws[2], Inf,
      {
        # Set up spatial basis function model matrices:
        coarse_mm <- t(apply_basis(photo_sites, grid_coarse, bws[1], gaussian))
        med_mm <- t(apply_basis(photo_sites, grid_med, bws[2], gaussian))
        
        # Filter out the columns with zero influence
        coarse_filter <- colSums(coarse_mm) > 0.001
        med_filter <- colSums(med_mm) > 0.001
        
        coarse_mm <- coarse_mm[, coarse_filter]
        med_mm <- med_mm[, med_filter]
        
        # X <- cbind(coarse_mm, med_mm, fine_mm)
        basis <- as.data.frame(coarse_mm)
        basis <- cbind(basis, med_mm)
        names(basis) <- paste0("X", 1:ncol(basis))
        basis$cnt=y
    
        glm(cnt ~ ., family='poisson', data=basis, weights = wt)$deviance
        
        # seal_density_model$deviance
      }
  )
}

optim(c(1200, 800), fn=nll, y=seal_counts, wt=rep(0.12*0.08, nrow(basis)),
      apply_basis=apply_basis, grid_coarse=grid_coarse, grid_med=grid_med, 
      photo_sites=photo_sites, gaussian=gaussian)

#########################
# Prediction
# Set up spatial basis function prediction matrices:
coarse_pm <- t(apply_basis(pred[,1:2], grid_coarse, r_coarse, gaussian))
med_pm <- t(apply_basis(pred[,1:2], grid_med, r_med, gaussian))
# fine_pm <- t(apply_basis(pred[,1:2], grid_fine, r_fine, gaussian))

coarse_pm <- coarse_pm[,coarse_filter]
med_pm <- med_pm[,med_filter]
# fine_pm <- fine_pm[,fine_filter]

X_pred <- data.frame(coarse_pm)
X_pred <- cbind(X_pred, med_pm)
names(X_pred) <- paste0("X", 1:ncol(X_pred))

link_pred <- predict(seal_density_model, X_pred, weights=0.053*0.053)
sum(0.053*0.053 / 0.12/0.08*exp(link_pred))

# compute the prediction standard error:
cc <- colSums(apply(cbind(1, X_pred), 2, function(c) c * exp(link_pred))) * (nrow(pred) * 0.053 * 0.053 - nrow(basis * 0.080 * 0.120)) / nrow(pred)
pred_var <- (t(as.matrix(cc)) %*% A %*% as.matrix(cc))
pois_var <- (sum(0.053*0.053 / 0.12/0.08*exp(link_pred)) - sum(basis$cnt))

# Compute the total standard error of our estimate:
total_std_err <- over_dispersion * sqrt(pred_var + pois_var)
