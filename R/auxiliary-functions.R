# this function filters pts, keeping those within the raster
filter_close <- function(pts, raster, radius_sq) {
  apply(pts, 1, function(p) any((p[1] - raster[,1])**2 + 
                                  (p[2] - raster[,2])**2 < radius_sq))
}

#define the bisquare basis function
bisquare <- function(dist, range) {
  ifelse(dist > range, 0, (1-(dist/range)**2)**2)
}

# gaussian kernel function
gaussian <- function(dist, range) {
  dnorm(dist, mean=0, sd=range)
}

# bisquare kernel function where we pass in a squared distance, in order to avoid computing a couple of square roots.
bisquare_sq <- function(dist_sq, range_sq) {
  ifelse(dist_sq > range_sq, 0, (1-(dist_sq/range_sq))**2)
}

# calculate the value taken by each basis function at ech point in pts
apply_basis <- function(pts, knots, range, basis_fn) {
  apply(pts, 1, function(p) basis_fn(sqrt((p[1]-knots[,1])**2 + (p[2] - knots[,2])**2),
                                     range))
}

# calculate the value taken by each basis function at ech point in pts (use the squared distance to avoid some square roots)
apply_basis_sq <- function(pts, knots, range, basis_fn) {
  r2 <- range**2
  apply(pts, 1, function(p) basis_fn((p[1]-knots[,1])**2 + (p[2] - knots[,2])**2,
                                     r2))
}
