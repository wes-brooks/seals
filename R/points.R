require("deldir")

# sample the points
n <- rlambda(30)
x <- runif(n)
y <- runif(n)

# create a tesselation
dd_tesselation <- deldir(x, y, rw=c(0,1,0,1))
dd <- dd_tesselation$delsgs

# generate the CAR precision matrix
precision <- matrix(0, n, n)
for (r in 1:nrow(dd$delsgs)) {
  precision[dd$ind1[r], dd$ind2[r]] <- precision[dd$ind2[r], dd$ind1[r]] <- 1
}

mydata <- read.csv("~/Downloads/1470736/S2_Dataset/dbay_20040622_seal_locs.csv")
pred <- read.csv("~/Downloads/1470736/S1_Dataset/20040622_nonpup_predict.csv")
