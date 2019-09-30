
valid_photos <- filter_close(photo_sites, pred[,1:2], 4*sqrt(rad2))
plot(seal_locs[,c('ET_X', 'ET_Y')])
> par(new=TRUE)
> plot(photo_sites[valid_photos,], col='red', xlim=range(seal_locs$ET_X),ylim=range(seal_locs$ET_Y))

photo_sites <- photo_sites[valid_photos,]
