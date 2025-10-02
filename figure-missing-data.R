library(animint2)
library(data.table)
pixel_dt <- fread("data-2025-09-26/hicream_chr19_50000.tsv")[
, Cluster := factor(clust)][]

grid_dt <- do.call(CJ, lapply(pixel_dt[, .(region1, region2)], unique))[region1<=region2]
setkey(pixel_dt, region1, region2)
setkey(grid_dt, region1, region2)
join_dt <- pixel_dt[grid_dt]
NA_dt <- join_dt[is.na(clust)]

gg <- ggplot()+
  ggtitle(sprintf(
    "%d missing pixels out of %d possible in grid",
    nrow(NA_dt), nrow(grid_dt)))+
  geom_point(aes(
    region1, region2),
    shape=1,
    data=NA_dt)+
  coord_equal()+
  geom_abline(slope=1,intercept=0)
png("figure-missing-data.png", width=5, height=7, units="in", res=200)
print(gg)
dev.off()

