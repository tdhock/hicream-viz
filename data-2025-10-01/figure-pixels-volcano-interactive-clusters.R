library(animint2)
library(data.table)
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")[, Cluster := factor(clust)]

(pixels_for_two_clusters <- pixel_dt[clust %in% c(166,1099,556) ])
get_polygons <- function(DT)DT[, {
  wide <- dcast(.SD, region1 ~ region2, length, fill=0)
  m <- as.matrix(wide[,-1])
  clust_id_mat <- cbind(0, rbind(0, m, 0), 0)
  exp_one <- function(x)c(min(x)-1, x, max(x)+1)
  band_list <- isoband::isobands(
    exp_one(as.integer(colnames(m))),
    exp_one(wide$region), clust_id_mat, 0.5, 1.5)
  xy <- c('x','y')
  circ_diff_vec <- function(z)diff(c(z,z[1]))
  circ_diff_dt <- function(DT, XY)DT[
  , paste0("d",XY) := lapply(.SD, circ_diff_vec), .SDcols=XY]
  out <- as.data.table(band_list[[1]][c(xy,"id")])
  for(XY in list(xy, paste0("d",xy))){
    circ_diff_dt(out, XY)
  }
  out[
  , both_zero := ddx==0 & ddy==0
  ][!c(both_zero[.N], both_zero[-.N]), .(id, region1=y, region2=x)]
}, by=Cluster][]

bands_for_two_clusters <- get_polygons(pixels_for_two_clusters)

ggplot()+
  theme_bw()+
  geom_tile(aes(
    region2, region1, fill=Cluster),
    data=pixels_for_two_clusters,
    alpha=0.5,
    color='grey50')+
  geom_polygon(aes(
    region2, region1, fill=Cluster, group=paste(id, Cluster)),
    data=bands_for_two_clusters,
    color='black',
    alpha=0.5)+
  geom_point(aes(
    region2, region1),
    data=bands_for_two_clusters)+
  geom_point(aes(
    region2, region1),
    data=pixel_dt[region1==48141 & region2==48190],
    fill="violet",
    color="black",
    size=5,
    shape=21)+
  scale_y_reverse()+
  coord_equal()

all_clusters <- get_polygons(pixel_dt)
some_pixels <- pixel_dt[region2<region1+50][region2<48200][
, neg.log10.p := -log10(p.value)][
, sign.neg.log10.p := sign(logFC)*neg.log10.p]
some_clusters <- all_clusters[unique(some_pixels$Cluster), on="Cluster"]

ggplot()+
  coord_equal()+
  scale_y_reverse()+
  geom_tile(aes(
    region2, region1, fill=p.value),
    data=some_pixels,
    color=NA)+
  scale_fill_gradient(low="white",high="black")

ggplot()+
  coord_equal()+
  scale_y_reverse()+
  geom_tile(aes(
    region2, region1, fill=neg.log10.p),
    data=some_pixels,
    color=NA)+
  scale_fill_gradient(low="white",high="black")

ggplot()+
  coord_equal()+
  scale_y_reverse()+
  geom_tile(aes(
    region2, region1, fill=sign.neg.log10.p),
    data=some_pixels,
    color=NA)+
  scale_fill_gradient2()

setkey(some_pixels[, r1r2 := paste0(region1,"-",region2)], region1, region2, r1r2)
(corner_dt <- setkey(some_pixels[, data.table(
  .SD,
  corner1=c(region1-0.5, region1-0.5, region1+0.5, region1+0.5),
  corner2=c(region2-0.5, region2+0.5, region2+0.5, region2-0.5)
)], region1, region2, r1r2))

ggplot()+
  coord_equal()+
  scale_y_reverse()+
  geom_polygon(aes(
    corner2, corner1, fill=sign.neg.log10.p, group=r1r2),
    data=corner_dt,
    color=NA)+
  geom_polygon(aes(
    region2, region1, group=paste(id,Cluster)),
    fill="black",
    alpha=0.1,
    color="black",
    data=some_clusters)+
  scale_fill_gradient2()

r1r2_xy_mat <- rbind(
  c(0.5, -0.5),
  c(0.5, 0.5))
corner_dt[, c("x","y") := as.data.table(cbind(corner1, corner2) %*% r1r2_xy_mat)]
some_clusters[, c("x","y") := as.data.table(cbind(region1, region2) %*% r1r2_xy_mat)]
ggplot()+
  coord_equal()+
  geom_polygon(aes(
    x, y, fill=sign.neg.log10.p, group=r1r2),
    data=corner_dt,
    color=NA)+
  geom_polygon(aes(
    x, y, group=paste(id,Cluster)),
    fill="black",
    alpha=0.1,
    color="black",
    data=some_clusters)+
  scale_fill_gradient2()

ggplot()+
  geom_point(aes(
    logFC, neg.log10.p),
    data=some_pixels)

myround <- function(x,value=1)round(x*value)/value
FC.rounder <- 4
some_pixels[, let(
  round_logFC = myround(logFC,FC.rounder),
  round_neg.log10.p = myround(neg.log10.p)
)][, let(
  volcano_bin = sprintf("logFC=%s -log10(p)=%s", round_logFC, round_neg.log10.p),
  relative_logFC = logFC-round_logFC,
  relative_neg.log10.p = neg.log10.p-round_neg.log10.p
)][]
dcast(some_pixels, . ~ ., list(min,max,Nvalues=function(x)length(unique(x))), value.var=patterns("round",cols=names(some_pixels)))
some_pixel_heat <- some_pixels[, .(
  pixels=.N
), by=.(round_logFC,round_neg.log10.p,volcano_bin)]
range(some_pixel_heat$pixels)
max.neg.log10.p.at.zero.logFC <- some_pixels[round_logFC==0, max(round_neg.log10.p)]
color.scale.dt <- unique(some_pixels[round_neg.log10.p>max.neg.log10.p.at.zero.logFC, .(
  sign_round_neg.log10.p=sort(sign(logFC)*round_neg.log10.p)
)])[, logFC := sign(sign_round_neg.log10.p)/FC.rounder/4][]

viz <- animint(
  title="Hi-C pixel clusters and volcano plot",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-01/figure-pixels-volcano-interactive-clusters.R",
  pixels=ggplot()+
    ggtitle("Genomic interaction plot, select cluster")+
    theme_bw()+
    theme_animint(width=1200, colspan=3, last_in_row=TRUE)+
    scale_x_continuous("Position on chr19 (kb)")+
    scale_y_continuous("Distance from position")+
    geom_polygon(aes(
      x, y, fill=sign.neg.log10.p, group=r1r2),
      data=corner_dt,
      alpha=1,
      color="transparent")+
    geom_polygon(aes(
      x, y, group=paste(id,Cluster)),
      fill="transparent",
      alpha=1,
      color_off="grey",
      clickSelects="Cluster",
      color="black",
      data=some_clusters)+
    scale_fill_gradient2(
      "-log10(P value)*sign"),
  volcanoHeat=ggplot()+
    ggtitle("Volcano, select tile")+
    theme_bw()+
    geom_point(aes(
      logFC, abs(sign_round_neg.log10.p),
      color=sign_round_neg.log10.p),
      data=color.scale.dt)+
    scale_color_gradient2(
      "-log10(P value)*sign")+
    geom_tile(aes(
      round_logFC, round_neg.log10.p,
      fill=log10(pixels)),
      color="grey",
      size=0.5,
      data=some_pixel_heat)+
    geom_point(aes(
      logFC, neg.log10.p),
      showSelected="Cluster",
      data=some_pixels,
      color="green")+
    geom_tile(aes(
      round_logFC, round_neg.log10.p,
      tooltip=paste(volcano_bin, "pixels=", pixels)),
      fill="transparent",
      clickSelects="volcano_bin",
      data=some_pixel_heat)+
    scale_y_continuous(
      "-log10(P value)")+
    scale_fill_gradient(low="white", high="black"),
  volcanoZoom=ggplot()+
    ggtitle("Volcano zoom, select cluster")+
    theme_bw()+
    geom_point(aes(
      relative_logFC, relative_neg.log10.p),
      clickSelects="Cluster",
      showSelected="volcano_bin",
      size=4,
      color="black",
      color_off="grey",
      fill="white",
      data=some_pixels),
  first=list()
)
if(FALSE){
  animint2pages(viz, "2025-10-01-HiC-pixels-volcano-clusters", chromote_sleep_seconds=5)
}
viz
