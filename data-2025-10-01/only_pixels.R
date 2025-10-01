library(animint2)
library(data.table)
add_ij <- function(DT)setkey(DT[, ij := paste0(i,",",j)], i, j, ij)
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")
pixel_dt[1]
pixel_region_min <- pixel_dt[, min(region1,region2)]
add_ij(pixel_dt[, let(
  i = region2-pixel_region_min+1L,
  j = region1-pixel_region_min+1L
)])
two_pixels <- pixel_dt[1:2]
two_pixels

some_pixels <- pixel_dt[region2<region1+50][region2<48200][
, neg.log10.p := -log10(p.value)][
, sign.neg.log10.p := sign(logFC)*neg.log10.p]

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
  region1=rep(region1, 4),
  region2=rep(region2, 4),
  r1r2=rep(r1r2, 4),
  sign.neg.log10.p=rep(sign.neg.log10.p, 4),
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
  scale_fill_gradient2()

r1r2_xy_mat <- rbind(
  c(0.5, -0.5),
  c(0.5, 0.5))
corner_dt[, c("x","y") := as.data.table(cbind(corner1, corner2) %*% r1r2_xy_mat)]

ggplot()+
  coord_equal()+
  geom_polygon(aes(
    x, y, fill=sign.neg.log10.p, group=r1r2),
    data=corner_dt,
    color=NA)+
  scale_fill_gradient2()


ggplot()+
  geom_point(aes(
    logFC, neg.log10.p),
    data=some_pixels)

viz <- animint(
  ggplot()+
    theme_bw()+
    theme_animint(width=1200, last_in_row=TRUE)+
    geom_polygon(aes(
      x, y, fill=sign.neg.log10.p, group=r1r2),
      data=corner_dt,
      clickSelects="r1r2",
      alpha=1,
      alpha_off=1,
      color="black",
      color_off="transparent")+
    scale_fill_gradient2(),
  ggplot()+
    theme_bw()+
    theme_animint(width=1200)+
    geom_point(aes(
      logFC, neg.log10.p),
      clickSelects="r1r2",
      color="black",
      color_off="grey",
      fill="white",
      data=some_pixels)+
    geom_point(aes(
      logFC, neg.log10.p),
      showSelected="r1r2",
      data=some_pixels))

pixels_for_two_clusters <- pixel_dt[clust %in% unique(clust)[1:2] ]
bands_for_two_clusters <- pixels_for_two_clusters[, {
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
  out <- as.data.table(band_list[[1]][xy])
  for(XY in list(xy, paste0("d",xy))){
    circ_diff_dt(out, XY)
  }
  out[
  , both_zero := ddx==0 & ddy==0][
  , keep := !c(both_zero[.N], both_zero[-.N])]
}, by=Cluster][]



