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
  title="Hi-C pixels linked to volcano plot",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-01/figure-pixels-volcano-interactive.R",
  ggplot()+
    theme_bw()+
    theme_animint(width=1200, last_in_row=TRUE)+
    scale_x_continuous("Position on chr19 (kb)")+
    scale_y_continuous("Distance from position")+
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
      data=some_pixels),
  selector.types=list(
    r1r2="multiple"),
  first=list(
    r1r2=some_pixels[region1==region2, r1r2][1:5]))
if(FALSE){
  animint2pages(viz, "2025-10-01-HiC-pixels-volcano", chromote_sleep_seconds=5)
}
