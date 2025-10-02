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
  title="Hi-C pixels linked to zoomable volcano plot",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-01/figure-pixels-volcano-interactive-zoom.R",
  pixels=ggplot()+
    ggtitle("Genomic interaction plot, select pixels")+
    theme_bw()+
    theme_animint(width=1200, colspan=3, last_in_row=TRUE)+
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
      showSelected="r1r2",
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
    ggtitle("Volcano zoom, select pixels")+
    theme_bw()+
    geom_point(aes(
      relative_logFC, relative_neg.log10.p),
      clickSelects="r1r2",
      showSelected="volcano_bin",
      size=4,
      color="black",
      color_off="grey",
      fill="white",
      data=some_pixels),
  selector.types=list(
    r1r2="multiple"),
  first=list(
    r1r2=some_pixels[region1==region2, r1r2][1:5]))
if(FALSE){
  animint2pages(viz, "2025-10-01-HiC-pixels-volcano-zoom", chromote_sleep_seconds=5)
}
viz
