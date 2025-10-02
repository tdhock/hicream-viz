library(animint2)
library(data.table)
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")

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
half.width <- 0.4
corners <- function(corner1, corner2)data.table(corner1, corner2)
(corner_dt <- setkey(some_pixels[, data.table(.SD, rbind(
  corners(region1-half.width, region2-half.width),
  corners(region1-half.width, region2+half.width),
  corners(region1+half.width, region2+half.width),
  corners(region1+half.width, region2-half.width),
  corners(region1-half.width, region2-half.width))
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

selected.color <- "green"
viz <- animint(
  title="Hi-C pixel heat map with borders linked to zoomable volcano plot",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-02/figure-pixels-volcano-interactive-heat-borders.R",
  pixels=ggplot()+
    ggtitle("Genomic interaction plot, select pixels")+
    theme_bw()+
    theme_animint(width=1200, colspan=3, last_in_row=TRUE)+
    scale_x_continuous("Position on chr19 (kb)")+
    scale_y_continuous("Distance from position")+
    geom_polygon(aes(
      x, y, fill=logFC, color=neg.log10.p, group=r1r2),
      data=corner_dt)+
    scale_color_gradient(
      "-log10(P value)",
      low="white", high="black",
      breaks=rev(seq(0, 20, by=5)),
      guide=guide_legend(override.aes=list(fill="white")))+
    geom_polygon(aes(
      x, y, group=r1r2),
      data=corner_dt,
      clickSelects="r1r2",
      fill="transparent",
      alpha=1,
      alpha_off=1,
      color=selected.color,
      color_off="transparent")+
    scale_fill_gradient2(
      "log(fold change)",
      breaks=rev(seq(-6, 2, by=1))),
  volcanoHeat=ggplot()+
    ggtitle("Volcano, select tile")+
    theme_bw()+
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
      color=selected.color)+
    geom_tile(aes(
      round_logFC, round_neg.log10.p,
      tooltip=paste(volcano_bin, "pixels=", pixels)),
      fill="transparent",
      clickSelects="volcano_bin",
      data=some_pixel_heat)+
    scale_x_continuous(
      "log(fold change)")+
    scale_y_continuous(
      "-log10(P value)")+
    scale_fill_gradient(low="white", high="black"),
  volcanoZoom=ggplot()+
    ggtitle("Volcano zoom, select pixels")+
    scale_x_continuous(
      "Relative log(fold change) in tile")+
    scale_y_continuous(
      "Relative -log10(P value) in tile")+
    theme_bw()+
    geom_point(aes(
      relative_logFC, relative_neg.log10.p),
      clickSelects="r1r2",
      showSelected="volcano_bin",
      size=4,
      color="black",
      color_off="grey50",
      fill="green",
      fill_off="white",
      data=some_pixels)+
    geom_point(aes(
      relative_logFC, relative_neg.log10.p),
      clickSelects="r1r2",
      showSelected=c("volcano_bin","r1r2"),
      size=4,
      color="black",
      color_off="grey50",
      fill="green",
      fill_off="white",
      data=some_pixels),
  selector.types=list(
    r1r2="multiple"),
  out.dir="figure-pixels-volcano-interactive-heat-borders",
  first=list(
    r1r2=some_pixels[region1==region2, r1r2][1:5]))
if(FALSE){
  animint2pages(viz, "2025-10-02-HiC-pixels-heat-borders", chromote_sleep_seconds=5)
}
viz
