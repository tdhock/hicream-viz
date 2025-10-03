library(animint2)
library(data.table)
myround <- function(x,value=1)round(x*value)/value
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")[, let(
  Cluster=factor(clust),
  neg.log10.p = -log10(p.value)
)]
for(i in 1:2){
  region_i <- pixel_dt[[paste0("region",i)]]
  round_region_i <- myround(region_i, 1/100)
  set(pixel_dt, j=paste0("round_region",i), value=round_region_i)
  set(pixel_dt, j=paste0("rel_region",i), value=region_i-round_region_i)
}
pixel_dt[
, rel_regions := paste0(rel_region1,",",rel_region2)
][, .(corners=.N), by=rel_regions]
for(prefix in c("","round_")){
  r <- function(i)pixel_dt[[paste0(prefix,"region",i)]]
  value <- paste0(r(1),"-",r(2))
  set(pixel_dt, j=paste0(prefix,"r1r2"), value=value)
}
length(unique(pixel_dt$r1r2))

region_tiles <- pixel_dt[, .(
  pixels=.N,
  mean.neg.log10.p=mean(neg.log10.p),
  mean.logFC=mean(logFC),
  clusters=length(unique(Cluster))
), by=.(round_r1r2,round_region1,round_region2)]
ggplot()+
  geom_tile(aes(
    round_region1, round_region2, fill=mean.logFC),
    data=region_tiles)+
  scale_fill_gradient2()+
  coord_equal()

cluster_dt <- pixel_dt[, {
  wide <- dcast(.SD, region1 ~ region2, length, fill=0)
  m <- as.matrix(wide[,-1])
  clust_id_mat <- cbind(0, rbind(0, m, 0), 0)
  exp_one <- function(x)c(min(x)-1, x, max(x)+1)
  path_list <- contourLines(
    exp_one(wide$region),
    exp_one(as.integer(colnames(m))),
    clust_id_mat, levels=0.5)
  xy <- c('x','y')
  circ_diff_vec <- function(z)diff(c(z,z[1]))
  circ_diff_dt <- function(DT, XY)DT[
  , paste0("d",XY) := lapply(.SD, circ_diff_vec), .SDcols=XY]
  out <- data.table(id=seq_along(path_list))[
  , path_list[[id]][xy]
  , by=id]
  for(XY in list(xy, paste0("d",xy))){
    circ_diff_dt(out, XY)
  }
  out[
  , both_zero := ddx==0 & ddy==0
  ][!c(both_zero[.N], both_zero[-.N]), .(id, region1=x, region2=y)]
}, by=Cluster][]

r1r2_xy_mat <- rbind(
  c(0.5, -0.5),
  c(0.5, 0.5))
corners <- function(corner1, corner2)data.table(corner1, corner2)
set_xy <- function(DT, prefix){
  geti <- function(i)DT[[paste0(prefix,i)]]
  DT[, paste0(
    prefix, "_", c("x","y")
  ) := as.data.table(
    cbind(geti(1), geti(2)) %*% r1r2_xy_mat
  )]
}
get_corners <- function(r1, r2, half.width){
  set_xy(rbind(
    corners(r1-half.width, r2-half.width),
    corners(r1-half.width, r2+half.width),
    corners(r1+half.width, r2+half.width),
    corners(r1+half.width, r2-half.width),
    corners(r1-half.width, r2-half.width)
  ), "corner")
}
expand <- 45
pixel_xy <- pixel_dt[, data.table(.SD, get_corners(region1, region2, expand/100))]
set_xy(pixel_xy, "round_region")
for(xy in c("x","y")){
  xy_val <- pixel_xy[[paste0("corner_",xy)]]
  rxy_val <- pixel_xy[[paste0("round_region_",xy)]]
  set(pixel_xy, j=paste0("rel_",xy), value=xy_val-rxy_val)
}
region_tiles_xy <- region_tiles[, data.table(
  .SD, get_corners(round_region1, round_region2, expand)
)]

ggplot()+
  geom_polygon(aes(
    corner_x, corner_y, fill=mean.logFC, group=round_r1r2),
    data=region_tiles_xy,
    color="black")+
  scale_fill_gradient2()+
  theme_bw()

two_tiles <- pixel_xy[round_r1r2 %in% unique(round_r1r2)[c(1,40)] ]
ggplot()+
  facet_wrap("round_r1r2")+
  geom_polygon(aes(
    rel_x, rel_y, fill=logFC,
    color=neg.log10.p,
    group=rel_regions),
    data=two_tiles)+
  scale_color_gradient(low="white",high="black")+
  scale_fill_gradient2()+
  theme_bw()

animint(
  out.dir="figure-pixels-volcano-interactive-clusters-zoom",
  pixelTiles=ggplot()+
    geom_polygon(aes(
      corner_x, corner_y, fill=clusters, group=round_r1r2),
      data=region_tiles_xy,
      clickSelects="round_r1r2",
      color="black")+
    scale_fill_gradient(low="white",high="black")+
    theme_bw(),
  pixelZoom=ggplot()+
    geom_polygon(aes(
      rel_x, rel_y, fill=logFC, group=r1r2),
      data=pixel_xy,
      showSelected="round_r1r2",
      color="black")+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()
)
## TODO.

ggplot()+
  geom_point(aes(
    logFC, neg.log10.p),
    data=some_pixels)

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
  title="Zoomable Hi-C pixel heat map with borders linked to zoomable volcano plot",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-02/figure-pixels-volcano-interactive-clusters-zoom.R",
  pixels=ggplot()+
    ggtitle("Genomic interaction plot, select pixels")+
    theme_bw()+
    theme_animint(width=1200, colspan=3, last_in_row=TRUE)+
    scale_x_continuous("Position on chr19 (kb)")+
    scale_y_continuous("Distance from position")+
    geom_polygon(aes(
      x, y, fill=logFC, color=neg.log10.p, group=r1r2),
      size=1,
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
  animint2pages(viz, "2025-10-02-HiC-pixels-heat-clusters-zoom", chromote_sleep_seconds=5)
}
viz
