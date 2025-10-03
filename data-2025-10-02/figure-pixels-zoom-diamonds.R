library(animint2)
devtools::load_all("~/R/animint2")
library(data.table)
myround <- function(x,value=1)round(x*value)/value
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")[, let(
  Cluster=factor(clust),
  neg.log10.p = -log10(p.value)
)]
denom <- 50
for(i in 1:2){
  region_i <- pixel_dt[[paste0("region",i)]]
  round_region_i <- myround(region_i, 1/denom)
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
expand.prop <- expand/100
pixel_xy <- pixel_dt[, data.table(.SD, get_corners(region1, region2, expand.prop))]
set_xy(pixel_xy, "round_region")
for(xy in c("x","y")){
  xy_val <- pixel_xy[[paste0("corner_",xy)]]
  rxy_val <- pixel_xy[[paste0("round_region_",xy)]]
  set(pixel_xy, j=paste0("rel_",xy), value=xy_val-rxy_val)
}
region_tiles_xy <- region_tiles[, data.table(
  .SD, get_corners(round_region1, round_region2, diff(unique(pixel_dt$round_region1))[1]*expand.prop)
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

viz.common <- animint(
  out.dir="figure-pixels-zoom-common",
  title="Hi-C pixels zoom using diamonds",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-02/figure-pixels-zoom.R",
  pixelTiles=ggplot()+
    geom_polygon(aes(
      corner_x, corner_y, fill=mean.logFC, group=round_r1r2),
      data=region_tiles_xy,
      clickSelects="round_r1r2",
      color="black")+
    scale_fill_gradient2()+
    theme_bw()+
    theme_animint(height=300),
  pixelZoom=ggplot()+
    geom_polygon(aes(
      rel_x, rel_y,
      fill=logFC,
      color=neg.log10.p,
      group=rel_regions),
      data=pixel_xy,
      showSelected="round_r1r2")+
    scale_fill_gradient2()+
    scale_color_gradient(
      low="white",high="black")+
      ##guide=guide_legend(override.aes=list(fill="white")))+
    theme_bw()+
    theme_animint(height=800, width=800))
animint2dir(viz.common)

viz.no.common <- animint(
  out.dir="figure-pixels-zoom-no-common",
  title="Hi-C pixels zoom using diamonds",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-02/figure-pixels-zoom.R",
  pixelTiles=ggplot()+
    geom_polygon(aes(
      corner_x, corner_y, fill=mean.logFC, group=round_r1r2),
      data=region_tiles_xy,
      clickSelects="round_r1r2",
      color="black")+
    scale_fill_gradient2()+
    theme_bw()+
    theme_animint(height=300),
  pixelZoom=ggplot()+
    geom_polygon(aes(
      rel_x, rel_y,
      fill=logFC,
      color=neg.log10.p,
      group=rel_regions),
      data=pixel_xy,
      showSelected="round_r1r2")+
    scale_fill_gradient2()+
    scale_color_gradient(
      low="white",high="black",
      guide=guide_legend(override.aes=list(fill="white")))+
    theme_bw()+
    theme_animint(height=800, width=800))
viz.no.common

if(FALSE){
  animint2pages(viz, "2025-10-02-HiC-pixels-zoom-diamonds", chromote_sleep_seconds=5)
}
viz
