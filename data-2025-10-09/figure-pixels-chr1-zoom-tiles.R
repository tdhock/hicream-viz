library(data.table)
rmax <- 3200
pixel_dt <- fread("hicream_chr1_50000.tsv")[, let(
  Cluster=factor(clust),
  neg.log10.p = -log10(p.value)
)][
  ##region1<=rmax & region2<=rmax
]
library(animint2)

## every region1,2 bin is 50kb.
r1r2_xy_mat <- rbind(
  c(0.5, -0.5),
  c(0.5, 0.5))
set_xy <- function(DT, prefix){
  geti <- function(i)DT[[paste0(prefix,i)]]
  DT[, paste0(
    prefix, "_", c("x","y")
  ) := as.data.table(
    cbind(geti(1), geti(2)) %*% r1r2_xy_mat
  )][]
}
get_corners <- function(r1, r2, half.width){
  corners <- function(d1, d2)data.table(corner1=r1+d1, corner2=r2+d2)
  set_xy(rbind(
    corners(-half.width, -half.width),
    corners(-half.width, half.width),
    corners(half.width, half.width),
    corners(half.width, -half.width),
    corners(-half.width, -half.width)
  ), "corner")
}
expand <- 45
expand.prop <- expand/100
pixel_xy <- set_xy(
  pixel_dt[, data.table(.SD, get_corners(region1, region2, expand.prop))],
  "region")

myround <- function(x, bin_size=1, offset=0)round((x+offset)/bin_size)*bin_size
off_list <- list(x=20, y=-25)
for(xy in names(off_list)){
  region_xy <- pixel_xy[[paste0("region_",xy)]]
  corner_xy <- pixel_xy[[paste0("corner_",xy)]]
  round_region_xy <- myround(region_xy, 50, off_list[[xy]])
  set(pixel_xy, j=paste0("round_region_",xy), value=round_region_xy)
  set(pixel_xy, j=paste0("rel_region_",xy), value=region_xy-round_region_xy)
  set(pixel_xy, j=paste0("rel_corner_",xy), value=round(corner_xy-round_region_xy,2))
}
pixel_xy[, let(
  rel_regions = paste0(rel_region_x,",",rel_region_y),
  round_regions = paste0(round_region_x,",",round_region_y)
)][, .(corners=.N), by=rel_regions]
region_tiles <- pixel_xy[, .(
  pixels=.N,
  mean.neg.log10.p=mean(neg.log10.p),
  mean.logFC=mean(logFC),
  clusters=length(unique(Cluster))
), by=.(round_regions,round_region_x,round_region_y)]

ggplot()+
  geom_tile(aes(
    round_region_x, round_region_y, fill=mean.logFC),
    data=region_tiles)+
  scale_fill_gradient2()+
  coord_equal()

two_tiles <- pixel_xy[round_regions %in% unique(round_regions)[1:2] ]
ggplot()+
  geom_polygon(aes(
    rel_corner_x, rel_corner_y,
    fill=logFC,
    color=neg.log10.p,
    group=rel_regions),
    data=two_tiles,
    showSelected="round_regions")+
  scale_fill_gradient2()+
  scale_color_gradient(
    low="white",high="black",
    guide=guide_legend(override.aes=list(fill="white")))+
  theme_bw()+
  facet_wrap("round_regions")

devtools::load_all("~/R/animint2")
system.time({
  viz.common <- animint(
    out.dir="figure-pixels-chr1-zoom-tiles",
    title="Hi-C pixels chr1 zoom using tiles",
    source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-09/figure-pixels-chr1-zoom-tiles.R",
    pixelTiles=ggplot()+
      geom_tile(aes(
        round_region_x, round_region_y, fill=mean.logFC),
        data=region_tiles,
        clickSelects="round_regions",
        color="black")+
      scale_fill_gradient2()+
      theme_bw()+
      theme_animint(height=300),
    pixelZoom=ggplot()+
      geom_polygon(aes(
        rel_corner_x, rel_corner_y,
        fill=logFC,
        color=neg.log10.p,
        group=rel_regions),
        data=pixel_xy,
        showSelected="round_regions")+
      scale_fill_gradient2()+
      scale_color_gradient(
        low="white",high="black",
        guide=guide_legend(override.aes=list(fill="white")))+
      theme_bw()+
      theme_animint(height=800, width=800))
  print(viz.common)
  system("du -ms figure-pixels-chr1-zoom-tiles")
})


if(FALSE){
  animint2pages(viz.common, "2025-10-09-HiC-pixels-chr1-zoom-tiles", chromote_sleep_seconds=5)
}

