library(data.table)
rmax <- 1000
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
set_xy(pixel_dt, "region")
pixel_corner_dt <- pixel_dt[, data.table(.SD, get_corners(region1, region2, expand.prop))]

get_boundaries <- function(DT){
  dcast_input <- rbind(
    DT[, .(
      region1=seq(min(region1), max(region1)),
      region2=min(region2)-1L,
      value=0
    )],
    DT[, .(
      region1=min(region1)-1L,
      region2=seq(min(region2), max(region2)),
      value=0
    )],
    DT[, .(region1, region2, value=1)])
  wide <- dcast(dcast_input, region1 ~ region2, fill=0)
  m <- as.matrix(wide[,-1])
  clust_id_mat <- cbind(rbind(m, 0), 0)
  exp_one <- function(x)c(x, max(x)+1)
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
  set_xy(out[
  , both_zero := ddx==0 & ddy==0
  ][!c(both_zero[.N], both_zero[-.N]), .(id, region1=x, region2=y)], "region")
}
global_cluster_dt <- pixel_dt[, get_boundaries(.SD), by=Cluster]

myround <- function(x, bin_size=1, offset=0)round((x+offset)/bin_size)*bin_size
off_list <- list(x=20, y=-24)
for(xy in names(off_list)){
  round_fun <- function(rxy)myround(rxy, 50, off_list[[xy]])
  region_xy <- pixel_corner_dt[[paste0("region_",xy)]]
  corner_xy <- pixel_corner_dt[[paste0("corner_",xy)]]
  pixel_dt[, paste0("round_region_",xy) := round_fun(get(paste0("region_",xy)))]
  round_region_xy <- round_fun(region_xy)
  set(pixel_corner_dt, j=paste0("round_region_",xy), value=round_region_xy)
  set(pixel_corner_dt, j=paste0("rel_region_",xy), value=region_xy-round_region_xy)
  set(pixel_corner_dt, j=paste0("rel_corner_",xy), value=round(corner_xy-round_region_xy,2))
  clust_region <- global_cluster_dt[[paste0("region_",xy)]]
  clust_round <- round_fun(clust_region)
  set(global_cluster_dt, j=paste0("round_region_",xy), value=clust_round)
  set(global_cluster_dt, j=paste0("rel_region_",xy), value=clust_region-clust_round)
}
add_round_regions <- function(DT)DT[
, round_regions := paste0(round_region_x,",",round_region_y)
][]
add_round_regions(global_cluster_dt)
add_round_regions(pixel_dt)
add_round_regions(pixel_corner_dt)[, let(
  rel_regions = paste0(rel_region_x,",",rel_region_y)
)][, .(corners=.N), by=rel_regions]

DT <- pixel_dt[Cluster==73 & round_regions=="600,350"]
## TODO investigate why there are Cluster 73 polygons adjacent to each other here.

local_cluster_dt <- pixel_dt[
, get_boundaries(.SD)
, by=.(Cluster,round_region_x,round_region_y,round_regions)]
region_tiles <- pixel_corner_dt[, .(
  pixels=.N,
  mean.neg.log10.p=mean(neg.log10.p),
  mean.logFC=mean(logFC),
  clusters=length(unique(Cluster))
), by=.(round_regions,round_region_x,round_region_y)]
for(xy in names(off_list)){
  local_cluster_dt[, paste0("rel_region_",xy) := get(paste0("region_",xy))-get(paste0("round_region_",xy))][]
  region_tiles[, (xy) := get(paste0("round_region_",xy))-off_list[[xy]]][]
}

ggplot()+
  geom_tile(aes(
    x, y, fill=mean.logFC),
    data=region_tiles)+
  geom_polygon(aes(
    region_x, region_y, group=paste(Cluster,id)),
    data=global_cluster_dt,
    fill="white",
    color="black",
    alpha=0.5)+
  scale_fill_gradient2()

corner_range <- dcast(
  pixel_corner_dt,
  . ~ .,
  list(min,max),
  value.var=c("rel_corner_x","rel_corner_y"))
pixel_dt[, lnl.p := log10(neg.log10.p)][, let(
  round_logFC = myround(logFC,0.25),
  round_neg.log10.p = myround(neg.log10.p, 5),
  round_lnl.p = myround(lnl.p, 0.1)
)][, let(
  volcano_bin = sprintf("logFC=%s -log10(p)=%s", round_logFC, round_neg.log10.p),
  lnl_bin = sprintf("logFC=%s l-l(p)=%s", round_logFC, round_lnl.p),
  relative_logFC = logFC-round_logFC,
  relative_neg.log10.p = neg.log10.p-round_neg.log10.p,
  relative_lnl.p = lnl.p-round_lnl.p
)][]
dcast(pixel_dt, . ~ ., list(min,max,Nvalues=function(x)length(unique(x))), value.var=patterns("round",cols=names(pixel_dt)))

volcano_heat_lnl <- pixel_dt[, .(
  pixels=.N
), by=.(round_logFC,round_lnl.p,lnl_bin)]
volcano_heat_lnl[order(pixels)]
ggplot()+
  geom_tile(aes(
    round_logFC, round_lnl.p, fill=log10(pixels)),
    data=volcano_heat_lnl)+
  geom_text(aes(
    round_logFC, round_lnl.p, label=round(log10(pixels))),
    color="red",
    size=5,
    data=volcano_heat_lnl)+
  scale_fill_gradient(low="white",high="black")+
  theme_bw()

volcano_heat <- pixel_dt[, .(
  pixels=.N
), by=.(round_logFC,round_neg.log10.p,volcano_bin)]
volcano_heat[order(pixels)]
ggplot()+
  geom_tile(aes(
    round_logFC, round_neg.log10.p, fill=log10(pixels)),
    data=volcano_heat)+
  geom_text(aes(
    round_logFC, round_neg.log10.p, label=round(log10(pixels))),
    color="red",
    size=5,
    data=volcano_heat)+
  scale_fill_gradient(low="white",high="black")+
  theme_bw()

## TODO break up huge log10(p)=0 counts into negative space bins.

## TODO select cluster by size (number of pixels / tiles).

## TODO use green points instead of polygon on interaction distance heat map.

## TODO on volcanoHeat, use geom_point(showSelected=c("Cluster","round_regions"),chunk_vars="round_regions") instead of showSelected="Cluster" -- show only pixels in currently selected interaction tile.

count_by_Cluster <- dcast(
  pixel_dt,
  Cluster ~ .,
  list(min, max),
  value.var=c("logFC", "neg.log10.p")
)[
  pixel_dt[, .(
    pixels=.N,
    tiles=length(unique(round_regions))
  ), keyby=Cluster]
]
count_by_Cluster_tile <- pixel_dt[, .(
  displayed_pixels=.N
), keyby=.(Cluster, round_regions)][
  count_by_Cluster
][, let(
  label=sprintf(
    "Cluster %s: %d/%d pixels shown, logFC %.1f to %.1f, log10(p) %.1f to %.1f",
    Cluster, displayed_pixels, pixels,
    logFC_min, logFC_max,
    neg.log10.p_min, neg.log10.p_max),
  rel_x = corner_range[, (rel_corner_x_max+rel_corner_x_min)/2],
  rel_y = corner_range$rel_corner_y_max+1
)]
first_list <- as.list(pixel_dt[which.max(neg.log10.p), .(volcano_bin, Cluster)])
first_list$round_regions <- local_cluster_dt[Cluster==first_list$Cluster, round_regions][1]

table(count_by_Cluster$tiles)
hist(log10(count_by_Cluster$pixels))

cluster_heat_dt <- count_by_Cluster[, let(
  log2.tiles=log2(tiles),
  log10.pixels=log10(pixels)
)][, let(
  round.log2.tiles=myround(log2.tiles, 0.25),
  round.log10.pixels=myround(log10.pixels, 0.25)
)][, let(
  size_bin=sprintf(
    "log2(tiles)=%.2f log10(pixels)=%.2f",
    round.log2.tiles, round.log10.pixels),
  rel.log2.tiles=log2.tiles-round.log2.tiles,
  rel.log10.pixels=log10.pixels-round.log10.pixels
)][, .(
  clusters=.N
), by=.(size_bin, round.log2.tiles, round.log10.pixels)]
ggplot()+
  geom_tile(aes(
    round.log10.pixels, round.log2.tiles,
    fill=log10(clusters)),
    data=cluster_heat_dt)+
  scale_fill_gradient(low="white",high="black")+
  theme_bw()

count_by_Cluster[
, rel_row := 1:.N
, by=.(size_bin, rel.log10.pixels)]
ggplot()+
  geom_point(aes(
    rel.log10.pixels, rel_row),
    data=count_by_Cluster)+
  facet_wrap("size_bin")

cluster.color <- "green"
viz.common <- animint(
  out.dir="figure-pixels-chr1-zoom-tiles-clusters-sizes",
  title="Hi-C pixels chr1 clusters zoom using tiles and sizes",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-09/figure-pixels-chr1-zoom-tiles-clusters-sizes.R",
  first=first_list,
  volcanoHeat=ggplot()+
    ggtitle("Volcano plot summary")+
    geom_tile(aes(
      round_logFC, round_neg.log10.p, fill=log10(pixels)),
      color="grey",
      size=1,
      data=volcano_heat)+
    geom_point(aes(
      logFC, neg.log10.p),
      data=pixel_dt,
      color="black",
      fill=cluster.color,
      chunk_vars="Cluster",
      showSelected="Cluster")+
    geom_tile(aes(
      round_logFC, round_neg.log10.p),
      clickSelects="volcano_bin",
      fill="transparent",
      data=volcano_heat)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()+
    theme_animint(width=300, height=300),
  genomeZoom=ggplot()+
    ggtitle("Genomic interaction zoom")+
    geom_polygon(aes(
      rel_corner_x, rel_corner_y,
      fill=logFC,
      color=neg.log10.p,
      group=rel_regions),
      data=pixel_corner_dt,
      showSelected="round_regions")+
    geom_polygon(aes(
      rel_region_x, rel_region_y, group=paste(Cluster, id)),
      fill="transparent",
      color=cluster.color,
      alpha=1,
      alpha_off=0.2,
      clickSelects="Cluster",
      showSelected="round_regions",
      data=local_cluster_dt)+
    geom_text(aes(
      rel_x, rel_y, label=label, group=1),
      size=20,
      data=count_by_Cluster_tile,
      chunk_vars="round_regions",
      showSelected=c("Cluster","round_regions"))+
    scale_fill_gradient2()+
    scale_color_gradient(
      low="white",high="black",
      guide=guide_legend(override.aes=list(fill="white")))+
    theme_bw()+
    theme_animint(height=800, width=800, rowspan=3),
  clusterHeat=ggplot()+
    ggtitle("Cluster size summary")+
    geom_tile(aes(
      round.log10.pixels, rel_row,
      fill=log10(clusters)),
      color="grey",
      data=cluster_heat_dt)+
    geom_point(aes(
      log10.pixels, log2.tiles),
      color="black",
      fill=cluster.color,
      showSelected="Cluster",
      data=count_by_Cluster)+
    geom_tile(aes(
      round.log10.pixels, rel_row),
      fill="transparent",
      color="black",
      clickSelects="size_bin",
      data=cluster_heat_dt)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()+
    theme_animint(width=300, height=300, last_in_row=TRUE),
  volcanoZoom=ggplot()+
    ggtitle("Volcano plot zoom")+
    geom_point(aes(
      relative_logFC, relative_neg.log10.p),
      data=pixel_dt,
      size=4,
      showSelected="volcano_bin",
      chunk_vars="volcano_bin",
      fill="white",
      color=cluster.color,
      color_off="black",
      clickSelects="Cluster")+
    theme_bw()+
    theme_animint(width=300, height=300),
  clusterZoom=ggplot()+
    ggtitle("Cluster size zoom")+
    geom_point(aes(
      rel.log10.pixels, rel.log2.tiles),
      size=4,
      fill="white",
      color=cluster.color,
      color_off="black",
      showSelected="size_bin",
      clickSelects="Cluster",
      data=count_by_Cluster)+
    scale_fill_gradient(low="white",high="black")+
    theme_bw()+
    theme_animint(width=300, height=300, last_in_row=TRUE),
  genomeSummary=ggplot()+
    ggtitle("Genomic interaction summary")+
    geom_tile(aes(
      x, y, fill=mean.logFC),
      data=region_tiles,
      color="transparent")+
    geom_polygon(aes(
      region_x, region_y, group=paste(Cluster,id)),
      data=global_cluster_dt,
      chunk_vars="Cluster",
      fill="transparent",
      color=cluster.color,
      showSelected="Cluster")+
    geom_tile(aes(
      x, y),
      data=region_tiles,
      clickSelects="round_regions",
      fill="transparent",
      color="black")+
    scale_x_continuous(
      "Genomic bin")+
    scale_y_continuous(
      "Interaction distance")+
    scale_fill_gradient2()+
    theme_bw()+
    theme_animint(width=300, height=300)
)
viz.common

if(FALSE){
  animint2pages(viz.common, "2025-10-09-HiC-pixels-chr1-zoom-tiles-clusters-sizes", chromote_sleep_seconds=5)
}

