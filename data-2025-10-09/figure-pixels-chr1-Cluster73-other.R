library(data.table)
remotes::install_github("animint/animint2@fix/252-polygon-holes-subgroup")
##remotes::install_github("animint/animint2")

rmax <- 1000
pixel_dt <- fread("hicream_chr1_50000.tsv")[, let(
  Cluster=factor(clust),
  neg.log10.p = -log10(p.value)
)][
  region1<=rmax & region2<=rmax
]
library(animint2)

expand <- 45
expand.prop <- expand/100
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
set_xy(pixel_dt, "region")

myround <- function(x, bin_size=1, offset=0)round((x+offset)/bin_size)*bin_size
off_list <- list(x=20, y=-24)
for(xy in names(off_list)){
  round_fun <- function(rxy)myround(rxy, 50, off_list[[xy]])
  pixel_dt[, paste0("round_region_",xy) := round_fun(get(paste0("region_",xy)))]
}
add_round_regions <- function(DT)DT[
, round_regions := paste0(round_region_x,",",round_region_y)
][]
add_round_regions(pixel_dt)

show_rr <- "600,350"
show_pixels <- pixel_dt[round_regions==show_rr]
ggplot()+
  geom_tile(aes(
    region1, region2, fill=Cluster),
    data=show_pixels)

local_clust <- show_pixels[, get_boundaries(.SD), by=Cluster]

viz <- animint(
  title="aes(subgroup) works for HiC data",
  source="https://github.com/tdhock/hicream-viz/blob/main/data-2025-10-09/figure-pixels-chr1-Cluster73-other.R",
  ggplot()+
    geom_polygon(aes(
      region1, region2,
      group=Cluster, subgroup=id,
      fill=Cluster, tooltip=Cluster),
      size=1,
      data=local_clust,
      color="black"))

if(FALSE){
  animint2pages(viz, "2026-05-22-HiC-pixels-chr1-Cluster73-subgroup", chromote_sleep_seconds=5)
}
