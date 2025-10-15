library(data.table)
rmax <- 1000
pixel_dt <- fread("hicream_chr1_50000.tsv")[, let(
  Cluster=factor(clust),
  neg.log10.p = -log10(p.value)
)][
  region1<=rmax & region2<=rmax
]
library(ggplot2)

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
show_Cluster <- 73
DT <- pixel_dt[Cluster==show_Cluster & round_regions==show_rr]
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
ggplot()+
  geom_tile(aes(
    region1, region2, fill=factor(value)),
    data=dcast_input)

wide <- dcast(dcast_input, region1 ~ region2, fill=0)
m <- as.matrix(wide[,-1])
isoband::plot_iso(m, 0.5, 1.5)

clust_id_mat <- cbind(rbind(m, 0), 0)
clust_id_dt <- data.table(
  row=as.integer(row(clust_id_mat)),
  col=as.integer(col(clust_id_mat)),
  value=as.numeric(clust_id_mat))
ggplot()+
  theme_bw()+
  geom_tile(aes(
    col, row, fill=factor(value)),
    data=clust_id_dt)+
  scale_y_reverse()+
  coord_equal()+
  scale_fill_manual(values=c(
    "1"="grey50",
    "0"="white"))

exp_one <- function(x)c(x, max(x)+1)
path_list <- contourLines(
  exp_one(wide$region),
  exp_one(as.integer(colnames(m))),
  clust_id_mat, levels=0.5)
band_list <- isoband::isobands(
  exp_one(as.integer(colnames(m))),
  exp_one(wide$region), clust_id_mat, 0.5, 1.5)
xy <- c('x','y')
out_list <- list(
  contourLines=data.table(id=seq_along(path_list))[
  , path_list[[id]][xy]
  , by=id],
  isobands=with(band_list[[1L]], data.table(id, x=y, y=x)))
circ_diff_vec <- function(z)diff(c(z,z[1]))
circ_diff_dt <- function(DT, XY)DT[
, paste0("d",XY) := lapply(.SD, circ_diff_vec), .SDcols=XY]
out_dt_list <- list()
for(fun in names(out_list)){
  out <- out_list[[fun]]
  for(XY in list(xy, paste0("d",xy))){
    circ_diff_dt(out, XY)
  }
  out_dt_list[[fun]] <- out[
  , both_zero := ddx==0 & ddy==0
  ][!c(both_zero[.N], both_zero[-.N]), data.table(fun,id,x,y)]
}
(out_dt <- rbindlist(out_dt_list))

gg <- ggplot()+
  ggtitle("holes work with ggplot2 and aes(subgroup)")+
  geom_tile(aes(
    region2, region1, fill=factor(value)),
    color=NA,
    data=dcast_input)+
  geom_polygon(aes(
    y, x, group=1, subgroup=id),
    alpha=0.2,
    fill="black",
    color="black",
    data=out_dt)+
  scale_y_reverse()+
  facet_grid(. ~ fun, labeller=label_both)
png("figure-pixels-chr1-Cluster73-ggplot2.png", width=7, height=3, units="in", res=200)
print(gg)
dev.off()
