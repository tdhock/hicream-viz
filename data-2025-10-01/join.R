library(data.table)
add_ij <- function(DT)setkey(DT[, ij := paste0(i,",",j)], i, j, ij)
if(FALSE){
  old_poly_dt <- add_ij(fread("../data-2025-09-26/all_chr19_50000.tsv"))
  old_poly_dt[ij=='2,1' & nb_clust==1114]
}
all_poly_dt <- add_ij(fread("all_chr19_50000.tsv"))
all_poly_dt[1]
pixel_dt <- fread("../data-2025-09-26/hicream_chr19_50000.tsv")
pixel_dt[1]
pixel_clusters <- length(unique(pixel_dt$clust))
poly_dt <- all_poly_dt[nb_clust==pixel_clusters][, .SD[id==id[1]], by=.(i,j,ij)]
pixel_region_min <- pixel_dt[, min(region1,region2)]
add_ij(pixel_dt[, let(
  i = region2-pixel_region_min+1L,
  j = region1-pixel_region_min+1L
)])
two_pixels <- pixel_dt[1:2]
two_pixels
two_polys <- poly_dt[two_pixels]
two_polys

pixels_for_two_clusters <- pixel_dt[clust %in% unique(clust)[1:2] ]
polys_for_two_clusters <- poly_dt[pixels_for_two_clusters]

bands_for_two_clusters <- polys_for_two_clusters[, {
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



join_dt <- pixel_dt[poly_dt]
join_dt[is.na(chr)]
