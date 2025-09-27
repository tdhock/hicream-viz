library(data.table)
library(ggplot2)
hicream_dt <- fread("hicream_chr19_50000.tsv")[, Cluster := factor(clust)]
## region1 and region2 range: 48092 49257
(clust406 <- hicream_dt[clust==406]) #25 rows.
hicream_dt[, range(clust)] #1 to 1114
hicream_dt[, length(unique(clust))] #1 to 1114
hicream_dt[, .(count=.N), by=.(region1,region2)][, table(count)] ## one row for each combination of region1, region2.

ggplot()+
  geom_tile(aes(
    region1, region2, fill=logFC, color=-log10(p.value)),
    size=3,
    data=clust406)+
  scale_fill_gradient2()+
  scale_color_gradient(low="white",high="black")

ggplot()+
  geom_point(aes(
    logFC, -log10(p.value)),
    data=clust406)

lo <- 48000
hi <- 48100
inside <- function(x)lo < x & x < hi
hicream_zoom <- hicream_dt[inside(region1) & inside(region2)]

ggplot()+
  geom_tile(aes(
    region1, region2, fill=factor(clust)),
    color="black",
    data=hicream_zoom)

band_dt <- hicream_zoom[, {
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
band_simple <- band_dt[keep==TRUE]
lpos_dt <- dcast(hicream_zoom, Cluster ~ ., min, value.var=c("region1","region2"))
data_type_list <- list(
  cluster=band_dt,
  simplified=band_simple,
  pixel=hicream_zoom)
data_geom_list <- list()
data_size_list <- list()
for(data_type in names(data_type_list)){
  data_type_dt <- data_type_list[[data_type]]
  data_geom_list[[data_type]] <- data.table(
    data_type, data_type_dt)
  data_size_list[[data_type]] <- data.table(
    data_type, data_type_dt[, .(rows=.N), by=Cluster])
}
data_size <- rbindlist(data_size_list)
label_dt <- data_size[, .(
  label=paste(paste(c(
    paste("cluster", Cluster),
    paste(data_type, "rows=", rows)
  ), collapse="\n"), "corners=", rows[data_type=="pixel"]*4)
), by=Cluster][lpos_dt, on="Cluster"]
poly_dt <- with(data_geom_list, rbind(simplified, cluster))
set.seed(1)
gg <- ggplot()+
  geom_tile(aes(
    region1, region2,
    fill=Cluster,
    color=data_type,
    size=data_type),
    alpha=0.5,
    data=data_geom_list$pixel)+
  geom_polygon(aes(
    y, x,
    group=paste(Cluster, data_type),
    fill=Cluster,
    color=data_type,
    size=data_type),
    alpha=0.5,
    data=poly_dt)+
  scale_color_manual(values=c(
    pixel="black",
    cluster="black",
    simplified="red"))+
  geom_point(aes(
    y, x,
    size=data_type,
    color=data_type),
    data=poly_dt[.N:1])+
  coord_equal()+
  geom_label(aes(
    region1, region2, label=label),
    alpha=0.7,
    hjust=0,
    size=2.9,
    data=label_dt)+
  guides(fill="none")+
  scale_size_manual(values=c(
    simplified=1,
    cluster=2,
    pixel=0.5))+
  theme(legend.position=c(0.8,0.2))
print(gg)
png("isoband.png", width=5, height=5, units="in", res=200)
print(gg)
dev.off()
