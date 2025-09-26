library(data.table)
library(ggplot2)
hicream_dt <- fread("hicream_chr19_50000.tsv")
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

levs <- sort(unique(hicream_zoom$clust))
levs <- 814
levs <- 608
levs <- 406#big
levs <- 530#one pixel
hicream_zoom_wide <- dcast(hicream_zoom, region1 ~ region2, value.var="clust", fill=0)
clust_id_mat <- cbind(0, rbind(
  0, as.matrix(hicream_zoom_wide[,-1]), 0), 0)
reg_int <- hicream_zoom_wide[, c(min(region1)-1, region1, max(region1)+1)]
band_list <- isoband::isobands(reg_int, reg_int, clust_id_mat, levs-0.5, levs+0.5)
line_list <- isoband::isolines(reg_int, reg_int, clust_id_mat, levs)
names(band_list)
band_dt <- data.table(clust=levs)[, {
  as.data.table(band_list[[.I]])
}, by=clust]
line_dt <- data.table(clust=levs)[, {
  as.data.table(line_list[[.I]])
}, by=clust]
ggplot()+
  facet_grid(. ~ geom)+
  geom_tile(aes(
    region1, region2, fill=factor(clust)),
    color="black",
    data=hicream_zoom)+
  geom_polygon(aes(
    y, x, group=id, fill=factor(clust)),
    color="black",
    size=2,
    data=data.table(band_dt, geom="polygon"))+
  geom_path(aes(
    y, x, group=id),
    color="black",
    size=2,
    data=data.table(line_dt, geom="line"))+
  coord_equal()


all_dt <- fread("all_chr19_50000.tsv")
all_dt[, table(nb_clust)]

only1114 <- all_dt[nb_clust==1114]
all_dt[nb_clust==114]
range(all_dt$i)
range(all_dt$j)


hist(hicream_dt$logFC)
hicream_dt[, hist(-log10(p.value), 100)]
hicream_dt[p.value>0.001, hist(-log10(p.value))]

tile_dt <- hicream_dt[
, neg.log.p := -log10(p.value)
][, let(
  round_neg.log.p = round(neg.log.p),
  round_logFC = round(logFC)
)][, .(
  log10.pixels=log10(.N)
), by=.(round_neg.log.p, round_logFC)]

ggplot()+
  geom_tile(aes(
    round_logFC, round_neg.log.p,
    fill=log10.pixels),
    data=tile_dt)+
  scale_fill_gradient(low="white",high="black")
    
old[id==406]
range(all_dt$id)
all_dt[id==406]
