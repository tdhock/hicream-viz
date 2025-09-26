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

band_dt <- hicream_zoom[, {
  wide <- dcast(.SD, region1 ~ region2, length, fill=0)
  m <- as.matrix(wide[,-1])
  clust_id_mat <- cbind(0, rbind(0, m, 0), 0)
  exp_one <- function(x)c(min(x)-1, x, max(x)+1)
  band_list <- isoband::isobands(
    exp_one(as.integer(colnames(m))),
    exp_one(wide$region), clust_id_mat, 0.5, 1.5)
  as.data.table(band_list[[1]][c("x","y")])
}, by=clust]
ggplot()+
  facet_grid(. ~ geom)+
  geom_tile(aes(
    region1, region2, fill=factor(clust)),
    color="black",
    data=hicream_zoom)+
  geom_polygon(aes(
    y, x, group=clust, fill=factor(clust)),
    color="black",
    size=2,
    data=data.table(band_dt, geom="polygon"))


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
