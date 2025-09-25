library("animint2")
library("data.table")
options(scipen=999)
chromosome <- 19
resolution <- 50000
out_elements_fake <- fread(paste0("fake_chr", chromosome, "_", resolution, ".tsv"))
out_elements_all <- fread(paste0("all_chr", chromosome, "_", resolution, ".tsv"))
df_crit <- fread(paste0("df_crit", chromosome, "_", resolution, ".tsv"))

(out_elements_some <- out_elements_all[y<100])
dcast(
  out_elements_some,
  nb_clust ~ .,
  list(length, n_uniq=function(x)length(unique(id))),
  value.var="id"
)[, rows_per_id := id_length/id_n_uniq][]
ggplot() + 
  geom_polygon(
    data = out_elements_some,
    showSelected = "nb_clust", 
    aes(x = x, y = y, group= id, fill=IFF), 
    colour=NA) +
  facet_grid(nb_clust ~ .)+
  scale_fill_gradient(low="white",high="black")

divisor <- 10
out_elements_some[, let(
  round_x=round(x/divisor),
  round_y=round(y/divisor)
)][, let(
  rel_x=x-min(x),
  rel_y=y-min(y),
  tile_id=paste(round_x,round_y)
), by=.(nb_clust, round_x, round_y)]
out_elem_mean <- out_elements_some[, .(
  mean_IFF=mean(IFF)
), by=.(nb_clust, round_x, round_y, tile_id)]
ggplot()+
  geom_tile(aes(
    round_x, round_y, fill=mean_IFF),
    data=out_elem_mean)+
  scale_fill_gradient(low="white",high="black")+
  facet_grid(nb_clust ~ .)

out_first <- out_elements_some[, .SD[1], by=.(nb_clust, id)]
dcast(out_first, nb_clust + tile_id ~., list(min, max), value.var=c("rel_x","rel_y"))
viz <- animint(
  title="Hicream Exploration",
  scatter = ggplot() +
    ggtitle("1. Select number of clusters") +
    geom_point(aes(
      x=nb_clust, y=crit),
      shape=21,
      size=4,
      data=df_crit)  +
    make_tallrect(df_crit, "nb_clust") +
    scale_x_continuous("nb_clust", labels = as.character(df_crit$nb_clust), breaks = df_crit$nb_clust)+
    theme_animint(width=400, height=400),
  hic = ggplot() +
    ggtitle("3. Details of selected tile")+
    geom_polygon(
      data = out_elements_some,
      showSelected = c("nb_clust", "tile_id"),
      aes(rel_x, rel_y, group= id, fill=IFF),
      chunk_vars=character(),
      colour=NA) +
    scale_fill_gradient(low="white",high="black")+
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    ) + 
    theme_animint(width=400, height=400),
  tiles=ggplot()+
    ggtitle("2. Select tile")+
    theme_animint(width=1200, height=250)+
    geom_tile(aes(
      round_x, round_y, fill=mean_IFF),
      showSelected="nb_clust",
      clickSelects="tile_id",
      data=out_elem_mean)+
    scale_fill_gradient(low="white",high="black"),
  out.dir="plot_hicream",
  first=list(
    tile_id="31 1",
    nb_clust=614),
  source="https://github.com/tdhock/hicream-viz/blob/main/example_animint2_hicream/plot_2025-09-25.R"
)
if(FALSE){
  animint2pages(viz, "2025-09-25-hicream")
}


