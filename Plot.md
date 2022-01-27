# Spatial plot
```{r}
#unsupervised spatial
plot_spatial_map <-
  function(object = Seurat.object,
           group.by = "seurat_clusters",
           pt.size = 4.8) {
    my.cat <- unique(object$patientID)
    tmp.slice.name <- switch (
      my.cat,
      '2_5' = "slice1",
      '18_64' = "slice1.1",
      '2_8' = "slice1.2",
      'T4857' = "slice1.3",
      '1_1' = "slice1.4",
      '2_3' = "slice1.5"
    )
    xy.cor <-  object@images[[tmp.slice.name]]
    xy.cor <- xy.cor@coordinates
    my.meta <- object@meta.data[, group.by]
    plot.df <- as.data.frame(cbind(xy.cor, group.by = my.meta))
    p <-
      ggplot(plot.df, aes(
        x = imagecol,
        y = -imagerow,
        color = group.by
      )) +
      geom_point(size = pt.size)  + theme_void()
    if (any(grepl("noise", unique(my.meta), ignore.case = T))) {
      p <-
        p + scale_color_manual(
          values = c(as.character(jcolors::jcolors(palette = "pal7"))[-8], "#808080"),
          breaks = c(
            "Layer 1",
            "Layer 2",
            "Layer 3",
            "Layer 4",
            "Layer 5",
            "Layer 6",
            "White Matter",
            "Noise"
          )
        ) +
        labs(col = group.by, title = group.by)
    } else{
      p <-
        p + scale_color_manual(values = as.character(jcolors::jcolors(palette = "pal7")),
                               breaks = sort(unique(my.meta))) +
        labs(col = group.by, title = group.by)
    }
    return(p)
  }
# prepare object
table(sample6.combined$patientID)
obj.list <- SplitObject(sample6.combined, split.by = "patientID")
table(obj.list[[1]]$seurat_clusters)

p.2_5_un <-
  plot_spatial_map(object = obj.list[[1]], group.by = "seurat_clusters")
p.18_64_un <-
  plot_spatial_map(object = obj.list[[2]], group.by = "seurat_clusters")
p.2_8_un <-
  plot_spatial_map(object = obj.list[[3]], group.by = "seurat_clusters")
p.T4857_un <-
  plot_spatial_map(object = obj.list[[4]], group.by = "seurat_clusters")
p.1_1_un <-
  plot_spatial_map(object = obj.list[[5]], group.by = "seurat_clusters")
p.2_3_un <-
  plot_spatial_map(object = obj.list[[6]], group.by = "seurat_clusters")
# obj.list[[1]]@images[setdiff(name.list, i)]
# name.list <- seq(1:8)
# p.2_3 <- plot_spatial_map(object = obj.list[[6]], group.by = "Layer")
plot.un <-
  list(p.2_5_un, p.18_64_un, p.2_8_un, p.T4857_un, p.1_1_un, p.2_3_un)
names(plot.un) <- names(obj.list)

for (i in 1:6) {
  ggsave(
    plot = plot.un[[i]],
    filename = paste0(
      "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/",
      names(plot.un)[i],
      "_un.tiff"
    ),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}
```
