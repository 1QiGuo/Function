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

# Heatmap
```{r}
#load library
library(RColorBrewer)
library("pheatmap")
libary("ggplot2")

#load data
load("/users/PAS1475/liuzeyi/guoqi/output/sample6_clusters.RData")
#express data
setwd("/users/PAS1475/liuzeyi/guoqi/output/AD_Control_DEG")
deg <- read.csv("Layer marker for heatmap.csv")
expr_data <- as.data.frame(sample6.combined@assays$SCT@data[deg$Marker,])
expr_data<-rbind(expr_data,layer=sample6.combined$Layer)
expr_data<-as.data.frame(t(expr_data))
temp_data<-data.frame(rep(0,7))
for(j in 1:21){
  gene_data<-c()
  for(i in sort(unique(expr_data$layer))){
    temp_gene<-mean(as.numeric(expr_data[which(expr_data$layer==i),j]))
    gene_data<-c(gene_data,temp_gene)
  }
  temp_data<-cbind(temp_data,gene_data)
}
temp_data<-temp_data[,-1]
rownames(temp_data)<-sort(unique(expr_data$layer))
colnames(temp_data)<-deg$Marker
heatmap_data<-t(temp_data)
#plot
my_colour = list(
  gene_cat = c( Layer_1 = brewer.pal(n = 7, name = "Dark2")[1],
                Layer_2 = brewer.pal(n = 7, name = "Dark2")[2],
                Layer_3 = brewer.pal(n = 7, name = "Dark2")[3],
                Layer_4 = brewer.pal(n = 7, name = "Dark2")[4],
                Layer_5 = brewer.pal(n = 7, name = "Dark2")[5],
                Layer_6 = brewer.pal(n = 7, name = "Dark2")[6],
                White_Matter = brewer.pal(n = 7, name = "Dark2")[7]))
gene_cat_val <- as.numeric(table(deg$Layer))
gene_anno <- data.frame(gene_cat = c(rep("Layer_1",as.numeric(table(deg$Layer)[1])),
                                     rep("Layer_2",as.numeric(table(deg$Layer)[2])),
                                     rep("Layer_4",as.numeric(table(deg$Layer)[3])),
                                     rep("Layer_5",as.numeric(table(deg$Layer)[4])),
                                     rep("Layer_6",as.numeric(table(deg$Layer)[5])),
                                     rep("White_Matter",as.numeric(table(deg$Layer)[6]))
) )

p <- pheatmap(heatmap_data,
              gaps_row = cumsum(gene_cat_val),
              cluster_rows = F,
              cluster_cols = F,
              scale = "row",
              border_color = NA,
              annotation_colors = my_colour, 
              #annotation_row = gene_anno,
              labels_col =rownames(temp_data),
              color = colorRampPalette(c("blue","white","red"))(100))

ggsave(
  plot = p,
  filename = "/users/PAS1475/liuzeyi/guoqi/output/picture/integrating landscape/heatmap_revised.tiff",
  device = "tiff",
  dpi = 150,
  width = 7,
  height = 11,
  units = "in"
)
```

# Dot plot for enrichment analysis
```{r}
library(ggplot2)
setwd("/users/PAS1475/liuzeyi/guoqi/output/picture/module_4_enrichment_dot")
module4<-list()
for(i in 1:4){
  temp<-read_excel("No.10, 4 modules GO pathway.xlsx", sheet = i)
  temp_sort<-temp[with(temp, order(Count)), ]
  temp_sort$Description<-factor(temp_sort$Description,levels=(temp_sort$Description))
  count_sum<-as.numeric(unlist(strsplit(temp_sort$GeneRatio[1],split="/"))[2])
  temp_sort$GeneRatio_num<-c(temp_sort$Count)/count_sum
  module4[[i]]<- temp_sort
}
names(module4)<-excel_sheets(path = "No.10, 4 modules GO pathway.xlsx")
for(i in 1:4){
  g<-ggplot(module4[[i]], # you can replace the numbers to the row number of pathway of your interest
            aes(x = GeneRatio_num, y = Description)) + 
    geom_point(aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 14) +
    theme(axis.text = element_text(size = 14, face = "bold"),
    )+
    scale_colour_gradient(limits=NULL, low="red",high="blue") +
    ylab(NULL) +
    ggtitle("GO enrichment")
  ggsave(
    plot = g,
    filename = paste0(names(module4)[i],"_enrichment.tiff"),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
}
```
# bar plot
## geom_bar
highest is 1
cell components
```{r}
#data
cell.prop<-as.data.frame(prop.table(table(Idents(heart_all), heart_all$time)))
#plot
ggplot(cell.prop,aes(time,proportion,fill=cluster))+
  
  geom_bar(stat="identity",position="fill")+
  
  ggtitle("")+
  
  theme_bw()+
  
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))
```
## geom_col
highest number do not have limitation
umi propotion
```{r}
ggplot(cell.umi,aes(time,cell_umi,fill=cell_type))+
  geom_col(aes(fill=cell_type))
```
