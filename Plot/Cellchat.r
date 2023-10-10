#------------------------------------------------Take N0 for example
#load packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(igraph)
library(Hmisc)
library(hash)
library(tidyverse)
library(circlize)
library(ggsci)
library(gtools)
library(ComplexHeatmap)
setwd("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/")

#load cellchat data
CellChatDB <- CellChatDB.mouse 
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo

#load our lung data
lung <- readRDS("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Data/lung.rds")
split_lung <- SplitObject(lung, split.by = "orig.ident")
N0 <-split_lung$N0
DefaultAssay(N0) <- 'RNA'
Idents(N0)<- 'main.cluster'
levels(N0)
#count data
counts_N0 <- GetAssayData(N0, assay = "RNA", slot = "counts") # raw counts
#write.table(counts_N0, file = 'counts_N0.csv', sep =',')
#normalize counts for cell chat
counts_N0 <-normalizeData(counts_N0)
#meta data
meta_N0 <-N0@meta.data
#not all groups have Mixed Endo pop, so need to exclude mismatched levels
meta_N0$main.cluster <- droplevels(meta_N0$main.cluster, exclude = setdiff(levels(meta_N0$main.cluster),unique(meta_N0$main.cluster)))
#Please note that blanks are not allowed when plotting, such as "T cells". Hence, I replaced all blanks in cell types with dashes (-) ("T-cells"). 
#Feel free to switch to alternatives or modify the figure legend later as needed.

#switch blank " " to "-"
meta_N0$main.cluster.qi<-gsub(" ", "-", meta_N0$main.cluster)
table(meta_N0$main.cluster.qi)


#-------------the function to use CellChat
use_CellChat <- function(counts, meta,celltype_column) {
  # counts : the gene expression matrix
  # meta   : the meta data matched with the gene expression matrix
  # celltype_column.  : cell type column in the metadata
  counts <- as.matrix(counts) # transform data frame into matrix
  cellchat <- createCellChat(object = counts, meta = meta, group.by = celltype_column) # create a CellChat object
  levels(cellchat@idents) # check the cell types
  CellChatDB <- CellChatDB.mouse # set the database according to the organism
  dplyr::glimpse(CellChatDB$interaction) # show the structure of the database
  unique(CellChatDB$interaction$annotation) # check the annotation of interactions
  CellChatDB.use <- CellChatDB # use all databases
  cellchat@DB <- CellChatDB.use # set the used database in the object
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = T, nboot = 1, population.size = T,trim = 0.1,type = "truncatedMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat) # get the pairs of cell types or ligand receptors in the interactions
  cell_type_pair <- paste0(df.net$source, '_', df.net$target) # build the pairs of interaction cell types
  mean_exp <- sapply(1:nrow(df.net), function(x) {
    ct1 <- df.net$source[x] # source cell typeQ
    ct2 <- df.net$target[x] # target cell type
    cells1 <- meta[meta[[celltype_column]] == ct1,][[1]] # cells belonging to the source cell type
    cells2 <- meta[meta[[celltype_column]] == ct2,][[1]] # cells belonging to the target cell type
    gene_pair1 <- df.net$ligand[x] # genes corresponding to the ligand
    genes1<-gene_pair1
    #genes1 <- capitalize(tolower(strsplit(gene_pair1, split = '_')[[1]]))
    gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    genes2<-gene_pair2
    # genes1 <- capitalize(tolower(strsplit(gene_pair1, split = '_')[[1]]))
    # gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    # genes2 <- capitalize(tolower(strsplit(gene_pair2, split = '_')[[1]]))
    if (length(genes2) > 1 & genes2[2] == 'R2') {
      genes2[2] <- gsub(pattern = '1', replacement = '2', genes2[1])
    }
    cat("Ligand genes: ", genes1, "\n")
    cat("Receptor genes: ", genes2, "\n\n")
    mean1 <- mean(counts[genes1, cells1]) # the mean expression value of the ligand within the source cell type
    mean2 <- mean(counts[genes2, cells2]) # the mean expression value of the receptor within the source cell type
    avg <- (mean1 + mean2) / 2 # the mean expression value of the ligands and receptors
    return(avg)
  })
  df.net <- cbind(cell_type_pair, mean_exp, df.net) # add the mean expression and cell type pairs
  return(df.net)
}
#------------the function to save the adjacency matrix and L-R categories
save_data <- function(dir,net.df, prefix) {
  # dir.    : directory for the results
  # net.df  : the data frame to save the CCIs
  # prefix  : the prefix of the files to output
  setwd(dir)
  node_v1 <- vector() # the first node
  node_v2 <- vector() # the second node
  weight_v <- vector() # the edge weight
  
  for (i in 1:nrow(net.df)) {
    cci <- strsplit(net.df$interaction_name_2[i], split = ' ')[[1]] # split the interaction name
    ct.pair <- strsplit(net.df$cell_type_pair[i], split = '_')[[1]] # get the cell type pair
    ct1 <- ct.pair[1]
    ct2 <- ct.pair[2]
    prefix1 <- ct1
    prefix2 <- ct2
    node1 <- paste0(prefix1, '_', cci[1])
    id2 <- cci[length(cci)]
    weight <- net.df$prob[i]
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])))
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1))
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2)
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) {
    return()
  }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  write.csv(x = as.data.frame(as.matrix(adj)), file = paste0(prefix, '_mat.csv'))
  
  Genes <- names(V(g))
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    return(x[1])
  })
  category <- data.frame(Genes = Genes, ID = ID)
  write.csv(x = category, file = paste0(prefix, '_cat.csv'))
}
#-------------the function to plot a chord diagram for a single datasete
genChord <- function(mat, cat, prefix) {
  # mat    : the adjacency matrix for the CCIs
  # cat   : the category for each gene, i.e., cell type
  # prefix : the prefix of the figure and legend files
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  # read data
  graph_adj <- read.csv(mat, row.names = 1)
  graph_module <- read.csv(cat, row.names = 1)
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  graph_module <- graph_module %>%
    mutate(color = as_factor(ID))%>%  
    mutate(color = fct_recode(color ,
                              "#5A5156FF" = "Monocytes",
                              "#E4E1E3FF" = "T-Cells",
                              "#F6222EFF" = "Type-1-Neutrophils",
                              "#FE00FAFF" = "Type-2-Neutrophils",
                              "#16FF32FF" = "B-Cells",
                              "#3283FEFF" = "General-Capillary",
                              "#FEAF16FF" = "Capillary-Aerocytes",
                              "#B00068FF" = "Myofibroblasts",
                              "#1CFFCEFF" = "Fibroblasts",
                              "#90AD1CFF" = "Vein-Endothelia",
                              "#2ED9FFFF" = "Basophils",
                              "#DEA0FDFF" = "Mesothelial",
                              "#AA0DFEFF" = "Interstitial-Macrophages",
                              "#F8A19FFF" = "Pericytes",
                              "#325A9BFF" = "Club-Cells",
                              "#C4451CFF" = "NK-Cells",
                              "#1C8356FF" = "Alveolar-Macrophages", 
                              "#85660DFF" = "Col14a1+-Fibroblasts",
                              "#B10DA1FF" = "Col13a1+-Fibroblasts"
                              #"#FBE426FF" = "Mixed Endothelia" #exclude Mixed endothelia because there is no mixed endothelia in results
    ))
  
  # filter graph by weights
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  raw_edges <-  as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>%
    arrange(V3)
  
  # Normalize the weight score to quantiles
  # quartiles_weight <- quantcut(edges$V3, 4)
  # levels(quartiles_weight) <- c(1:4)
  # edges$V3 <- as.integer(quartiles_weight)
  
  
  nodes <-  unique(c(edges$V1, edges$V2))
  
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  # diagram colors
  grid_col <- graph_module %>%
    dplyr::filter(Genes %in% nodes) %>% 
    dplyr::select(Genes, color) %>%
    mutate(color = as.character(color)) %>%
    deframe()
  
  #https://vuetifyjs.com/en/styles/colors/#material-colors
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15,0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  # we go back to the first track and customize sector labels
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      this_node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color) %>%
        as.character()
      
      circos.rect(
        xlim[1],
        0,
        xlim[2],
        1,
        col = this_node_text_color,
        border = NA
      )
      
      circos.text(
        mean(xlim),
        2,
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col= this_node_text_color
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]],
                c(0, 1),
                link[[2]],
                c(0, 1),
                col = col_fun(link[[3]]),
                border = NA)
  }
  
  
  # plot legend
  #lgd <- Legend(title ="Score", col_fun = col_fun)
  #grid.draw(lgd)
  #
  #png(
  #  paste(paste0(prefix, "_legend.png"), sep = ""),
  #  width = 1000,
  #  height = 1000,
  #  res = 300
  #)
  #grid.draw(lgd)
  #dev.off()
  
}


#----------------------run functions
id <- 0
df_list <- c() # the list to hold the CellChat object
# run CellChat for all datasets
net.N0 <- use_CellChat(counts_N0, meta_N0,"main.cluster.qi")
id <- id + 1
df_list[[id]] <- 'net.N0'
# convert the CCIs into the format to plot chord diagrams
dir<-"/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cellchat_plot"
for (i in 1:length(df_list)) {
  condition <- df_list[[i]] # the experiment condition
  prefix <- paste0(gsub(pattern = '^net.', replacement = '', df_list[[i]]))
  save_data(dir,get(df_list[i][[1]]), prefix)
}
# Plot figure
# mat.list <- Sys.glob("../data/codes_for_CellChat/all_cells/*_mat.csv") # the list of adjacency matrix files
# cat.list <- Sys.glob("../data/codes_for_CellChat/all_cells/*_cat.csv") # the list of the gene category files, i.e., cell types
mat.list <- Sys.glob("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cellchat_plot/*_mat.csv") # the list of adjacency matrix files
cat.list <- Sys.glob("/fs/ess/PCON0022/guoqi/Lung_sci_katherine/Results/Cellchat_plot/*_cat.csv") # the list of the gene category files, i.e., cell types
for (i in 1:length(mat.list)) {
  print(mat.list[i])
  a<-genChord(mat = mat.list[i], cat = cat.list[i],
           prefix = strsplit(x = mat.list[i], split = '_')[[1]][1])
}
# save plot a ...
