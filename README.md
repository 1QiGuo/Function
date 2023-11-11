# Function
## Enrichment analysis with three function

input: a dataframe which must has gene names as rownames
output: a dataframe consisting of GO, KEGG, reactome analysis.
```{r}
RunPathway_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  genes.use <- rownames(Seurat.DEGs)
  #Obtain a GO object
  GO_pathway <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  if (is.null(GO_pathway)) {
    GO_simplied_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_GO"))
  } else{
    #GO_res <- simplify(GO_pathway)
    GO_simplied_res <-
      GO_pathway@result[GO_pathway@result$pvalue < 0.05, ]
    dim(GO_simplied_res)
    if (dim(GO_simplied_res)[1] == 0) {
      GO_simplied_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_GO"))
    }
  }
  #convert gene type for pathway annotation method-enrichPathway based on reactome and kegg
  gene.convert <-
    bitr(
      genes.use,
      fromType = "SYMBOL",
      toType = c("ENSEMBL", "ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
  reactome_analysis <-
    enrichPathway(
      gene = gene.convert$ENTREZID,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  if (is.null(reactome_analysis)) {
    reactome_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_reactome"))
  } else{
    reactome_res <-
      reactome_analysis@result[reactome_analysis@result$pvalue < 0.05, ]
    dim(reactome_res)
    if (dim(reactome_res)[1] == 0) {
      reactome_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_reactome"))
    }
  }
  
  KEGG_pathway <-
    enrichKEGG(
      gene = gene.convert$ENTREZID,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  #KEGG_pathway@result[1:5,]
  if (is.null(KEGG_pathway)) {
    KEGG_res <- as.data.frame(cbind(NO_results = "NO_result_for_KEGG"))
  } else{
    KEGG_res <- KEGG_pathway@result[KEGG_pathway@result$pvalue < 0.05, ]
    dim(KEGG_res)
    if (dim(KEGG_res)[1] == 0) {
      KEGG_res <- as.data.frame(cbind(NO_results = "NO_result_for_KEGG"))
    } else{
      for (j in 1:nrow(KEGG_res)) {
        tmp_name <- KEGG_res$geneID[j]
        tmp_name <- unlist(strsplit(tmp_name, "/"))
        index_entizID <- gene.convert$ENTREZID %in% tmp_name
        converted_geneSymbol <- gene.convert$SYMBOL[index_entizID]
        final_name <- paste0(converted_geneSymbol, collapse = "/")
        KEGG_res$geneID[j] <- final_name
      }
    }
  }
  Go_KEGG.bind <-
    gdata::cbindX(GO_simplied_res, KEGG_res, reactome_res)
  final.data <- gdata::cbindX(Seurat.DEGs, Go_KEGG.bind)
  return(final.data)
}
```

GO enrichment analysis
Note: input is character format
```{r}
GOenrichment_human <- function(Seurat.DEGs = NULL) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs
  #Obtain a GO object
  GO_pathway <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  
  if (is.null(GO_pathway)) {
    GO_simplied_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_GO"))
  } else{
    #GO_res <- simplify(GO_pathway)
    GO_simplied_res <-
      GO_pathway@result[GO_pathway@result$pvalue < 0.05, ]
    dim(GO_simplied_res)
    if (dim(GO_simplied_res)[1] == 0) {
      GO_simplied_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_GO"))
    }
  }
  return(GO_simplied_res)
}
```

showing specific categroies
input is a dataframe which consist of selected categroies
```{r}
library(multienrichjam)
multi<-enrichDF2enrichResult(module4[["Pink"]])#go object
#plot
p<-dotplot(multi,font.size=20)
#Note: if the plot have decimals count number(range of reference, do not stand for any real meaning) and you don't want them.
#solve count decimals problem
min = min(target_data$Count)
max = max(target_data$Count)
# 10 is length
step = ceiling((max-min)/5)
p <- p + scale_size_continuous(breaks=seq(min, max, by=step))
p<-p+ theme(legend.key.size = unit(1.0, 'cm'),
           legend.title = element_text(size=18), #change legend title font size
           legend.text = element_text(size=18))
```
## Fisher test
input: (character)-gene set1, gene set2, dataset1 containing all gene set1, dataset2 containing all gene set2. (union)
output: (dataframe)-Intersect number of gene set1 and gene set2. odd_ratio(a correlation between group1 and group2)
### generate a more comprehensive consequence
```{r}
Fishertest<-function(gene1,gene2,dataset1,dataset2){
  a<-length(intersect(gene1,gene2))
  b<-length(intersect(setdiff(dataset2,gene2),gene1))
  c<-length(intersect(setdiff(dataset1,gene1),gene2))
  d<-length(intersect(setdiff(dataset1,gene1),setdiff(dataset2,gene2)))
  data_fisher<-data.frame(gene2=c(a,c),dataset2=c(b,d))
  rownames(data_fisher)<-c("gene1","dataset1")
  re<-fisher.test(data_fisher)
  expect<-length(intersect(dataset1,dataset2))*(length(gene1)/length(dataset1))*(length(gene2)/length(dataset2))
  result<-data.frame(Intersection_num=a,Expection=expect,Odd_ratio=re[["estimate"]][["odds ratio"]],P.value=re$p.value)
  return(result)
}
```

### better to show
```{r}
Fishertest_show<-function(gene1,gene2,dataset1,dataset2){
  a<-length(intersect(gene1,gene2))
  b<-length(intersect(setdiff(dataset2,gene2),gene1))
  c<-length(intersect(setdiff(dataset1,gene1),gene2))
  d<-length(intersect(setdiff(dataset1,gene1),setdiff(dataset2,gene2)))
  data_fisher<-data.frame(gene2=c(a,c),dataset2=c(b,d))
  rownames(data_fisher)<-c("gene1","dataset1")
  re<-fisher.test(data_fisher)
  expect<-length(intersect(dataset1,dataset2))*(length(gene1)/length(dataset1))*(length(gene2)/length(dataset2))
  result<-c(round(a,1),round(re[["estimate"]][["odds ratio"]],4),re$p.value)
  return(result)
}
```
## adjust p value
adjust p value can only caculate a group of p value at a time and generate a group of adj_p value
```{r}
p.adj of p in a row<-p.adjust(result_fisher$P.value[which(result_fisher$Module==names(eight_modules)[i])],method = "fdr",ncol(result_fisher)
```

## WGCNA
library("WGCNA")
cutoffs need to choose in person-cex1,sft$powerEstimate
```{r}
#obtain highly variable genes
H_2_8 <- HVG_2_8@assays$SCT@var.features
H_2_3 <- HVG_2_3@assays$SCT@var.features
H_1_1 <- HVG_1_1@assays$SCT@var.features
H_T4857 <- HVG_T4857@assays$SCT@var.features
H_18_64 <- HVG_18_64@assays$SCT@var.features
H_2_5 <- HVG_2_5@assays$SCT@var.features
All_HVG <- unique(c(H_2_5, H_18_64, H_T4857, H_1_1, H_2_3, H_2_8))

#obtain data matrix(sct-continuous data)
all_matrix_hvg <- sample_6@assays$SCT@data
all_matrix_hvg <-
  sample_6@assays$SCT@data[match(All_HVG, rownames(all_matrix_hvg)), ]
#row:spot,col:gene.
all_matrix_hvg <- t(as.matrix(all_matrix_hvg))
dim(all_matrix_hvg)
gene.names = colnames(all_matrix_hvg)
# caculate soft threshold-β，
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))  # in practice this should include powers up to 20.
sft = pickSoftThreshold(
  all_matrix_hvg,
  dataIsExpr = TRUE,
  powerVector = powers,
  corFnc = cor,
  corOptions = list(use = 'p'),
  networkType = "unsigned"
)
sft$powerEstimate
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
# SFT index as a function of different powers
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.9, col = "red")
# Mean connectivity as a function of different powers
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
softPower = 2
# #automatic module detection
# mergingThresh = 0
bwnet = blockwiseModules(
  all_matrix_hvg,
  corType = "spearman",
  maxBlockSize = 5000,
  mergeCutHeight = 0.25,
  networkType = "unsigned",
  power = softPower,
  minModuleSize = 30,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMFileBase = "maleMouseTOM-blockwise",
  verbose = 3
)
TOM = TOMsimilarityFromExpr(
  all_matrix_hvg,
  networkType = "unsigned",
  TOMType = "unsigned",
  power = softPower
)
colnames(TOM) = rownames(TOM) = gene.names

dissTOM = 1 - TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")
plot(geneTree,
     xlab = "",
     sub = "",
     cex = 0.3)


# Set the minimum module size
minModuleSize = 20
# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(
  dendro = geneTree,
  method = "tree",
  minClusterSize = minModuleSize,
  deepSplit = 3
)
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

#discard the unassigned genes, and focus on the rest
restGenes = (dynamicColors != "grey")
diss1 = 1 - TOMsimilarityFromExpr(all_matrix_hvg[, restGenes], power = softPower)

colnames(diss1) = rownames(diss1) = gene.names[restGenes]
hier1 = flashClust(as.dist(diss1), method = "average")
plotDendroAndColors(
  hier1,
  dynamicColors[restGenes],
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
#set the diagonal of the dissimilarity to NA
diag(diss1) = NA
#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7, 7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
# extract module
module_colors = setdiff(unique(dynamicColors), "grey")
for (color in module_colors) {
  module = gene.names[which(dynamicColors == color)]
  write.table(
    module,
    paste("WGCNA/", "module_", color, ".txt", sep = ""),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}
```
## closeness centrality for modules
input:count matrix, a dataframe of two column. First is gene name and second is module name. 
```{r}
library(igraph)
module.cat<-data.frame(id=eight_modules_id,genes=eight_modules_genes)
counts<-sample_6@assays$SCT@counts

graph_module <- module.cat
gene_from_8Module <- c()

for (i in 1:8){
genes <- graph_module$genes[graph_module$id ==colors[i]]
all_matrix_DEG <- counts[genes,]
dim(all_matrix_DEG)
all_matrix_DEG <- t(as.matrix(all_matrix_DEG))
similarity.mat <- TOMsimilarityFromExpr(all_matrix_DEG, power = 3)
graph_adj <- adjacency.fromSimilarity( similarity.mat,power = 3,type = "unsigned")
colnames(graph_adj) <- genes
rownames(graph_adj) <- genes
graph_igraph <- graph_from_adjacency_matrix(graph_adj,mode = "undirected", weighted = T)
closebess_centrality <- closeness(graph_igraph)
#data.frame
sorted_centrality <- as.data.frame(cbind(Gene = names(closebess_centrality), 
                                         cloness.centrality = closebess_centrality))
sorted_centrality <- sorted_centrality[order(as.numeric(sorted_centrality$cloness.centrality),decreasing = T),]
write.csv(sorted_centrality,  file = paste0("/users/PAS1475/liuzeyi/guoqi/output/WGCNA/closeness centrality/",colors[i],".csv"),row.names = F)
closebess_centrality_gene_names <- names(sort(closebess_centrality,decreasing = T)[1:15])
#15
gene_from_8Module <- c(gene_from_8Module,closebess_centrality_gene_names)
}
```

# match
```{r}
match(x,y)
#return the location of x in y
```
