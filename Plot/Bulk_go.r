GOenrichment_mouse <- function(Seurat.DEGs = NULL) {
  library(org.Mm.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs
  #Obtain a GO object
  GO_bp <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  
  if (is.null(GO_bp)) {
    GO_simplied_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_GO"))
  } else{
    GO_simplied_res <-
      GO_bp@result[GO_bp@result$pvalue < 0.05, ]
    dim(GO_simplied_res)
    if (dim(GO_simplied_res)[1] == 0) {
      GO_simplied_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_GO"))
    }
  }
  return(GO_simplied_res)
}
GO_csv_plot<-function(deg_list){
  GO_up_list<-list()
  GO_down_list<-list()
  for(i in seq_along(deg_list)){
    sample_name <- names(deg_list)[i]
    up_genes <- rownames(deg_list[[i]][deg_list[[i]]$log2FoldChange > 0, ])
    GO_results_up<-GOenrichment_mouse(up_genes)
    GO_results_up_f<-GO_results_up[GO_results_up$p.adjust < 0.05,]
    GO_up_list[[i]]<-GO_results_up_f
    #up
    names(GO_up_list)[i]<-c(paste("BP_0.05_cluster",sample_name,"up",sep = "_"))
    write.csv(GO_results_up_f,paste("GO_BP_0.05",sample_name,"up.csv",sep = "_"))
    #down
    down_genes <- rownames(deg_list[[i]][deg_list[[i]]$log2FoldChange < 0, ])
    GO_results_down<-GOenrichment_mouse(down_genes)
    GO_results_down_f<-GO_results_down[GO_results_down$p.adjust < 0.05,]
    GO_down_list[[i]]<-GO_results_down_f
    names(GO_down_list)[i]<-c(paste("BP_0.05_cluster",sample_name,"_down",sep = "_"))
    write.csv(GO_results_down_f,paste("GO_BP_0.05_cluster",sample_name,"down.csv",sep = "_"))
    print(i)
  }
  #visualization
  for(i in seq_along(deg_list)){
    #downregulate
    data<-GO_down_list[[i]]
    if (is.data.frame(data) && nrow(data) > 0) {
    temp_sort<-data[with(data, order(Count,decreasing = T)), ]
    temp_sort$Description<-factor(temp_sort$Description,levels=(temp_sort$Description))
    count_sum<-as.numeric(unlist(strsplit(temp_sort$GeneRatio[1],split="/"))[2])
    temp_sort$GeneRatio_num<-c(temp_sort$Count)/count_sum
    library(ggplot2)
    p<-ggplot(temp_sort[1:10,], # you can replace the numbers to the row number of pathway of your interest
              aes(x = GeneRatio_num, y = Description)) + 
      geom_point(aes(size = Count, color = p.adjust)) +
      theme_bw(base_size = 20) +
      theme(axis.text = element_text(size = 10, face = "bold"),
      )+
      scale_colour_gradient(limits=NULL, low="red",high="blue") +
      ylab(NULL) +
      ggtitle("GO enrichment")
    #save
    ggsave(
      plot = p,
      filename = paste0(
        names(GO_down_list)[i],"_top10.tiff"
      ),
      device = "tiff",
      dpi = 150,
      width = 15,
      height = 10,
      units = "in")
    }
    #upregulate
    data<-GO_up_list[[i]]
    if (is.data.frame(data) && nrow(data) > 0) {
    temp_sort<-data[with(data, order(Count,decreasing = T)), ]
    temp_sort$Description<-factor(temp_sort$Description,levels=(temp_sort$Description))
    count_sum<-as.numeric(unlist(strsplit(temp_sort$GeneRatio[1],split="/"))[2])
    temp_sort$GeneRatio_num<-c(temp_sort$Count)/count_sum
    library(ggplot2)
    p<-ggplot(temp_sort[1:10,], # you can replace the numbers to the row number of pathway of your interest
              aes(x = GeneRatio_num, y = Description)) + 
      geom_point(aes(size = Count, color = p.adjust)) +
      theme_bw(base_size = 20) +
      theme(axis.text = element_text(size = 10, face = "bold"),
      )+
      scale_colour_gradient(limits=NULL, low="red",high="blue") +
      ylab(NULL) +
      ggtitle("GO enrichment")
    #save
    ggsave(
      plot = p,
      filename = paste0(
        names(GO_up_list)[i],"_top10.tiff"
      ),
      device = "tiff",
      dpi = 150,
      width = 15,
      height = 10,
      units = "in")
    }
  }
}

# treeplot
treeplot_go<-function(deg_df,cluster,direc,fontsize,shownumb){
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  deg_df_cluster<-deg_df[which(deg_df$cluster==cluster),]
  if(direc=="up"){
    deg_df_cluster_up<-deg_df_cluster[which(deg_df_cluster$avg_log2FC>0),]
    genes<-deg_df_cluster_up$gene
  }else{
    deg_df_cluster_down<-deg_df_cluster[which(deg_df_cluster$avg_log2FC<0),]
    genes<-deg_df_cluster_down$gene
  }
  GO_bp <-
    enrichGO(
      gene = genes,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
    )
  enrichres2 <- pairwise_termsim(GO_bp) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
  plot<-treeplot(enrichres2,fontsize = fontsize,showCategory=shownumb)
  ggsave(
      plot = plot,
      filename = paste(
        "tree_go",cluster,direc,".tiff",sep = "_"
      ),
      device = "tiff",
      dpi = 150,
      width = 18,
      height = 8,
      units = "in")
}

# treemap plot
library("GOSemSim")
MmGO <- godata('org.Mm.eg.db', ont="BP")
cluster0_up_go<-read.csv("./GO_BP_0.05_cluster_0_up.csv")
cluster1_up_go<-read.csv("./GO_BP_0.05_cluster_1_up.csv")
cluster9_up_go<-read.csv("./GO_BP_0.05_cluster_9_up.csv")
cluster_0_1<-mgoSim(cluster0_up_go$ID[1:5], cluster1_up_go$ID[1:5], semData=MmGO, measure="Wang", combine="BMA")
cluster_0_9<-mgoSim(cluster0_up_go$ID, cluster9_up_go$ID, semData=MmGO, measure="Wang", combine="BMA")
cluster_1_9<-mgoSim(cluster1_up_go$ID, cluster9_up_go$ID, semData=MmGO, measure="Wang", combine="BMA")
cluster_up<-list()
cluster_down<-list()
for(i in 0:11){
  cluster_up[[i+1]]<-read.csv(paste0("./GO_BP_0.05_cluster_",i,"_up.csv"))  
  cluster_down[[i+1]]<-read.csv(paste0("./GO_BP_0.05_cluster_",i,"_down.csv"))
  names(cluster_up)[i+1]<-paste0("cluster_",i,"_up")
  names(cluster_down)[i+1]<-paste0("cluster_",i,"_down")
}
cluster_up_sim<-data.frame(matrix(0,nrow = 12,ncol = 12))
rownames(cluster_up_sim)<-0:11
colnames(cluster_up_sim)<-0:11

for(i in 0:11){
  a<-cluster_up[[i+1]]
  for(j in i:11){
  b<-cluster_up[[j+1]]
  cluster_sim<-mgoSim(a$ID, b$ID, semData=MmGO, measure="Wang", combine="BMA")
  cluster_up_sim[i+1,j+1]<-cluster_sim
  }
}

BiocManager::install("rrvgo")
install.packages("tm")
library("tm")
library("rrvgo")

setwd("/Users/guo.2101/BMBL Dropbox/BMBL Project/Qi-Guo_Projects/Collaborations-Popovich/Dr. Popovich_Dr. Ma_Qi/1.Spinalcord_atlas_Qi/3-Annotation/Healthy_microglia_subcluster_03162025/Reduced_GO_vis/")
for(i in 1:length(cluster_up)){
 outputname_up=names(cluster_up)[i]
 reduceGO_vis(cluster_up[[i]],outputname_up)
 outputname_down=names(cluster_down)[i]
 reduceGO_vis(cluster_down[[i]],outputname_down)
}
reduceGO_vis<-function(go_df,output_name){
  simMatrix <- calculateSimMatrix(go_df$ID,
                                  orgdb="org.Mm.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames(go_df$Count, go_df$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Mm.eg.db")
  pdf(paste0(output_name,".pdf"))
  treemapPlot(reducedTerms)
  dev.off()
}  
