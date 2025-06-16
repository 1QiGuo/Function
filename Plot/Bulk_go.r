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
