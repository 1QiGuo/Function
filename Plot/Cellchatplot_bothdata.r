#reason for the bug, there is e.g.Tgfbr1_2
faith_final_data<-qread("./faith_final_data.qs")
meta_f<-faith_final_data@meta.data
split_faith <- SplitObject(faith_final_data, split.by = "orig.ident")
#preprocessing data from each condition
meta_list<-list()
count_list<-list()
for(i in 1:6){
  object <-split_faith[[i]]
  DefaultAssay(object) <- 'RNA'
  Idents(object)<- 'celltype'
  #count data
  counts_object <- GetAssayData(object, assay = "RNA", slot = "counts") # raw counts
  #normalize counts for cell chat
  counts_object <-NormalizeData(counts_object)
  meta_object <-object@meta.data
  meta_list[[i]]<-meta_object
  count_list[[i]]<-counts_object
}
id <- 0
df_list <- c() # the list to hold the CellChat object

net.1 <- use_CellChat(count_list[[1]], meta_list[[1]],"celltype")

#function
{ 
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
  mean_exp<-c()
  mean_exp <- sapply(1:nrow(df.net), function(x) {
    ct1 <- df.net$source[x] # source cell typeQ
    ct2 <- df.net$target[x] # target cell type
    cells1 <- rownames(meta[meta[[celltype_column]] == ct1,]) # cells belonging to the source cell type
    cells2 <- rownames(meta[meta[[celltype_column]] == ct2,]) # cells belonging to the target cell type
    #character "_" split;then;character including "-", do not change anything
    gene_pair1 <- df.net$ligand[x] # genes corresponding to the ligand
    if(grepl("_",gene_pair1)){
      #split
      gene_pair1<-strsplit(gene_pair1, split = '_')[[1]]
      #do not change anything if "-" existing
      genes1_unchange<-gene_pair1[grepl("-",gene_pair1)]
      genes1_lower<-capitalize(tolower(gene_pair1[!grepl("_",gene_pair1)]))
      genes1<-c(genes1_unchange,genes1_lower)
    }else{
      genes1=gene_pair1
    }
    gene_pair2 <- df.net$receptor[x] # genes corresponding to the receptor
    if(grepl("_",gene_pair2)){
      #split
      gene_pair2<-strsplit(gene_pair2, split = '_')[[1]]
      #do not change anything if "-" existing
      genes2_unchange<-gene_pair2[grepl("-",gene_pair2)]
      genes2_lower<-capitalize(tolower(gene_pair2[!grepl("-",gene_pair2)]))
      if (length(genes2_lower) > 1 & genes2_lower[2] == 'R2') {
        genes2_lower[2] <- gsub(pattern = '1', replacement = '2', genes2_lower[1])
      }
      genes2<-c(genes2_unchange,genes2_lower)
    }else{
      genes2=gene_pair2
    }
    cat("Ligand genes: ", genes1, "\n")
    cat("Receptor genes: ", genes2, "\n\n")
    #for receptor that not present but their interactor present (e.g.Tgfbr)
    if(any(!genes2 %in% rownames(counts))){#check whether not present, if not, perform the following code
      genes2=getreceptorfrominteractorname2(df.net,x)
    }
    #for ligand that not present but their interactor present
    if(any(!genes1 %in% rownames(counts))){#check whether not present, if not, perform the following code
      genes1=getligandfrominteractorname2(df.net,x)
    }
    mean1 <- mean(counts[genes1, cells1]) # the mean expression value of the ligand within the source cell type
    mean2 <- mean(counts[genes2, cells2]) # the mean expression value of the receptor within the source cell type
    avg <- (mean1 + mean2) / 2 # the mean expression value of the ligands and receptors
    return(avg)
  })
  df.net <- cbind(cell_type_pair, mean_exp, df.net) # add the mean expression and cell type pairs
  return(df.net)
}

#define function
getligandfrominteractorname2<-function(df.net,i){
  gene_pair<-df.net$interaction_name_2[i]
  gene_pair <- gsub(" ", "", gene_pair)
  gene_pair<-strsplit(gene_pair, split = '-')[[1]]
  #for ligand
  if(grepl("\\+",gene_pair[1])){
    gene_pair[1]<-gsub("[()]", "", gene_pair[1])
    gene_complex<-strsplit(gene_pair[1], split = '\\+')[[1]]
    genes1<-gene_complex
  }
  else{
    genes1=gene_pair[1]
  }
  return(genes1)
}
getreceptorfrominteractorname2<-function(df.net,i){
  gene_pair<-df.net$interaction_name_2[i]
  gene_pair <- gsub(" ", "", gene_pair)
  gene_pair<-strsplit(gene_pair, split = '-')[[1]]
  #for receptor
  if(grepl("\\+",gene_pair[2])){
    gene_pair[2]<-gsub("[()]", "", gene_pair[2])
    gene_complex<-strsplit(gene_pair[2], split = '\\+')[[1]]
    genes2<-gene_complex
  }
  else{
    genes2=gene_pair[2]
  }
  return(genes2)
}
}
