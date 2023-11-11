# Convert gene symbol to ensemble (mouse)

The function getBM is the basic function that is used to convert gene names in a species.
```{R}
library(biomaRt)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library("org.Mm.eg.db")
ensembl_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

z <- getBM(c("ensembl_gene_id","mgi_symbol"), values = gene_hugo, mart = ensembl_mart)
#create gene list
gene_list<-list()
for(i in 1:length(health_list_1108)){
  gene_hugo<-rownames(health_list_1108[[i]])
  z <- getBM(c("ensembl_gene_id","mgi_symbol"), filters = "mgi_symbol", values = gene_hugo,mart = ensembl_mart)
  gene_list[[i]]<-z$ensembl_gene_id
}
intersect_gene<-Reduce(intersect, gene_list)
```

#  converting gene names of mouse and human

input: (character)gene names
output: (dataframe)mouse names and human names
The function getLDS is the basic function that is used to convert gene names in across species.
```{r}
library(biomaRt)
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
geneMn <- rownames(mtx)
geneHs <- getLDS(attributes = "mgi_symbol",filters="mgi_symbol",values= geneMn,mart=mouse,attributesL = "hgnc_symbol",martL=human,uniqueRows = TRUE)
```
