# Pipiline

```{r}
a<-micro_obj %>% NormalizeData() %>% FindVariableFeatures()%>% ScaleData() %>% RunPCA (verbose=F)%>% FindNeighbors() %>% FindClusters()%>% RunUMAP(reduction = "pca", dims = 1:50) 
b<-a%>%DimPlot()
```
# Options

[example](https://blog.csdn.net/weixin_34233679/article/details/86265275)
```{r}
#adjust digit that shows in R
options(digits=10) #show 10
```
# shortcut key

Ctrl + I – Fixes line indentations. Ctrl + Shift + A – Does a complete reformat of the selected part of a code.

# apply
[example](https://www.jianshu.com/p/82a94d5bbbad)

# convert from seurat object to h5ad
```
v1
use_condaenv('/home/qiguo/.conda/envs/qi/envs/scvi-env/')
health_microglia<-qread("/bmbl_data/qiguo/SCI/atlas/integration/processeddata/health_integration_microglia_4qc_0411.qs")
health_microglia_v3 <- Convert_Assay(seurat_object = health_microglia, convert_to = "V3",assay = "RNA")
sceasy::convertFormat(health_microglia_v3, from="seurat", to="anndata",
                      outFile='/bmbl_data/qiguo/SCI/atlas/integration/processeddata/health_integration_microglia_4qc_0411.h5ad')
#v2
convert_purev5object_to_h5ad<-function(v5object,outputfile){
  options(Seurat.object.assay.version = "v3")
  count<-v5object@assays$RNA@layers$counts
  colnames(count)<-colnames(v5object)
  rownames(count)<-rownames(v5object)
  v3object <- CreateSeuratObject(
    count,
    assay = "RNA",
    min.cells = 0,
    min.features = 0
  )
  print(identical(colnames(v3object),colnames(v5object)))
  v3object@meta.data<-v5object@meta.data
  use_condaenv('/home/qiguo/.conda/envs/qi/envs/scvi-env/')
  #healthy_merge_4qc_v3_python <- Convert_Assay(seurat_object = healthy_merge_4qc, convert_to = "V3")
  sceasy::convertFormat(v3object, from="seurat", to="anndata",
                        outFile=outputfile)
  options(Seurat.object.assay.version = "v5")
}
```
