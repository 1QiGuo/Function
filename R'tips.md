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
