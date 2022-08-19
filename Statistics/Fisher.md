# Fisher
## Descriotion
Whether two factors have relationship
1. Mutate or no mutate 2. Wild type or experiment type
1. Interesting genes 2. Gene modules
服从超几何检验

## Example
input: (character)-gene set1, gene set2, dataset1 containing all gene set1, dataset2 containing all gene set2. (union) 
output: (dataframe)-Intersect number of gene set1 and gene set2. odd_ratio(a correlation between group1 and group2)
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

better to show
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

# Ajusted p value:
adjust p value can only caculate a group of p value at a time and generate a group of adj_p value
```{r}
p.adj of p in a row<-p.adjust(result_fisher$P.value[which(result_fisher$Module==names(eight_modules)[i])],method = "fdr",ncol(result_fisher)
```
