# Hypergeometric_test
## Description
超几何分布检验常用来对venn图两个圈overlap的显著性进行检验
## Example
设总共有29个人，其中11个吸烟者，18个非吸烟者，现从中随机抽取16个样本（在此实验中对应着肺癌病人），有10个是吸烟者，这样的事件是否显著？
```{r}
p-value=phyper(10-1, 11, 18, 16, lower.tail=F)=0.003135274
```
