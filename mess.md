# Apply

## Input is a list
### sapply
output is a vector
```{r}
a <- c("3_k","4_f","2_t","8_s")
b <- strsplit(a, "_")
c <- sapply(b, "[",1)
c
```
### lapply
output is a list
```{r}
d <- lapply(b, "[",1)
print(d)
```
