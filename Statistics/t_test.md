Do not use this if all assumptions have been satisfied, otherwise, wilcoxon.

# Samples are normal distribution

# Means are equal

# Sqares are equal

```{r}
o=lm(temp.df[1:6,i] ~ as.factor(c(rep(1,3), rep(2,3))))
plot(o)
```
