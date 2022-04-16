# Read xlsx
```{r}
library("readxl")
up_deg <- read_excel("No.4, AD vs nonAD GO for dot plot.xlsx", sheet = 1)
```

Second method- more direct
```{r}
library(openxlsx)
AD_DEGs_up <- read.xlsx("/users/PAS1475/liuzeyi/guoqi/output/picture/heatmap/Fig4 heatmap genelist.xlsx",sheet = 1)
```
# Write xlsx
```{r}
dyn.load(
  "/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.41.x86_64/jre/lib/amd64/server/libjvm.so"
)
library(xlsx, lib.loc = "/fs/scratch/PCON0022/liyang/lib")
#define first sheet and for function for other sheet
 write.xlsx2(list_all[[i]][[1]], paste0(names(list_all)[i],".xlsx"), sheetName = names(list_all[[1]])[1],
              col.names = TRUE, row.names = T, append = FALSE)
  for(j in 2:3){
    write.xlsx2(list_all[[i]][[j]], paste0(names(list_all)[i],".xlsx"), sheetName = names(list_all[[1]])[j],
                col.names = TRUE, row.names = T, append = TRUE)
  }
```

# Names of sheet in xlsx
```{r}
specific<-excel_sheets(path = "DEGs of WFS1AT8 v.s. level 123 in layer 2-specific-separate level.xlsx")
```
# qs
```{r}
###qs
library(qs)
#read
raw<-qread(file = 'F:/20180401.qs')
#save
qsave(raw, file = 'F:/20180401.qs')
#save many variables in R environment
library(rlang)
qsavem(LEV,offroad,p,plotdata,raw,result,file = 'F:/myopjects.qs')
```
