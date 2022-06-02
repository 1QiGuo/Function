# Change permissions
## change all document in this file
chmod /fs/project/PAS1475/guoqi/spatial_data/ -R 777  
chmod -R 777 *  
**Note: do not use blank when name folder**  

## compress and decompress file 
### compress
```{r}
gzip FileName
```

### decompress file
end with gz--gzip
```{r}
gzip -d FileName.gz
```
end with tar.gz-tar
```{r}
tar zxvf FileName.tar.gz
```
# Copy file
```{r}
cd originalpath
cp file targetpath
```

# R  
## start R in terminal
```{r}
module load R/4.1.0-gnu9.1
R
```
## quit R in terminal
```{r}
quit()
```
# load JOB on server
```{linux}
#upload
sbatch xxx.sh
#supervise
qstat jobname
```
