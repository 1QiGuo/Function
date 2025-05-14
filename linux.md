# Change permissions
## change all document in this file
chmod /fs/project/PAS1475/guoqi/spatial_data/ -R 777  
chmod -R 777 *  
**Note: do not use blank when name folder**  

# compress and decompress file 

## compress

```{r}
gzip FileName
zip -r my_folder.zip my_directory_to_zip
```

## decompress file
end with zip
```{r}
unzip FileName.zip
```
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
qstat -u guoqi
https://maveric-informatics.readthedocs.io/en/latest/OSC.html
```
# transfer file by remote server
```{r}
scp -r guoqi@192.148.247.179:/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/temp ./rawdata/temp
```
# Run sh script
```
/bin/bash .command.sh
```

# Find specific file 
```
specific file in this path
cd path
find ./ -iname projection*
#find specific file in all path
find . -type f -name "*h5Seurat*"
```
# Check object storage
```
du -ah
```

# Show all files including hidden
```
ls -a
```

# OSC

## create sh command
```{r}
#!/bin/bash
#SBATCH --job-name=directnet_examp
#SBATCH --time=20:50:59
#SBATCH --output="directnet_sh_out"
#SBATCH --account=PCON0022
#SBATCH --mem=150GB
#SBATCH --mail-type=BEGIN,END,FAIL

cd /fs/ess/PCON0022/guoqi/Yang/Stream/Directnet
start=`date +%s`
module load R/4.1.0-gnu9.1
Rscript /fs/ess/PCON0022/guoqi/Yang/Stream/Directnet/Directnet.r
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
```

# Change file name in batch
```
for file in *_expression*; do     mv -- "$file" "${file//_expression}"; done
for file in sim_*; do     new_name=$(echo "$file" | sed 's/_/-/g');     mv "$file" "$new_name"; done
for file in *.h5ad*; do
    newname="${file/.h5ad/_sim_sim_expression.h5ad}"
    mv -- "$file" "$newname"
done
for file in *.h5Seurat*; do
    newname="${file/.h5Seurat/_sim_sim_expression.h5Seurat}"
    mv -- "$file" "$newname"
done

for file in *_xy*; do     mv -- "$file" "${file//_xy}"; done
for file in sim_*; do     new_name=$(echo "$file" | sed 's/_/-/g');     mv "$file" "$new_name"; done
for file in *.csv*; do
    newname="${file/.csv/_sim_sim_xy.csv}"
    mv -- "$file" "$newname"
done
```

# VP03

## Check id and cancel

```
ps aux | grep 04_Integration_scvi.py
#id is the first number in the first row : 633763
kill id
```


## Enter a conda environment

```{r}
conda activate ...
```


## Submit the r script

```{r}
nohup Rscript --no-save 05_integrateallct.r > output_05allctintegrate_r.log 2>&1 &
```

## Submit the Python script

```{r}
nohup python ML_microglia_vali_1.2.py > output_1.2.log 2>&1 &
```

## Jobs

```{r}
jobs
jobs -ll
kills ...

top
```

## check the save

```{r}

du -sh */
```

## make a fake folder

```{r}
ln -s 'new folder' 'old folder'
ln -s /bmbl_data/qiguo/conda_envs/scvi-env /home/qiguo/.conda/envs/qi/envs/
```

# check GPU
```{r}
nvidia-smi
```
