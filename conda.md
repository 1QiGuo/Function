# Install anaconda in windows
It need to configure environment like  R

## Check if conda is installed in this computer successfully  
```{r}
#run in cmd
conda
conda --v
```
## About anaconda
It's like a tool for create conda environment.  
There are two app in anaconda that are used frequently.  
1.Jupyter notebook: a online python for writing code  
2.Anaconda promoter: Like cmd, used for create conda environment

## error
“An HTTP error occurred when trying to retrieve this URL.” when create a new conda environment  

You need to follow this [guide](https://blog.csdn.net/z124560745/article/details/106819527)

- Note that using "%homepath%" to locate the home path.

#  A new conda environment  
## Create a new conda
```{r}
#run in cmd
#create
conda create -n yourenvironmentname python=X.X
#activate
conda init
#close shell and reopen
conda activate Testname

e.g.:
conda create -n spagft python==3.8
conda init
#reopen
conda activate spagft
```  
## Quit the new conda environment
```{r}
conda deactivate
```
## Delete a conda
```{r}
conda remove -n condaname all
```
## Look up all conda environment  
```{r}
conda env list
```

## Add conda to jupyter
```{r}
conda install ipykernel
python -m ipykernel install --user --name commot2 --display-name "commot2"

#check all environment in jupyter
jupyter kernelspec list
```
Please refer to [this website](https://blog.csdn.net/weixin_47381639/article/details/119798672)

# Git to connect github repositories
1. Download git [Reference here](https://blog.csdn.net/Passerby_Wang/article/details/120767020)
2. Configure environment [Reference here](https://www.cnblogs.com/ldq678/p/13287924.html)
3. Reopen cmd or anaconda promoter
```{r}
git clone githubwebsite/URL
e.g.:
git clone https://github.com/jxLiu-bio/SpaGFT
```
## error
“ERROR: Exception:
Traceback (most recent call last):”-------Turn off VPN
