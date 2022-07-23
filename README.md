# STAE
STAE:Spatial Temporal Auto-Encoder
website url:https://cell.ownbox.cn/

## Overview
## Installation
### Setup
STAE is available from GitHub with:

```
#If you don't have devtools installed, please install it first

install.packages("devtools")

devtools::install_github("zenglab-regeneration/STAE")

```

### Depend

Some programs of our project require a python environment, so if you don't have a python environment, please follow the steps below.  
* install [anaconda](https://www.anaconda.com/ "anaconda")
* Create a python3.9 environment  

When you have anaconda installed, you need to create a python conda.
```
conda create -n testconda python = 3.9
```

## Examples from paper
### Dataset 
- Single-cell transcriptome gene expression matrix
- Spatial transcriptome gene expression matrix
- Single cell type data

### Environment settings


```
library('STAE')
help(package = 'DRIM')

#set the py conda
env_python_set("D:/anaconda/envs/testconda")

#Check the dependent environment for the program to run, and automatically install the missing python package
env_test()
```

