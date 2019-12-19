# scGRN

A R package designed for users to get gene regulatory network via Elastic net. 

## Introduction
Hi-C data is often used to analyze genome-wide chromatin organization like topologically associating domains, linearly contiguous regions of the genome that are associated in 3-D space. Thus we implement this pipeline which makes use of the hic data and user-provided enhancer list to discover latent gene regulatory network. The framework is depicted in the below figure.

![Pipeline](https://github.com/mying4/scGRN/blob/master/pipeline.png)


## Install scGRN R package
The devtools package provides install_github() that enables people to install packages from GitHub. You can use the following code to download the scGRN package.
```R
install.packages("devtools")
library(devtools)
install_github("daifengwanglab/tools/scGRN")
```

After the installation, you can directly use **library(scGRN)** to load the package.

### Functions in the package

There are three functions inside the package.

1. **scGRN_interaction** :
    * input:
        * hic_interaction: data frame containing the variables chr1,start1,end1,chr2,start2,end2
          or chr,start1,end1,start2,end2 
          chr1,ch2 or chr should have the following format 'chrD'(D represents digit 1-22 or X Y)                          
          start1,end1,start2,end2 should be integer type
          
        * enhancer: a data frame containing chr,start,end for enhancers
        
    * output: a data frame containing gene, gene_chr, promoter_start, promoter_end,               
      enh_chr,enh_start,enh_end
 
2. **scGRN_getTF** :
    * input: 
        * df: output data frame from **scGRN_interaction**
    
    * output: a data.table which contains transcription factors for each region

3. **scGRN_getNt** :
    * input: 
        * df: output data.table from **scGRN_getTF**
        * gexpr: gene expression data in which each row represents a gene and each column represents an observation
        
    * output:  a data frame containing TG, TF, promoter, enhancer and coef
    
### Example

Here we use hic data [Promoter-anchored chromatin loops](http://resource.psychencode.org/Datasets/Integrative/Promoter-anchored_chromatin_loops.bed) and enhancer list [DER-03a_hg19_PEC_enhancers](http://resource.psychencode.org/Datasets/Derived/DER-03a_hg19_PEC_enhancers.bed) as input and show how to use these functions.
For convenience, the file names are changed to hic_data.txt and enhancers.txt
