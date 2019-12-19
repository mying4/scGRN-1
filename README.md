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

Here we use hic data [Promoter-anchored chromatin loops](http://resource.psychencode.org/Datasets/Integrative/Promoter-anchored_chromatin_loops.bed), enhancer list [DER-03a_hg19_PEC_enhancers](http://resource.psychencode.org/Datasets/Derived/DER-03a_hg19_PEC_enhancers.bed)  and gene expression data [DER-02_PEC_Gene_expression_matrix_TPM](http://resource.psychencode.org/Datasets/Derived/DER-02_PEC_Gene_expression_matrix_TPM.txt).
For convenience, the file names are changed to hic_data.txt, enhancers.txt and gexpr.txt respectively.

```{r}
# format the column names
# change the column names  as standard
# change the chromosome format as standard
hic_data <- read.table('hic_data.txt',skip = 1, header = F, col.names = c('chr','start1','end1','start2','end2')) 
hic_data$chr <- paste('chr',hic_data$chr,sep = '')
enhancers <- read.table('enhancers.txt',header = T,stringsAsFactors = F)[,1:3]
colnames(enhancers) <- c('chr','start','end')
# it is a good option to use fread in data.table to speep up the following process
gexpr = read.table('gexpr.txt',header = T,row.names = 1)
```

Now use the functions in the package.
```{r}
df1 <- scGRN_interaction(hic_data,enhancers)
head(df1)
```
``` r
setwd('~/Desktop/scGRN/')
load('df1.Rda')
head(df1)
```

    ##          gene gene_chr promoter_start promoter_end enh_chr enh_start   enh_end
    ## 19475   NLRC3    chr16        3624893      3629893   chr16   3907323   3907699
    ## 5942    USP19     chr3       49155872     49160872    chr3  49434888  49435293
    ## 54258    ODF2     chr9      131255635    131260635    chr9 131320645 131321147
    ## 6275   FERMT2    chr14       53415316     53420316   chr14  53772807  53773548
    ## 26539 ZCCHC14    chr16       87454988     87459988   chr16  87573686  87575832
    ## 82796  BNIP3L     chr8       26245406     26250406    chr8  26484874  26485106
