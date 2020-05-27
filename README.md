# scGRN

A computational pipeline for predicting cell-type gene regulatory networks via multi-omics data integration and analyses. 

## Introduction
Understanding the gene regulatory mechanisms at the cellular resolution remains challenging. To address this, we have developed a computational pipeline to integrate single cell multi-omic data including chromatin interactions, epigenomics and transcriptomics, for predicting cell-type gene regulatory networks linking transcription factors, distal regulatory elements and target genes in the human brain (e.g., excitatory and inhibitory neurons, microglia, oligodendrocyte). The pipeline is depicted in the below flow chart.

![Pipeline](https://github.com/daifengwanglab/scGRN/blob/master/pipeline.png)


## Prerequisite
Our package needs the following packages to work.
> GenomicInteractions,biomaRt,GenomicRanges,GenomicFeatures,GenomeInfoDb,IRanges,S4Vectors,BSgenome.Hsapiens.UCSC.hg19,TxDb.Hsapiens.UCSC.hg19.knownGene,TFBSTools,JASPAR2018,glmnet,motifmatchr,
data.table,dplyr,parallel,doParallel,foreach

All these

## Install scGRN R package
The devtools package provides install_github() that enables people to install packages from GitHub. You can use the following code to download the scGRN package.

```R
install.packages("devtools")
library(devtools)
install_github("daifengwanglab/scGRN")
```

After the installation, you can directly use **library(scGRN)** to load the package.

## Functions in the package

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
    
## Example

To present how our pipeline works, we will use microglia interactome data [1], microglia enhancers [1], gene expression data [2] to illustrate the steps.

To begin with, we prepare and format the data in R .

```{r}
library(readxl)
interactome_data = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia interactome',skip = 2)[,1:6]

enhancers = read_xlsx("PLAC-seq promoter interactome map.xlsx",
                       sheet = 'Microglia enhancers',skip = 2,
                       col_names = c('chr','start','end'))
    
gexpr = read.table('DER-22_Single_cell_expression_raw_UMI.tsv',
                       header = T,row.names = 1,sep = '\t')
# Here we do some simple data processing without normalizing the data
# Select microglia cells 
# Remove genes that are not expressed (based on the gene expression data)
gexpr = gexpr[,grep('Micro',colnames(gexpr))]
gexpr = gexpr[rowSums(gexpr != 0) > 5,]
```

Now we will start using the functions in our package.

**Step1: find interaction**

```{r}
df1 <- scGRN_interaction(interactome_data,enhancers)
# For convenience, I use a subset of df1 and show relevant outputs.
df1 <- df1[sample(nrow(df1),20),]
head(df1)
```

    ##         gene gene_chr promoter_start promoter_end enh_chr enh_start   enh_end
    ## 70681  PICK1    chr22       38451253     38456253   chr22  38731974  38732426
    ## 59812  CBLL1     chr7      107381779    107386779    chr7 107515988 107516634
    ## 25923  GSK3A    chr19       42743546     42748546   chr19  42629613  42633317
    ## 4621    CHD3    chr17        7810319      7815319   chr17   7942326   7946725
    ## 62217 PICALM    chr11       85689772     85694772   chr11  85859249  85859502
    ## 28046  RAB43     chr3      128838148    128843148    chr3 128952814 128954805



**Step2: get TFs for each promoter and enhancer**

```{r}
df2 <- scGRN_getTF(df1)
head(df2,1)
```
   
    ##     gene                promoter                enhancer
    ## 1: PICK1 chr22:38451253-38456253 chr22:38731974-38732426
    ##                             promoter_TF                          enhancer_TF
    ## 1: IRF2,PPARG,RXRA,RREB1,NR1H2,REST,... CREB1,MZF1,PPARG,RXRA,RORA,NR1H2,...
 
**Step3: predict TF -> TG via Elastic net**

Last, what we need to do is to input the gene expression data and the output of function *scGRN_getTF*.

```{r}
df3 <- scGRN_getNt(df = df2, gexpr = gexpr)
head(df3)
```
    ##      TG    TF                enhancer                promoter     TFbs          coef
    ## 1 PICK1 CREB1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer   0.004038054
    ## 2 PICK1  MZF1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer   0.016943537
    ## 3 PICK1 PPARG chr22:38731974-38732426 chr22:38451253-38456253     both  -0.008975126
    ## 4 PICK1  RXRA chr22:38731974-38732426 chr22:38451253-38456253     both  -0.002222325
    ## 6 PICK1 NR1H2 chr22:38731974-38732426 chr22:38451253-38456253     both   0.002630548
    ## 7 PICK1   SP1 chr22:38731974-38732426 chr22:38451253-38456253 enhancer  -0.002453652


## Reference and source
1. Nott A, Holtman IR, Coufal NG, Schlachetzki JCM, Yu M, Hu R, et al. Brain cell type-specific enhancer-promoter interactome maps and disease-risk association. Science. 2019;366:1134–9. 

2. Wang D, Liu S, Warrell J, Won H, Shi X, Navarro FCP, et al. Comprehensive functional genomic resource and integrative model for the human brain. Science. 2018;362. 

