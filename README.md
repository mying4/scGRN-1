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
install_github("daifengwanglab/scGRN")
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
hic_data <- read.table('hic_data.txt',skip = 1, header = F, stringsAsFactors = F,col.names = c('chr','start1','end1','start2','end2')) 
hic_data$chr <- paste('chr',hic_data$chr,sep = '')
enhancers <- read.table('enhancers.txt',header = T,stringsAsFactors = F)[,1:3]
colnames(enhancers) <- c('chr','start','end')
# it is a good option to use fread in data.table to speep up the following process
gexpr = read.table('gexpr.txt',header = T,row.names = 1)
```

Now use the functions in the package.

**Step1: find interaction**
```{r}
df1 <- scGRN_interaction(hic_data,enhancers)
# Here I only use subset of df1 to demonstrate the steps.
df1 <- df1[sample(nrow(df1),10),]
head(df1)
```

    ##          gene gene_chr promoter_start promoter_end enh_chr enh_start   enh_end
    ## 18741   DUSP8    chr11        1590651      1595651   chr11   1579427   1581277
    ## 44189    POLN     chr4        2241361      2246361    chr4   2071734   2072014
    ## 88221 LDLRAD4    chr18       13639257     13644257   chr18  13560988  13561222
    ## 91073   GSDMD     chr8      144637583    144642583    chr8 144987814 144988401
    ## 5995    RAB10     chr2       26254229     26259229    chr2  25452258  25452926
    ## 73546 PITPNM2    chr12      123587547    123592547   chr12 123569147 123569555


**Step2: get TFs for each promoter and enhance**
```{r}
df2 <- scGRN_getTF(df1)
head(df2,1)
```

    ##    gene              promoter              enhancer
    ## 1 NLRC3 chr16:3624893-3629893 chr16:3907323-3907699
    ##                                                                                                                                                                                                                                                                                                                    promoter_TF
    ## 1 PPARG, RREB1, ESR1, NR1H2, RXRA, REST, CTCF, EWSR1-FLI1, RARA, ESR2, HSF1, MAFF, MEF2C, MAF, NFE2, NR2C2, SP2, TP63, ZNF263, IRF1, MEF2A, PAX5, SRF, TP53, MAFG, POU4F2, SPI1, T, NR3C1, NR3C2, EGR3, EGR4, GLIS2, SCRT1, SCRT2, ONECUT3, MEF2D, NFKB1, PKNOX1, PKNOX2, CREB3L1, BACH2, CTCFL, MAFK, PBX3, SIX2, PPARA, EGR1
    ##                                                                                                                                                                                                                                                                                                  enhancer_TF
    ## 1 GATA2, GATA3, FOXI1, SPI1, ETS1, STAT1, INSM1, FOXO3, SOX10, CEBPB, FOXP1, POU2F2, RFX5, STAT2, ZNF263, ZEB1, FOXP2, SREBF2, RFX2, ETV6, FOXL1, HIC2, SNAI2, MEIS1, MEIS2, MEIS3, POU1F1, POU2F1, POU3F1, POU3F2, POU3F3, POU3F4, POU5F1B, FIGLA, ID4, TCF3, TCF4, FOXD2, FOXO4, FOXO6, FOXP3, CTCFL, PBX2
 
**Step3: predict TF -> TG via Elastic net**

Notice that the gene_id in gene expression data is ensembl_id so we need to specify the gexpr_gene_id.
```{r}
df3 <- scGRN_getNt(df = df2, gexpr = gexpr,gexpr_gene_id = 'ensembl_gene_id')

# Now change the ensembl_gene_id to hgnc_symbol
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                        dataset="hsapiens_gene_ensembl",
                                        host="uswest.ensembl.org")
gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",values = unique(c(df3$TF,df3$TG)), mart = mart)
df3$TG <- gene_names$hgnc_symbol[match(df3$TG, gene_names$ensembl_gene_id)]
df3$TF <- gene_names$hgnc_symbol[match(df3$TF, gene_names$ensembl_gene_id)]

head(df3)
```
    ##           TG      TF               enhancer               promoter     TFbs      coef
    ## 1  RNF220   FOXI1 chr1:44988522-44988778 chr1:44868550-44873550 enhancer 159.21054
    ## 17 RNF220 BHLHE23 chr1:44988522-44988778 chr1:44868550-44873550 enhancer -11.72905
    ## 20 RNF220   OLIG3 chr1:44988522-44988778 chr1:44868550-44873550 enhancer 148.79653
    ## 26 RNF220    TLX1 chr1:44988522-44988778 chr1:44868550-44873550 promoter  27.85926
    ## 43 RNF220  POU4F2 chr1:44988522-44988778 chr1:44868550-44873550 promoter  24.35724

