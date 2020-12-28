####################################################
#
#GTEx VAlidation?
#
#
#####################################################


library(recount)
library(DESeq2)
library(tidyverse)

setwd("~/Documents/G/Work_Bench/Epidemiology/Projects/Sharad 2020 Age/")

##Download the TCgA data
###http://duffel.rail.bio/recount/v2/TCGA/rse_gene_prostate.Rdata

###read it in
load('~/Documents/xfer/Darryl/rse_gene_prostate.Rdata')
#age strata cuts



##REad in genes
gencode=rtracklayer::import("~/dn/Annotations/gencode/gencode.v25.annotation.gtf")

#gencode %>% dplyr::filter(type=="gene")
gencode=subset(gencode,type=="gene")
#select(gencode,c("gene_name","gene_id"))

g1=data.frame(gencode) %>% dplyr::select(gene_name,gene_id)


##Genes to replicate
rep=c("ERG",'MYCBP','MYCL','MYCNUT','MYCNOS','MYCN','MYCT1','MYC','MYCBP2-AS1','MYCBP2','MYCBP2-AS2','MYCBPAP','MYCLP2','MYCLP1',
      "VEGFA", "HLA-A","HLA-B", "ANXA2", "FOLH1", "NPY","NEDD4L", "MAOA", "TWIST2","TWIST1","TWISTNB",
      "LDHB", "ID4", "PSMA", "ID2")


subrep=g1 %>% filter(gene_name %in% rep)


##Fromt he data 
rse <- scale_counts(rse_gene)

#geochar = lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)
sub=read_tsv("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
## Specify design and switch to DESeq2 forma
#df=colData(rse_gene)

rse$ids_1=sapply(strsplit(df$sampid,"-"),function(x)paste0(x[1],"-",x[2]))

sub=sub %>% filter(SUBJID %in% ids_1)
rse$agecat=sub$AGE



rse$agecat1=as.factor(ifelse(rse$agecat %in% c("40-49","50-59"),0,
                                  ifelse(rse$agecat %in% c("60-69","70-79"),1,NA)))


rse$agecat2=as.factor(ifelse(rse$agecat %in% c("20-29","30-39"),0,
                             ifelse(rse$agecat %in% c("60-69","70-79"),1,NA)))



##Filter only the younger and older
rse_dt1=rse[,!is.na(rse$agecat1)]
rse_dt2=rse[,!is.na(rse$agecat2)]






##Run the test for Age cat1
dds <- DESeqDataSet(rse_dt1, ~ agecat1)
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)


res_cat1=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)

##Run the test for Age cat2
dds <- DESeqDataSet(rse_dt2, ~ agecat2)
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)


res_cat2=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)

