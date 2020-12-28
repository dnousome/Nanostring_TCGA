##################
#
#
#TCGA Validation on Young vs Older PCa
#By: Darryl Nousome
#
#
##################
library(recount)
library(DESeq2)
library(tidyverse)

setwd("~/Documents/G/Work_Bench/Epidemiology/Projects/Sharad 2020 Age/")

##Download the TCgA data
###http://duffel.rail.bio/recount/v2/TCGA/rse_gene_prostate.Rdata

###read it in
load('~/dn/Projects/PublicData/RNA-seq/rse_gene_prostate.Rdata')


###get the age strata cuts

##Fromt he data 
rse <- scale_counts(rse_gene)

#geochar = lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Specify design and switch to DESeq2 forma
df=colData(rse_gene)
s=data.frame(df) %>% select_if(~ !all(is.na(.)))

##Tumors
#cgc_sample_sample_type_code==1

## Filter out only at sample type code ==1

age_rse=rse[,rse$cgc_sample_sample_type_code %in% 1]

##Extract only the CA men
#
library(readxl)
ea=read_xlsx("~/dn/Projects/PublicData/TCGA_PRAD_ethnicityEIGEN.xlsx",skip=3) %>% dplyr::select(ID=1,Race=4) %>%
  dplyr::filter(Race=="EA")

age_rse$Gleason1=sapply(age_rse$xml_stage_event_gleason_grading,function(x)ifelse(substr(x,1,1)==1,substr(x,3,4),substr(x,2,3)))

##Recode
#age_rse[,age_rse$Gleason=="24"]$xml_stage_event_gleason_grading
#age_rse[,age_rse$Gleason=="54"]$xml_stage_event_gleason_grading

age_rse$Gleason=ifelse(age_rse$Gleason1==24,"6",
                        ifelse(age_rse$Gleason1==33,"6",
                               ifelse(age_rse$Gleason1==44,"8-10",
                                      ifelse(age_rse$Gleason1==53,'8-10',
                        ifelse(age_rse$Gleason1==35,'8-10',
                               ifelse(age_rse$Gleason1==54,'8-10',
                               ifelse(age_rse$Gleason1==45,'8-10',
                                            ifelse(age_rse$Gleason1==55,'8-10',age_rse$Gleason1))))))))


age_rse$Gleason=factor(age_rse$Gleason, levels = c("6", "34", "43","8-10"))

table(age_rse$Gleason,age_rse$Gleason1)

##


##Calculate the Gleason and the Agecategories
###Young (40-60 yrs old) 
#Old     ( 61-80 yrs old) 


#This we used in our current analysis
#Young (42-58 yrs old)
#Old (66-73 yrs Old)


age_rse$age_cat1=as.factor(ifelse(age_rse$xml_age_at_initial_pathologic_diagnosis>=40 &
                                   age_rse$xml_age_at_initial_pathologic_diagnosis<=60 ,0,
                      ifelse(age_rse$xml_age_at_initial_pathologic_diagnosis>=61 &
                               age_rse$xml_age_at_initial_pathologic_diagnosis<=80,1,NA)))

age_rse$age_cat2=as.factor(ifelse(age_rse$xml_age_at_initial_pathologic_diagnosis>=42 &
                                    age_rse$xml_age_at_initial_pathologic_diagnosis<=58 ,0,
                                  ifelse(age_rse$xml_age_at_initial_pathologic_diagnosis>=66 &
                                           age_rse$xml_age_at_initial_pathologic_diagnosis<=73,1,NA)))


age_rse1=age_rse[,!is.na(age_rse$age_cat1)]

###Extract the Age cat and only EA Men
age_rse2=age_rse[,!is.na(age_rse$age_cat2)]
age_rse2=age_rse2[,age_rse2$xml_bcr_patient_barcode %in% ea$ID]

###Load the annotation
##Gencode v25 is used
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
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



##Run the test for Age cat1
dds <- DESeqDataSet(age_rse1, ~ age_cat1 )
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)



res_cat1_unadj=data.frame(res) %>% rownames_to_column("Gene") %>%
    left_join(.,g1,by=c("Gene"="gene_id")) %>%
    arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)




##Run the test for Age cat1 adjusted
dds <- DESeqDataSet(age_rse1, ~ age_cat1 +Gleason)
dds <- DESeq(dds, test = "LRT", reduced = ~Gleason, fitType = "local")
res <- results(dds)



res_cat1_adj=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)




##Age cat 2
dds <- DESeqDataSet(age_rse2, ~ age_cat2 )
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)



res_cat2_unadj=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)




##Run the test for Age cat2 adjusted
dds <- DESeqDataSet(age_rse2, ~ age_cat2 +Gleason)
dds <- DESeq(dds, test = "LRT", reduced = ~Gleason, fitType = "local")
res <- results(dds)



res_cat2_adj=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)



###Stratify into the 2 WD/PD groups
##Select only the agecat 2 with EA men
g_6=age_rse2[,age_rse2$Gleason=="6" | age_rse2$Gleason=="34"]
g_8_10=age_rse2[,age_rse2$Gleason=="8-10" |age_rse2$Gleason=="43"]

##Age cat 2
g_6_c2=g_6[,!is.na(g_6$age_cat2)]
g_8_10_c2=g_8_10[,!is.na(g_8_10$age_cat2)]


###Rerun Split by Gleason
##Gleason 6 and 3+4
#dds <- DESeqDataSet(g_6, ~ age_cat1 )
#dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
#res <- results(dds)

#g6_res_cat1=data.frame(res) %>% rownames_to_column("Gene") %>%
#  left_join(.,g1,by=c("Gene"="gene_id")) %>%
#  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
#  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)


dds <- DESeqDataSet(g_6_c2, ~ age_cat2 )
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)

g6_res_cat2=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)





###4+3 and 8-10
#dds <- DESeqDataSet(g_8_10, ~ age_cat1 )
#dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
#res <- results(dds)

#g8_res_cat1=data.frame(res) %>% rownames_to_column("Gene") %>%
#  left_join(.,g1,by=c("Gene"="gene_id")) %>%
#  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
#  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)


dds <- DESeqDataSet(g_8_10_c2, ~ age_cat2 )
dds <- DESeq(dds, test = "LRT", reduced = ~1, fitType = "local")
res <- results(dds)

g8_res_cat2=data.frame(res) %>% rownames_to_column("Gene") %>%
  left_join(.,g1,by=c("Gene"="gene_id")) %>%
  arrange(pvalue) %>% filter(gene_name %in% rep) %>%
  select(GeneName=gene_name,Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)














final_res=list(res_cat2_unadj,res_cat2_adj)
g_res=list(g6_res_cat2,g8_res_cat2)
covar=colData(age_rse) %>% data.frame() %>% select(Gleason,age_cat1,age_cat2,xml_age_at_initial_pathologic_diagnosis)

saveRDS(final_res,"finalres.rds")
saveRDS(g_res,"gres.rds")
saveRDS(covar,"covar.rds")






##Create heatmap of the data and add annotations?
##Generate the output neccesary to creat the heatmaps

##Age cat 1
dds <- DESeqDataSet(age_rse2, ~ age_cat2 +Gleason)
vst_1=assay(vst(dds),blind=T)


vst_1_genes=vst_1[rownames(vst_1) %in% subrep$gene_id,]
subrep1=subrep %>% arrange(match(gene_id,rownames(vst_1_genes)))
rownames(vst_1_genes)=subrep$gene_name

vst_heat_annot=colData(age_rse2)[,c("Gleason","age_cat2")]
#vst_heat_annot$age_cat1=factor(vst_heat_annot$age_cat1,labels=c("40-60","61-80"))
vst_heat_annot$age_cat2=factor(vst_heat_annot$age_cat2,labels=c("42-58","66-73"))

vst_heat_annot=data.frame(vst_heat_annot)
names(vst_heat_annot)=c("Gleason","Young/Old Cutpoint")

hm_dt=list(vst_heat_annot,vst_1_genes)
saveRDS(hm_dt,"hm_tcga.rds")


#######################################




age=df[,"xml_age_at_initial_pathologic_diagnosis"]
cgc_sample_sample_type_code==1
cgc_case_age_at_diagnosis
hist(age)
sum(age<62)


sum(age>70)


rse1=rse[rse$xml_age_at_initial_pathologic_diagnosis<60|rse$xml_age_at_initial_pathologic_diagnosis>70]

rse1$Age_cat=
#ANXA2