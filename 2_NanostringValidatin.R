#####################################
#
#Nanostring Validation between Young/Old 
#Prepared by Darryl Nousome
#2020-11-06
#
####################################

library(readxl)
library(tidyverse)
library(haven)


##read in the nanostring?

setwd("~/Documents/G/Work_Bench/Epidemiology/Projects/Shared 2019 LTF/")

n1=read_xlsx("DATA/Nanostring Data 63 Samples.xlsx",sheet=1)
n2=read_xlsx("DATA/Nanostring Data 63 Samples.xlsx",sheet=2)
n3=read_xlsx("DATA/Nanostring Data 63 Samples.xlsx",sheet=3)


genesrep=c("ERG(Pan)","ERG1/ERG2/ERG3","ERG8","ERG_exon2_fusion","ERG_exon4_fusion","ERG_exon5_fusion", 
        "C-MYC" ,"MYCN", "MAOA", "ANXA2","VEGFA","VEGFR","VEGFR1", "TWIST1", "PSMA/FOLH1", "NPY")

gene_name=n3 %>% slice(-1) %>% filter(`Gene Name` %in% genesrep) %>% select(`Gene Name`)
tn_dt=n3  %>% slice(-1) %>% filter(`Gene Name` %in% genesrep) %>% select(-`Gene Name`) %>% mutate(across(everything(),as.numeric))



##Read the covariate data
covar=read_sas("DATA/ltf_nanostring_merge63.sas7bdat") %>% select(-LTF,-Tumor) %>% distinct() %>%
  mutate(Agecat1=cut(DXAGE,breaks = c(40,60,80))) %>%
  mutate(Agecat2=ifelse(DXAGE>=42 & DXAGE<=58,0,
                        ifelse(DXAGE>=66 & DXAGE<=73,1,NA))) %>%
  mutate(ID=paste0("FP",FP))


genes=tn_dt %>% select(covar$ID)

###Run the LM
res=lapply(1:15,function(x){
  summary(lm(unlist(genes[x,])~covar$Agecat1))$coefficients[2,]
})
res_all_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##CA only
res=lapply(1:15,function(x){
  summary(lm(unlist(genes[x,])[covar$RACE=="Caucasian"]~covar$Agecat1[covar$RACE=="Caucasian"]))$coefficients[2,]
})
res_ca_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##AA Only
res=lapply(1:15,function(x){
  tryCatch({
    summary(lm(unlist(genes[x,])[covar$RACE=="AfricanAmerican"]~covar$Agecat1[covar$RACE=="AfricanAmerican"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_aa_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)


###Age cat2
res=lapply(1:15,function(x){
  tryCatch({
  summary(lm(unlist(genes[x,])~covar$Agecat2))$coefficients[2,]
},error=function(e){
  dt=data.frame(NA,NA,NA,NA)
  names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
  dt})
  })

res_all_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##CA only
res=lapply(1:15,function(x){
  tryCatch({
    summary(lm(unlist(genes[x,])[covar$RACE=="Caucasian"]~covar$Agecat2[covar$RACE=="Caucasian"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})

res_ca_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)


##AA Only
res=lapply(1:15,function(x){
  tryCatch({
    summary(lm(unlist(genes[x,])[covar$RACE=="AfricanAmerican"]~covar$Agecat2[covar$RACE=="AfricanAmerican"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_aa_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)






###Run only the T/N and 
res=lapply(1:15,function(x){
  aa_y=genes[x,][covar$RACE=="AfricanAmerican" & covar$Agecat1 %in% "(40,60]"]
  aa_o=genes[x,][covar$RACE=="AfricanAmerican" & covar$Agecat1 %in% "(60,80]"]
  ca_y=genes[x,][covar$RACE=="Caucasian" & covar$Agecat1 %in% "(40,60]"]
  ca_o=genes[x,][covar$RACE=="Caucasian" & covar$Agecat1 %in% "(60,80]"]
  a_y=genes[x,][covar$Agecat1 %in% "(40,60]"]
  a_o=genes[x,][covar$Agecat1 %in% "(60,80]"]
  tibble(`Mean AA Young`=paste(round(mean(unlist(aa_y),na.rm=T),2),"\u00b1",round(sd(unlist(aa_y),na.rm=T),2)),
         `Mean AA Older`=paste(round(mean(unlist(aa_o),na.rm=T),2),"\u00b1",round(sd(unlist(aa_o),na.rm=T),2)),
         `Mean CA Young`=paste(round(mean(unlist(ca_y),na.rm=T),2),"\u00b1",round(sd(unlist(ca_y),na.rm=T),2)),
         `Mean CA Older`=paste(round(mean(unlist(ca_o),na.rm=T),2),"\u00b1",round(sd(unlist(ca_o),na.rm=T),2)),
         `Mean All Young`=paste(round(mean(unlist(a_y),na.rm=T),2),"\u00b1",round(sd(unlist(a_y),na.rm=T),2)),
         `Mean All Older`=paste(round(mean(unlist(a_o),na.rm=T),2),"\u00b1",round(sd(unlist(a_o),na.rm=T),2)))
  
})

res_mean_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,everything())

res=lapply(1:15,function(x){
  aa_y=genes[x,][covar$RACE=="AfricanAmerican" & covar$Agecat2 %in% 0]
  aa_o=genes[x,][covar$RACE=="AfricanAmerican" & covar$Agecat2 %in% 1]
  ca_y=genes[x,][covar$RACE=="Caucasian" & covar$Agecat2 %in% 0]
  ca_o=genes[x,][covar$RACE=="Caucasian" & covar$Agecat2 %in% 1]
  a_y=genes[x,][covar$Agecat2 %in% 0]
  a_o=genes[x,][covar$Agecat2 %in% 1]
  tibble(`Mean AA Young`=paste(round(mean(unlist(aa_y),na.rm=T),2),"\u00b1",round(sd(unlist(aa_y),na.rm=T),2)),
         `Mean AA Older`=paste(round(mean(unlist(aa_o),na.rm=T),2),"\u00b1",round(sd(unlist(aa_o),na.rm=T),2)),
         `Mean CA Young`=paste(round(mean(unlist(ca_y),na.rm=T),2),"\u00b1",round(sd(unlist(ca_y),na.rm=T),2)),
         `Mean CA Older`=paste(round(mean(unlist(ca_o),na.rm=T),2),"\u00b1",round(sd(unlist(ca_o),na.rm=T),2)),
         `Mean All Young`=paste(round(mean(unlist(a_y),na.rm=T),2),"\u00b1",round(sd(unlist(a_y),na.rm=T),2)),
         `Mean All Older`=paste(round(mean(unlist(a_o),na.rm=T),2),"\u00b1",round(sd(unlist(a_o),na.rm=T),2)))

  })

res_mean_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,everything())


finalres=list(res_all_cat1,res_ca_cat1,res_aa_cat1,res_all_cat2,res_ca_cat2,res_aa_cat2,res_mean_cat1,res_mean_cat2)
finalcovar=covar %>% select(FP,DXAGE,Agecat1,Agecat2,RACE,path_gleason,BCR)






saveRDS(finalres,"../Sharad 2020 Age/finalres_nano.rds")
saveRDS(finalcovar,"../Sharad 2020 Age/finalcovar_nano.rds")

bind_cols(gene_name,genes) %>%
  saveRDS(.,"../Sharad 2020 Age/hm_63.rds")















#####################Nanostring 40 Samples
##Only 29 are available
dt_40=read_xlsx("DATA/Nanostring AA-CA Jocelyn data.xlsx",sheet=2)
d=names(dt_40)[-1:-2]
id=gsub("fp","",d)

covar_40=read_xlsx("DATA/LTF_Nanostring Raw data n=90_01132020.xlsx",sheet=1) %>% filter(FP %in% id) %>% select(FP,RACE,DXAGE) %>%
  mutate(Agecat1=cut(DXAGE,breaks = c(40,60,80))) %>%
  mutate(Agecat2=ifelse(DXAGE>=42 & DXAGE<=58,0,
                        ifelse(DXAGE>=66 & DXAGE<=73,1,NA))) %>%
  mutate(ID=paste0("fp",FP))


gene_name=dt_40 %>% slice(-1) %>% filter(`Gene Name` %in% genesrep) %>% select(`Gene Name`)
tn_dt=dt_40  %>% slice(-1) %>% filter(`Gene Name` %in% genesrep) %>% 
  select(-`Gene Name`,-`Accession #`) %>% 
  mutate(across(everything(),as.numeric))

tn_dt=tn_dt %>% select(covar_40$ID)



###Run the LM
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])~covar_40$Agecat1))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_all_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##CA only
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])[covar_40$RACE=="Caucasian"]~covar_40$Agecat1[covar_40$RACE=="Caucasian"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_ca_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##AA Only
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])[covar_40$RACE=="AfricanAmerican"]~covar_40$Agecat1[covar_40$RACE=="AfricanAmerican"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_aa_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)


###Age cat2
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])~covar_40$Agecat2))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_all_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)

##CA only
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])[covar_40$RACE=="Caucasian"]~covar_40$Agecat2[covar_40$RACE=="Caucasian"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_ca_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)


##AA Only
res=lapply(1:16,function(x){
  tryCatch({
    summary(lm(unlist(tn_dt[x,])[covar_40$RACE=="AfricanAmerican"]~covar_40$Agecat2[covar_40$RACE=="AfricanAmerican"]))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','t value','Pr(>|t|)')
    dt})
})
res_aa_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,T=`t value`,P=`Pr(>|t|)`)



finalres_40=list(res_all_cat1,res_ca_cat1,res_aa_cat1,res_all_cat2,res_ca_cat2,res_aa_cat2)
finalcovar_40=covar_40 %>% select(FP,DXAGE,Agecat1,Agecat2,RACE)



saveRDS(finalres_40,"../Sharad 2020 Age/finalres_nano_40.rds")
saveRDS(finalcovar_40,"../Sharad 2020 Age/finalcovar_nano_40.rds")


bind_cols(gene_name,tn_dt) %>%
  saveRDS(.,"../Sharad 2020 Age/hm_40.rds")

