#####################################
#
#Nanostring 
#Prepared by Darryl Nousome
#2020-12-26
#
####################################

library(tidyverse)
library(readxl)
library(haven)
dt_40_tn=read_xlsx("../../Shared 2019 LTF/DATA/Nanostring AA-CA Jocelyn data.xlsx",sheet=2)
dt_40_all=read_xlsx("../../Shared 2019 LTF/DATA/Nanostring AA-CA Jocelyn data.xlsx",sheet=1)
d=names(dt_40_tn)[-1:-2]
id=gsub("fp","",d)

covar_40=read_sas("../../Shared 2019 LTF/DATA/nanostring_01152020.sas7bdat") %>% 
  filter(FP %in% id) %>% 
  arrange(match(FP,id)) %>%
  mutate(Agecat1=cut(DXAGE,breaks = c(40,60,80))) %>%
  mutate(Agecat2=ifelse(DXAGE>=42 & DXAGE<=58,0,
                        ifelse(DXAGE>=66 & DXAGE<=73,1,NA))) %>%
  mutate(ID=paste0("fp",FP))


gene_name=dt_40_tn %>% slice(-1) %>% select(`Gene Name`)
tn_dt=dt_40_tn  %>% slice(-1)  %>% 
  select(-`Gene Name`,-`Accession #`) %>% 
  mutate(across(everything(),as.numeric))

tn_dt=tn_dt %>% select(covar_40$ID)


##Run all

###Run the LM
res=lapply(1:165,function(x){
  tryCatch({
    summary(glm(unlist(covar_40$Agecat1)~unlist(tn_dt[x,]),family = "binomial"))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','z value','Pr(>|z|)')
    dt})
})
res_all_cat1=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,Z=`z value`,P=`Pr(>|z|)`) %>% 
  arrange(P)



res=lapply(1:165,function(x){
  tryCatch({
    summary(glm(unlist(covar_40$Agecat2)~unlist(tn_dt[x,]),family = "binomial"))$coefficients[2,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','z value','Pr(>|z|)')
    dt})
})
res_all_cat2=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,Z=`z value`,P=`Pr(>|z|)`) %>% 
  arrange(P)

###############BCR
res=lapply(1:165,function(x){
  tryCatch({
    summary(glm(unlist(covar_40$BCR)~unlist(tn_dt[x,]),family = "binomial"))$coefficients[2,]
    
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','z value','Pr(>|z|)')
    dt})
})
res_all_bcr=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,Z=`z value`,P=`Pr(>|z|)`) %>% 
  arrange(P)

##Do time to bcr
library(survival)
res=lapply(1:165,function(x){
  tryCatch({
    summary(coxph(Surv(covar_40$time_RP_BCR,unlist(covar_40$BCR))~unlist(tn_dt[x,])))$coefficients[1,]
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA,NA)
    names(dt)=c('coef','exp(coef)','se(coef)','z','Pr(>|z|)')
    dt})
})
res_all_bcr_cox=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate=`exp(coef)`,SE=`se(coef)`,Z=z,P=`Pr(>|z|)`) %>% 
  arrange(P)

##Run for Race
res=lapply(1:165,function(x){
  tryCatch({
    summary(glm(factor(unlist(covar_40$race_all),levels = c("Caucasian","African American"))~unlist(tn_dt[x,]),family = "binomial"))$coefficients[2,]
    
  },error=function(e){
    dt=data.frame(NA,NA,NA,NA)
    names(dt)=c('Estimate','Std. Error','z value','Pr(>|z|)')
    dt})
})
res_all_race=bind_rows(res) %>% mutate(Gene=gene_name$`Gene Name`)%>%
  select(Gene,Estimate,SE=`Std. Error`,Z=`z value`,P=`Pr(>|z|)`) %>% 
  arrange(P)



##Take means of all group
mean_cat1=lapply(1:165,function(y){
    Young=round(mean(unlist(tn_dt[y,covar_40$Agecat1=="(40,60]"])),2)
    Old=round(mean(unlist(tn_dt[y,covar_40$Agecat1=="(60,80]"])),2)
    data.frame(Young,Old)
    })

mean_cat1=bind_rows(mean_cat1) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,everything())%>% mutate(Difference=Old-Young)  %>%
  arrange(desc(abs(Difference)))


mean_cat2=lapply(1:165,function(y){
  Young=round(mean(unlist(tn_dt[y,covar_40$Agecat2 %in%"0"])),2)
  Old=round(mean(unlist(tn_dt[y,covar_40$Agecat2 %in% "1"])),2)
  data.frame(Young,Old)
})

mean_cat2=bind_rows(mean_cat2) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,everything()) %>% mutate(Difference=Old-Young)  %>%
  arrange(desc(abs(Difference)))


mean_bcr=lapply(1:165,function(y){
  BCRNeg=round(mean(unlist(tn_dt[y,covar_40$BCR %in%"0"])),2)
  BCRPos=round(mean(unlist(tn_dt[y,covar_40$BCR %in% "1"])),2)
  data.frame(BCRPos,BCRNeg)
})

mean_bcr=bind_rows(mean_bcr) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,everything()) %>% mutate(Difference=BCRPos-BCRNeg)  %>%
  arrange(desc(abs(Difference)))


mean_race=lapply(1:165,function(y){
  AA=round(mean(unlist(tn_dt[y,covar_40$race_all %in% "African American"])),2)
  CA=round(mean(unlist(tn_dt[y,covar_40$race_all %in% "Caucasian"])),2)
  data.frame(AA,CA)
})

mean_race=bind_rows(mean_race) %>% mutate(Gene=gene_name$`Gene Name`) %>%
  select(Gene,everything()) %>% mutate(Difference=AA-CA)  %>%
  arrange(desc(abs(Difference)))

res_out=list(res_all_cat1,res_all_cat2,res_all_bcr,res_all_bcr_cox,res_all_race)
res_means=list(mean_cat1,mean_cat2,mean_bcr,mean_race)
write_rds(res_out,"../Data/results_race_bcr_age.rds")
write_rds(res_means,"../Data/results_mean_race_bcr_age.rds")

write_rds(covar_40,"../Data/covar_race_bcr_age.rds")
