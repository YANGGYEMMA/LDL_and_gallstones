library(data.table)
library(stringr)
library(readxl)
library(dplyr)
library(tibble)
library(scales)
library(ieugwasr)
library(TwoSampleMR)
library(MendelianRandomization)
library(metafor)
library(coloc)
library(mrclust)

setwd("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Files")

####Gene_gs (UKBB)####
results<-data.frame()

#summary statistics#
overall<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_gene_overall.csv")
overall<-overall[,!"V1"]

male<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_gene_male.csv")
male<-male[,!"V1"]

female<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_gene_female.csv")
female<-female[,!"V1"]

summary<-merge(overall,male,by="rsid")
summary<-merge(summary,female,by="rsid")

gene_extract<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/emma_lipid_gall.txt")
gs<-merge(summary,gene_extract,by="rsid")

#overall#
table<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eur_gene-specific (r2=0.1).csv")

for (i in unique(table$drug)) {
  gene<-table[table$drug %in% i,]
  gene_gs<-gs[gs$rsid %in% gene$rsid,]
  
  exposure<-format_data(gene,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                        se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
  outcome<-format_data(gene_gs,type="outcome",snp_col = 'rsid',beta_col = "beta_overall",
                       se_col = "se_overall",effect_allele_col = "alternative_alleles",other_allele_col = "first_allele",
                       eaf_col = "eaf")
  data<-harmonise_data(exposure,outcome,action=1) #no palindromic snps with intermediate af
  data2<-data[data$mr_keep%in% TRUE,]
  
  #ivw#
  correl<-ifelse(nrow(data2)>1,TRUE,FALSE) #get ivw and mr egger estimates with correlations
  mrinput<-dat_to_MRInput(data2,get_correlations = correl,pop = "EUR")
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]],correl = correl)
  dat_tab<-data.frame('overall',i,"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data2)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',i,"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]],correl = correl)
    dat_tab<-data.frame('overall',i,"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
  }
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","table","gs")])}

#male#
for (i in unique(table$drug)) {
  gene<-table[table$drug %in% i,]
  gene_gs<-gs[gs$rsid %in% gene$rsid,]
  
  exposure<-format_data(gene,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                        se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
  outcome<-format_data(gene_gs,type="outcome",snp_col = 'rsid',beta_col = "beta_male",
                       se_col = "se_male",effect_allele_col = "alternative_alleles",other_allele_col = "first_allele",
                       eaf_col = "eaf")
  data<-harmonise_data(exposure,outcome,action=1) #no palindromic snps with intermediate af
  data2<-data[data$mr_keep%in% TRUE,]
  
  #ivw#
  correl<-ifelse(nrow(data2)>1,TRUE,FALSE) #get ivw and mr egger estimates with correlations
  mrinput<-dat_to_MRInput(data2,get_correlations = correl,pop = "EUR")
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]],correl = correl)
  dat_tab<-data.frame('male',i,"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data2)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('male',i,"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]],correl = correl)
    dat_tab<-data.frame('male',i,"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
  }
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","table","gs")])}

#female#
for (i in unique(table$drug)) {
  gene<-table[table$drug %in% i,]
  gene_gs<-gs[gs$rsid %in% gene$rsid,]
  
  exposure<-format_data(gene,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                        se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
  outcome<-format_data(gene_gs,type="outcome",snp_col = 'rsid',beta_col = "beta_female",
                       se_col = "se_female",effect_allele_col = "alternative_alleles",other_allele_col = "first_allele",
                       eaf_col = "eaf")
  data<-harmonise_data(exposure,outcome,action=1) #no palindromic snps with intermediate af
  data2<-data[data$mr_keep%in% TRUE,]
  
  #ivw#
  correl<-ifelse(nrow(data2)>1,TRUE,FALSE) #get ivw and mr egger estimates with correlations
  mrinput<-dat_to_MRInput(data2,get_correlations = correl,pop = "EUR")
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]],correl = correl)
  dat_tab<-data.frame('female',i,"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data2)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('female',i,"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]],correl = correl)
    dat_tab<-data.frame('female',i,"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
  }
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","table","gs")])}

rm(gs,table)
####Gene_gs (FINNGEN)####
results<-data.frame()

gs<-fread("finngen_R8_CHOLELITH_BROAD.gz") 

table<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eur_gene-specific (r2=0.1).csv")

#replace notfound snps with proxies#
proxies<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Proxies/proxies_finngen.csv")
proxies<-proxies[proxies$proxy %in% table$rsid,]

setnames(proxies,"rsID","rsid")
setnames(proxies,"pvalue_GC","pval")

drug<-table[,c("rsid","drug")]
proxies<-merge(proxies,drug,by.x ="proxy",by.y = "rsid")

table<-table[!table$rsid %in% proxies$proxy,]

table<-subset(table,select=-id)
proxies<-subset(proxies,select=-proxy)
table<-rbind(table,proxies)

#analyses#
for (i in unique(table$drug)) {
  gene<-table[table$drug %in% i,]
  gene_gs<-gs[gs$rsids %in% gene$rsid,]
  
  exposure<-format_data(gene,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                        se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
  outcome<-format_data(gene_gs,type="outcome",snp_col = 'rsids',beta_col = "beta",
                       se_col = "sebeta",effect_allele_col = "alt",other_allele_col = "ref",
                       eaf_col = "af_alt",pval_col = "pval")
  data<-harmonise_data(exposure,outcome,action=1) #no palindromic snps with intermediate af
  data2<-data[data$mr_keep%in% TRUE,]
  
  #ivw#
  correl<-ifelse(nrow(data2)>1,TRUE,FALSE) #get ivw and mr egger estimates with correlations
  mrinput<-dat_to_MRInput(data2,get_correlations = correl,pop = "EUR")
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]],correl = correl)
  dat_tab<-data.frame('overall',i,"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data2)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',i,"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]],correl = correl)
    dat_tab<-data.frame('overall',i,"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
  }
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","table","gs")])}

rm(table,gs)
gc()
####Gene_gs (BBJ)####
results<-data.frame()

gs<-fread("GWASsummary_Cholelithiasis_Japanese_SakaueKanai2020.auto.txt.gz")

table<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eas_gene-specific (r2=0.1).csv")

#replace notfound snps with proxies#
proxies<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Proxies_eas/proxies_bbj.csv")
proxies<-proxies[proxies$proxy %in% table$rsid,]

setnames(proxies,"rsID","rsid")
setnames(proxies,"pvalue_GC","pval")

drug<-table[,c("rsid","drug")]
proxies<-merge(proxies,drug,by.x ="proxy",by.y = "rsid")

table<-table[!table$rsid %in% proxies$proxy,]

table<-subset(table,select=-id)
proxies<-subset(proxies,select=-proxy)
table<-rbind(table,proxies)

#analyses#
for (i in unique(table$drug)) {
  gene<-table[table$drug %in% i,]
  gene_gs<-gs[gs$SNPID %in% gene$rsid,]
  
  exposure<-format_data(gene,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                        se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                        eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
  outcome<-format_data(gene_gs,type="outcome",snp_col = 'SNPID',beta_col = "BETA",
                       se_col = "SE",effect_allele_col = "Allele2",other_allele_col = "Allele1",
                       eaf_col = "AF_Allele2",pval_col = "p.value")
  data<-harmonise_data(exposure,outcome,action=1) #no palindromic snps with intermediate af
  data2<-data[data$mr_keep%in% TRUE,]
  
  #ivw#
  correl<-ifelse(nrow(data2)>1,TRUE,FALSE) #get ivw and mr egger estimates with correlations
  mrinput<-dat_to_MRInput(data2,get_correlations = correl,pop = "EAS")
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]],correl = correl)
  dat_tab<-data.frame('overall',i,"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data2)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',i,"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]],correl = correl)
    dat_tab<-data.frame('overall',i,"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
  }
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2","table","gs")])}

rm(table,gs)
gc()
####Gene_gs (META)####
results3<-data.frame()

#FINNGEN#
tab1<-read_excel("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/Revise.xlsx", 
                 sheet = "Gene_gs (FINNGEN)")

#UKBB#
tab2<-read_excel("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/Revise.xlsx", 
                 sheet = "Gene_gs (UKBB)")
tab2<-tab2[tab2$sex %in% "overall",]

#BBJ#
tab3<-read_excel("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/Revise.xlsx", 
                 sheet = "Gene_gs (BBJ)")

#meta-analysis#
for (i in unique(tab1$exp) ){
  
  for (j in c("IVW","Weighted median","MR Egger")) {
    combined1<-tab1[tab1$exp %in% i &tab1$method %in% j,]
    combined2<-tab2[tab2$exp %in% i &tab2$method %in% j,]
    combined3<-tab3[tab3$exp %in% i &tab3$method %in% j,]
    combined<-rbind(combined1,combined2,combined3)
    
    if(nrow(combined)!=0){
      combined$beta<-as.numeric(combined$beta)
      combined$se<-as.numeric(combined$se)
      
      meta<-rma(yi=combined$beta,sei=combined$se,method = "FE")
      model<-ifelse(meta[["QEp"]]<0.05,"REML","FE")
      meta<-rma(yi=combined$beta,sei=combined$se,method = model)
      
      dat_tab<-data.frame("overall",i,"gs",j,
                          meta[["beta"]],meta[["se"]],meta[["ci.lb"]],meta[["ci.ub"]],
                          meta[["pval"]],meta[["k"]],meta[["QEp"]])
      colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","n","pval(heterogeneity)")
      results3<-rbind(results3,dat_tab)}
    
    rm(list=ls()[! ls() %in% c('variants',"results","results2","results3"
                               ,"tab1","tab2","tab3","i")])}
}
rm(tab1,tab2,tab3,i)
gc()
####Coloc (UKBB)####
results<-data.frame()
results2<-data.frame()

for (i in c("pcsk9","apob","abcg58")) {
  
  data_ldl<-read.csv(paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/coloc_ldl_",i,".csv",sep = ""))
  data_gs<-read.csv(paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/coloc_ukbb_",i,".csv",sep = ""))  
  
  data_gs<-na.omit(data_gs)
  data_gs<-data_gs[data_gs$MAF>0.001,]
  
  common<-intersect(data_ldl$snp,data_gs$snp)
  
  data_ldl<-data_ldl[data_ldl$snp %in% common,]
  data_ldl<-data_ldl[!duplicated(data_ldl$snp),]
  
  data_gs<-data_gs[data_gs$snp %in% common,]
  data_gs<-data_gs[!duplicated(data_gs$snp),]
  data_gs$N<-30547+336742
  
  ##check dataset#
  data_ldl<-as.list(data_ldl)
  data_ldl$type<-"quant"
  data_ldl$sdY<-1
  check_dataset(data_ldl)
  
  data_gs<-as.list(data_gs)
  data_gs$type<-"cc"
  data_gs$s<-30547/(30547+336742)
  check_dataset(data_gs)
  
  my.res<-coloc.abf(data_ldl, data_gs, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
  saveRDS(my.res,paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/my.res_ukbb_",i,".rds",sep = ""))
  
  data<-data.frame(i,"gs",1e-05,
                   my.res[["summary"]][["nsnps"]],my.res[["summary"]][["PP.H0.abf"]],
                   my.res[["summary"]][["PP.H1.abf"]],my.res[["summary"]][["PP.H2.abf"]],
                   my.res[["summary"]][["PP.H3.abf"]],my.res[["summary"]][["PP.H4.abf"]])
  colnames(data)<-c("exp","out","p12","snp","H0","H1","H2",'H3',"H4")
  data$conditional_H4<-data$H4/(data$H2+data$H3+data$H4) # the probability of colocalization conditional on the presence of a causal variant for the outcome
  results<-rbind(results,data)
  
  data2<-data.frame(i,"gs",1e-05,
                    my.res[["results"]][["snp"]],my.res[["results"]][["SNP.PP.H4"]])
  colnames(data2)<-c("exp","out","p12","snp","SNP.PP.H4")
  data2<-data2[data2$SNP.PP.H4%in% max(data2$SNP.PP.H4),]
  results2<-rbind(results2,data2)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2")])
  gc()}
####Coloc (FINNGEN)####
results<-data.frame()
results2<-data.frame()

for (i in c("hmgcr","apob","abcg58")) {
  
  data_ldl<-read.csv(paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/coloc_ldl_",i,".csv",sep = ""))
  data_gs<-read.csv(paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/coloc_finngen_",i,".csv",sep = ""))  
  
  data_gs<-data_gs[data_gs$MAF>0.001,]
  
  common<-intersect(data_ldl$snp,data_gs$snp)
  
  data_ldl<-data_ldl[data_ldl$snp %in% common,]
  data_ldl<-data_ldl[!duplicated(data_ldl$snp),]
  
  data_gs<-data_gs[data_gs$snp %in% common,]
  data_gs<-data_gs[!duplicated(data_gs$snp),]
  
  pos<-data_ldl[,c("snp","position")]
  data_gs<-merge(data_gs,pos,by="snp")
  data_gs$N<-34461+301383
  
  ##check dataset#
  data_ldl<-as.list(data_ldl)
  data_ldl$type<-"quant"
  data_ldl$sdY<-1
  check_dataset(data_ldl)
  
  data_gs<-as.list(data_gs)
  data_gs$type<-"cc"
  data_gs$s<-34461/(34461+301383)
  check_dataset(data_gs)
  
  my.res<-coloc.abf(data_ldl, data_gs, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
  saveRDS(my.res,paste("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/my.res_finngen_",i,".rds",sep = ""))
  
  data<-data.frame(i,"gs",1e-05,
                   my.res[["summary"]][["nsnps"]],my.res[["summary"]][["PP.H0.abf"]],
                   my.res[["summary"]][["PP.H1.abf"]],my.res[["summary"]][["PP.H2.abf"]],
                   my.res[["summary"]][["PP.H3.abf"]],my.res[["summary"]][["PP.H4.abf"]])
  colnames(data)<-c("exp","out","p12","snp","H0","H1","H2",'H3',"H4")
  data$conditional_H4<-data$H4/(data$H2+data$H3+data$H4) # the probability of colocalization conditional on the presence of a causal variant for the outcome
  results<-rbind(results,data)
  
  data2<-data.frame(i,"gs",1e-05,
                    my.res[["results"]][["snp"]],my.res[["results"]][["SNP.PP.H4"]])
  colnames(data2)<-c("exp","out","p12","snp","SNP.PP.H4")
  data2<-data2[data2$SNP.PP.H4%in% max(data2$SNP.PP.H4),]
  results2<-rbind(results2,data2)
  
  rm(list=ls()[! ls() %in% c('variants',"results","results2")])
  gc()}
####MRclust (UKBB)####
results<-data.frame()

#summary statistics#
overall<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_ldl_overall.csv")
overall<-overall[,!"V1"]

male<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_ldl_male.csv")
male<-male[,!"V1"]

female<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/summary_statistics_ldl_female.csv")
female<-female[,!"V1"]

summary<-merge(overall,male,by="rsid")
summary<-merge(summary,female,by="rsid")

gene_extract<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/emma_lipid_gall.txt")
gs<-merge(summary,gene_extract,by="rsid")

ldl_clump<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eur_ldl_clump.csv")
ldl_gs<-gs[gs$rsid %in% ldl_clump$rsid,]

#format data#
exposure<-format_data(ldl_clump,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                      se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                      eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
outcome<-format_data(ldl_gs,type="outcome",snp_col = 'rsid',beta_col = "beta_overall",
                     se_col = "se_overall",effect_allele_col = "alternative_alleles",other_allele_col = "first_allele",
                     eaf_col = "eaf")
data<-harmonise_data(exposure,outcome,action=1) #check palindromic snps by hand
data2<-data[data$mr_keep%in% TRUE,]

#per reduction in ldl#
data2$beta.exposure<--data2$beta.exposure

#ivw#
mrinput<-dat_to_MRInput(data2)
mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
dat_tab<-data.frame('overall',"ldl","gs","IVW",mr@Estimate,
                    mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
results<-rbind(results,dat_tab)

mr@Heter.Stat

if(nrow(data2)>2){
  #weighted median#
  mr<-MendelianRandomization::mr_median(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","Weighted median",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  #egger#
  mr<-MendelianRandomization::mr_egger(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","MR Egger",mr@Estimate,mr@StdError.Est,
                      mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                      mr@Pvalue.Int)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)}

# #create data#
# bx = data2$beta.exposure
# by = data2$beta.outcome
# bxse = data2$se.exposure
# byse = data2$se.outcome
# ratio_est<-data2$beta.outcome/data2$beta.exposure
# ratio_est_se<-data2$se.outcome/abs(data2$beta.exposure)
# 
# #mrclust#
# res_em = mr_clust_em(theta = ratio_est,
#                      theta_se = ratio_est_se,
#                      bx = bx,
#                      by = by,
#                      bxse = bxse,
#                      byse = byse,
#                      obs_names = data2$SNP
# )
# saveRDS(res_em,"C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (UKBB).rds")

res_em<-readRDS("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (UKBB).rds")
res<-res_em$results$best

#remove snps with probability <0.8 and cluster with <4 snps#
res<-res[res$probability>=0.8,]

res$cluster_snp<-0
for (k in unique(res$cluster_class)) {
  res[res$cluster_class %in% k,]$cluster_snp<-nrow(res[res$cluster_class %in% k,])
}
res<-res[res$cluster_snp>=4,]

#orientate cluster by effect size#
oc<-sort(unique(res$cluster_class))
oc<-oc[!oc %in% c("Junk","Null")]
om<-sort(unique(res[res$cluster_class %in% oc,]$cluster_mean))

for (j in 1:length(oc)) {
  res[res$cluster_mean==om[j],]$cluster_class<-j}

res$cluster<-res$cluster_class
res[res$cluster=="Null",]$cluster<-length(oc)+1
res[res$cluster=="Junk",]$cluster<-length(oc)+2
res$cluster<-as.numeric(res$cluster)

#cluster-specific#
res<-res[!res$cluster_class %in% c("Junk","Null"),]

for (i in sort(unique(res$cluster_class))) {
  data3<-data2[data2$SNP %in% res[res$cluster_class==i,]$observation,]
  
  #ivw#
  mrinput<-dat_to_MRInput(data3)
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
  dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data3)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)}}

rm(list=ls()[! ls() %in% c('variants',"results","results2")])
gc()
####MRclust (FINNGEN)####
results<-data.frame()

gs<-fread("finngen_R8_CHOLELITH_BROAD.gz") 

ldl_clump<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eur_ldl_clump.csv")

#replace notfound snps with proxies#
proxies<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Proxies/proxies_finngen.csv")
proxies<-proxies[proxies$proxy %in% ldl_clump$rsid,]

setnames(proxies,"rsID","rsid")
setnames(proxies,"pvalue_GC","pval")

ldl_clump<-ldl_clump[!ldl_clump$rsid %in% proxies$proxy,]

ldl_clump<-subset(ldl_clump,select=-id)
ldl_clump<-ldl_clump[,-1]
proxies<-subset(proxies,select=-proxy)
proxies<-proxies[,-1]

ldl_clump<-rbind(ldl_clump,proxies)

#analyses#
ldl_gs<-gs[gs$rsids %in% ldl_clump$rsid,]

#format data#
exposure<-format_data(ldl_clump,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                      se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                      eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
outcome<-format_data(ldl_gs,type="outcome",snp_col = 'rsids',beta_col = "beta",
                     se_col = "sebeta",effect_allele_col = "alt",other_allele_col = "ref",
                     eaf_col = "af_alt",pval_col = "pval")
data<-harmonise_data(exposure,outcome,action=1) #check palindromic snps by hand
data2<-data[data$mr_keep%in% TRUE,]

#per reduction in ldl#
data2$beta.exposure<--data2$beta.exposure

#ivw#
mrinput<-dat_to_MRInput(data2)
mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
dat_tab<-data.frame('overall',"ldl","gs","IVW",mr@Estimate,
                    mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
results<-rbind(results,dat_tab)

mr@Heter.Stat

if(nrow(data2)>2){
  #weighted median#
  mr<-MendelianRandomization::mr_median(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","Weighted median",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  #egger#
  mr<-MendelianRandomization::mr_egger(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","MR Egger",mr@Estimate,mr@StdError.Est,
                      mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                      mr@Pvalue.Int)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)}

# #create data#
# bx = data2$beta.exposure
# by = data2$beta.outcome
# bxse = data2$se.exposure
# byse = data2$se.outcome
# ratio_est<-data2$beta.outcome/data2$beta.exposure
# ratio_est_se<-data2$se.outcome/abs(data2$beta.exposure)
# 
# #mrclust#
# res_em = mr_clust_em(theta = ratio_est,
#                      theta_se = ratio_est_se,
#                      bx = bx,
#                      by = by,
#                      bxse = bxse,
#                      byse = byse,
#                      obs_names = data2$SNP
# )
# saveRDS(res_em,"C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (FINNGEN).rds")

res_em<-readRDS("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (FINNGEN).rds")
res<-res_em$results$best

#remove snps with probability <0.8 and cluster with <4 snps#
res<-res[res$probability>=0.8,]

res$cluster_snp<-0
for (k in unique(res$cluster_class)) {
  res[res$cluster_class %in% k,]$cluster_snp<-nrow(res[res$cluster_class %in% k,])
}
res<-res[res$cluster_snp>=4,]

#orientate cluster by effect size#
oc<-sort(unique(res$cluster_class))
oc<-oc[!oc %in% c("Junk","Null")]
om<-sort(unique(res[res$cluster_class %in% oc,]$cluster_mean))

for (j in 1:length(oc)) {
  res[res$cluster_mean==om[j],]$cluster_class<-j}

res$cluster<-res$cluster_class
res[res$cluster=="Null",]$cluster<-length(oc)+1
res[res$cluster=="Junk",]$cluster<-length(oc)+2
res$cluster<-as.numeric(res$cluster)

#cluster-specific#
res<-res[!res$cluster_class %in% c("Junk","Null"),]

for (i in sort(unique(res$cluster_class))) {
  data3<-data2[data2$SNP %in% res[res$cluster_class==i,]$observation,]
  
  #ivw#
  mrinput<-dat_to_MRInput(data3)
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
  dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data3)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)}}

rm(list=ls()[! ls() %in% c('variants',"results","results2")])
gc()
####MRclust (BBJ)####
results<-data.frame()

gs<-fread("GWASsummary_Cholelithiasis_Japanese_SakaueKanai2020.auto.txt.gz")

ldl_clump<-fread("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/eas_ldl_clump.csv")

#replace notfound snps with proxies#
proxies<-read.csv("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Proxies_eas/proxies_bbj.csv")
proxies<-proxies[proxies$proxy %in% ldl_clump$rsid,]

setnames(proxies,"rsID","rsid")
setnames(proxies,"pvalue_GC","pval")

ldl_clump<-ldl_clump[!ldl_clump$rsid %in% proxies$proxy,]

ldl_clump<-subset(ldl_clump,select=-id)
ldl_clump<-ldl_clump[,-1]
proxies<-subset(proxies,select=-proxy)
proxies<-proxies[,-1]

ldl_clump<-rbind(ldl_clump,proxies)

#analyses#
ldl_gs<-gs[gs$SNPID %in% ldl_clump$rsid,]

#format data#
exposure<-format_data(ldl_clump,type="exposure",snp_col = 'rsid',beta_col = "EFFECT_SIZE",
                      se_col = "SE",effect_allele_col = "ALT",other_allele_col = "REF",
                      eaf_col = "POOLED_ALT_AF",pval_col = "pval",samplesize_col = "N")
outcome<-format_data(ldl_gs,type="outcome",snp_col = 'SNPID',beta_col = "BETA",
                     se_col = "SE",effect_allele_col = "Allele2",other_allele_col = "Allele1",
                     eaf_col = "AF_Allele2",pval_col = "p.value")
data<-harmonise_data(exposure,outcome,action=1) #check palindromic snps by hand
data2<-data[data$mr_keep%in% TRUE,]

#per reduction in ldl#
data2$beta.exposure<--data2$beta.exposure

#ivw#
mrinput<-dat_to_MRInput(data2)
mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
dat_tab<-data.frame('overall',"ldl","gs","IVW",mr@Estimate,
                    mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
results<-rbind(results,dat_tab)

mr@Heter.Stat

if(nrow(data2)>2){
  #weighted median#
  mr<-MendelianRandomization::mr_median(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","Weighted median",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  #egger#
  mr<-MendelianRandomization::mr_egger(mrinput[[1]])
  dat_tab<-data.frame('overall',"ldl","gs","MR Egger",mr@Estimate,mr@StdError.Est,
                      mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                      mr@Pvalue.Int)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)}

# #create data#
# bx = data2$beta.exposure
# by = data2$beta.outcome
# bxse = data2$se.exposure
# byse = data2$se.outcome
# ratio_est<-data2$beta.outcome/data2$beta.exposure
# ratio_est_se<-data2$se.outcome/abs(data2$beta.exposure)
# 
# #mrclust#
# res_em = mr_clust_em(theta = ratio_est,
#                      theta_se = ratio_est_se,
#                      bx = bx,
#                      by = by,
#                      bxse = bxse,
#                      byse = byse,
#                      obs_names = data2$SNP
# )
# saveRDS(res_em,"C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (BBJ).rds")

res_em<-readRDS("C:/Users/Yang Guoyi/Desktop/Files (YANG Guoyi)/Statins and gallstones/R/Revise/res_em (BBJ).rds")
res<-res_em$results$best

#remove snps with probability <0.8 and cluster with <4 snps#
res<-res[res$probability>=0.8,]

res$cluster_snp<-0
for (k in unique(res$cluster_class)) {
  res[res$cluster_class %in% k,]$cluster_snp<-nrow(res[res$cluster_class %in% k,])
}
res<-res[res$cluster_snp>=4,]

#orientate cluster by effect size#
oc<-sort(unique(res$cluster_class))
oc<-oc[!oc %in% c("Junk","Null")]
om<-sort(unique(res[res$cluster_class %in% oc,]$cluster_mean))

for (j in 1:length(oc)) {
  res[res$cluster_mean==om[j],]$cluster_class<-j}

res$cluster<-res$cluster_class
res[res$cluster=="Null",]$cluster<-length(oc)+1
res[res$cluster=="Junk",]$cluster<-length(oc)+2
res$cluster<-as.numeric(res$cluster)

#cluster-specific#
res<-res[!res$cluster_class %in% c("Junk","Null"),]

for (i in sort(unique(res$cluster_class))) {
  data3<-data2[data2$SNP %in% res[res$cluster_class==i,]$observation,]
  
  #ivw#
  mrinput<-dat_to_MRInput(data3)
  mr<-MendelianRandomization::mr_ivw(mrinput[[1]])
  dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","IVW",mr@Estimate,
                      mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
  colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
  results<-rbind(results,dat_tab)
  
  if(nrow(data3)>2){
    #weighted median#
    mr<-MendelianRandomization::mr_median(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","Weighted median",mr@Estimate,
                        mr@StdError,mr@CILower,mr@CIUpper,mr@Pvalue,mr@SNPs,NA)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)
    
    #egger#
    mr<-MendelianRandomization::mr_egger(mrinput[[1]])
    dat_tab<-data.frame('overall',paste("cluster",i,sep = ""),"gs","MR Egger",mr@Estimate,mr@StdError.Est,
                        mr@CILower.Est,mr@CIUpper.Est,mr@Pvalue.Est,mr@SNPs,
                        mr@Pvalue.Int)
    colnames(dat_tab)<-c("sex","exp","out","method","beta","se","lo","up","pval","snps","pval(intercept)")
    results<-rbind(results,dat_tab)}}

rm(list=ls()[! ls() %in% c('variants',"results","results2")])
gc()