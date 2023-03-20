#!/usr/bin/env Rscript
rm(list=ls())

setwd("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/")



##packages
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(parallel)
library(ggplot2)
library(reshape2)
library(adegenet)




### NB re named  rmR02SSEFsamples_gapped0.9_maxhet0.6_singleSNP_mac3_DP3_minGQ20_maf0.05_rmindels_biallelic_rmINDIV90%_missingSNP0.6_rmINDIV60%_missingSNP0.7_rmINDIV50% to: mac3_DP3_minGQ20_maf0.05_rmindels_biallelic_rmINDIV90%_missingSNP0.6_rmINDIV60%_missingSNP0.7_rmINDIV50%
# program coudlnt' deal with long name


#loading the data file
#filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.vcf")



#setting populations

pop.data.HG <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_805samples_april2021_HEADER_twoareas.txt", sep = "\t", header = TRUE)
dim(pop.data.HG)

pop.data.HG_2<-pop.data.HG

pop.data.HG_2$pop<-as.factor(pop.data.HG_2$pop)
levels(pop.data.HG_2$pop) # 47 pops 


########
# edit pop map file - remove individuals that were siblings and removed!  -copied from R output for sibs removed and also vcf output of list of individuals with over 30% missing data
##############

#sibgs
#rmsibs<-read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/Rmsibs/populationssnps_ibd_individualstoremovesibs2nddgree.txt")


## individuals missing data

# 90%
#rmINV <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.indv", header = TRUE)

# 70%
#rmINV2 <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.indv", header = TRUE)

# 50%
#rmINV3 <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.indv", header = TRUE)

# 40%
rmINV1 <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/225INDIV_723SNPS_40%missINDV.indv", header = TRUE)
dim(rmINV1) #487   1

#30%
rmINV2 <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/225INDIV_723SNPS_40%missINDV_30%missINDV.indv", header = TRUE)
dim(rmINV2) # 117   1

# join rmINV1 and rmINV2

(rmINV5<-merge(rmINV1,rmINV2))

rmINV5<-dplyr::bind_rows(rmINV1,rmINV2)
dim(rmINV5) #604   1

# R02-SS-EF SAMPLES
#rmINV4 <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/REMOVE_R02-SS-EF-01.txt", header = TRUE)


rmINV5 
rmINV5$INDV

# 90%
#rm row names and column names
#write.table(rmINV$INDV,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINVlist <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.txt")

#70%
#write.table(rmINV2$INDV,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINV2list <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.txt")

#50%
#write.table(rmINV3$INDV,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINV3list <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.txt")



#30%
write.table(rmINV5$INDV,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/225INDIV_723SNPS_30%missINDV.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
rmINV5list <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/225INDIV_723SNPS_30%missINDV.txt")



#join sibs and missing indiv together
#ALLindtoremove<-rbind(rmINVlist,rmINV2list,rmINV3list,rmINV5list)
ALLindtoremove<-rbind(rmINV5list)


#listALLindtoremove<-c(ALLindtoremove)  # do it this way if I want to keep the pops with their original numbers - then use %in% listALLindtoremove$V1 below instead of ALL

#remove these indiv from pop data map list
pop.data.HG_rm <- pop.data.HG[!pop.data.HG$sample.id %in% ALLindtoremove$V1, ]
dim(pop.data.HG_rm) #225   7

head(pop.data.HG_rm)

## add one cluster

pop.data.HG_rm$onecluster<-rep("one")


# three clusters


#pop.data.HG_rm$threeclusters<-pop.data.HG_rm$two_areas


#pop.data.HG_rm$threeclusters[pop.data.HG_rm$threeclusters=="swBC"] <- "CoastalBC"
#pop.data.HG_rm$threeclusters[pop.data.HG_rm$threeclusters=="Northwest"] <- "CoastalBC"

## add two clusters - haida gwaii in one, rest in another


#pop.data.HG_rm <- pop.data.HG[ ! pop.data.HG$sample.id %in% rmsibs$V1, ]
#pop.data.HG_rm

#save shortened file
write.table(pop.data.HG_rm,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/popmap.225_samples_header_region_rmINdiv_sorted.txt",quote=FALSE,sep = "\t",row.names=FALSE)
rm(pop.data.HG)

pop.data.HG<-pop.data.HG_rm 
dim(pop.data.HG) #  225   7

#open shortened file
pop.data.HG <- read.table("F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/popmap.225_samples_header_region_rmINdiv_sorted_3clusters.txt", sep = "\t", header = TRUE)


# count pops and get list
str(pop.data.HG$pop)

pop.data.HG2<-pop.data.HG

pop.data.HG2$pop<-as.factor(pop.data.HG2$pop)

levels(pop.data.HG2$pop)

ponds_kept_225INDIV_723SNPS<-levels(pop.data.HG2$pop)




write.table(ponds_kept_225INDIV_723SNPS,"F:/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/ponds_kept_225INDIV_723SNPS.txt",quote=FALSE,sep = "\t",row.names=FALSE,col.names = F)
