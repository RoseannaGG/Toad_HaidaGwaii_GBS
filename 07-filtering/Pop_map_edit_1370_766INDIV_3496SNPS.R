#!/usr/bin/env Rscript
rm(list=ls())

setwd("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/")



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
#filtered.VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.vcf")



#setting populations
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/1370_Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_HEADER_sorted.txt", sep = "\t", header = TRUE)

########
# edit pop map file - remove individuals that were siblings and removed!  -copied from R output for sibs removed and also vcf output of list of individuals with over 30% missing data
##############

#sibgs
#rmsibs<-read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/populationssnps_ibd_individualstoremovesibs2nddgree.txt")


## individuals missing data

# 90%
#rmINV <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.indv", header = TRUE)

# 70%
#rmINV2 <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.indv", header = TRUE)

# 50%
#rmINV3 <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.indv", header = TRUE)

# 40%
rmINV1 <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_40%missINDV.indv", header = TRUE)
dim(rmINV1) #487   1

#30%
rmINV2 <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_40%missINDV_30%missINDV.indv", header = TRUE)
dim(rmINV2) # 117   1

# join rmINV1 and rmINV2

(rmINV5<-merge(rmINV1,rmINV2))

rmINV5<-dplyr::bind_rows(rmINV1,rmINV2)
dim(rmINV5) #604   1

# R02-SS-EF SAMPLES
#rmINV4 <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/REMOVE_R02-SS-EF-01.txt", header = TRUE)


rmINV5 
rmINV5$INDV

# 90%
#rm row names and column names
#write.table(rmINV$INDV,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINVlist <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS90%missing.txt")

#70%
#write.table(rmINV2$INDV,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINV2list <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS70%missing.txt")

#50%
#write.table(rmINV3$INDV,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
#rmINV3list <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/371INDIV_2392SNPS.listindivmissingSNPS50%missing.txt")



#30%
write.table(rmINV5$INDV,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_30%missINDV.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
rmINV5list <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_30%missINDV.txt")



#join sibs and missing indiv together
#ALLindtoremove<-rbind(rmINVlist,rmINV2list,rmINV3list,rmINV5list)
ALLindtoremove<-rbind(rmINV5list)


#listALLindtoremove<-c(ALLindtoremove)  # do it this way if I want to keep the pops with their original numbers - then use %in% listALLindtoremove$V1 below instead of ALL

#remove these indiv from pop data map list
pop.data_rm <- pop.data[!pop.data$sample.id %in% ALLindtoremove$V1, ]
dim(pop.data_rm) #766   7


## add one cluster

pop.data_rm$onecluster<-rep("one")


#pop.data_rm <- pop.data[ ! pop.data$sample.id %in% rmsibs$V1, ]
#pop.data_rm

#save shortened file
write.table(pop.data_rm,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt",quote=FALSE,sep = "\t",row.names=FALSE)
#rm(pop.data)

pop.data<-pop.data_rm 
dim(pop.data) #  766   7

#open shortened file
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt", sep = "\t", header = TRUE)

# count pops and get list
str(pop.data$pop)

pop.data2<-pop.data

pop.data2$pop<-as.factor(pop.data2$pop)

levels(pop.data2$pop)

ponds_kept_766INDIV_3496SNPS<-levels(pop.data2$pop)

write.table(ponds_kept_766INDIV_3496SNPS,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ponds_kept_766INDIV_3496SNPS.txt",quote=FALSE,sep = "\t",row.names=FALSE,col.names = F)
