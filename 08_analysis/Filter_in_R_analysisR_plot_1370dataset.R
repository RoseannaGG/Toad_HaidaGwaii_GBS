
### R file for analysis of MS Haida Gwaii dataset, author Roseanna Gamlen-Greene, March 2023, roseanna.gamlen.greene@gmail.com ###
## metadata on https://github.com/RoseannaGG/Toad_HaidaGwaii_GBS ##


rm(list=ls())

library(vcfR)
library(adegenet)
library(hierfstat)
library(dartR)
library(ggplot2)
library(spaa)
library(fields)
library(vegan)
library(StAMPP)
library(ggmap)
library(car)
library(VennDiagram)
library(radiator)
library(pegas)
library(RColorBrewer)
library(SNPRelate)
library(reshape2)
library(svglite)
library(dplyr)


setwd("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/")

sink("R_ref_Aboreas_PCA_RGGoutput_766INDIV_3496SNPS_editlocinamesinR_08_03_23.txt", split = TRUE)



###############################################
############# INPUT DATA ################
###################################################

############ CHANGE LOCI NAMES IN VCF #######


VCF_original<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS.vcf")
head(VCF_original_@fix[,'ID'])

VCF_original

"***** Object of Class vcfR *****
766 samples
65 CHROMs
3,496 variants
Object size: 56.3 Mb
0 percent missing data
*****        *****         *****"

VCF_original_test<-VCF_original

#loci names
VCF_original_test@fix[,'ID']

#remove :+ at end of name
VCF_original_test@fix[,'ID']<-gsub("\\:\\+", "", VCF_original_test@fix[,'ID'])
VCF_original_test@fix[,'ID']

#remove :- at end of name
VCF_original_test@fix[,'ID']<-gsub("\\:\\-", "", VCF_original_test@fix[,'ID'])
VCF_original_test@fix[,'ID']

# add ANBO as a prefix
VCF_original_test@fix[,'ID']<-paste("ANBO", sep = '_', VCF_original_test@fix[,'ID'])

VCF_original_test@fix[,'ID']
head(VCF_original_test@fix[,'ID'])

str(VCF_original_test)


"***** Object of Class vcfR *****
766 samples
65 CHROMs
3,496 variants
Object size: 56.3 Mb
0 percent missing data
*****        *****         *****"


#vcfR::write.vcf(VCF_original_test, _original_ = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf.gz")

# then uncompress manually in folder and upload below


#VCF_original_test<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf")


####################################
##### LOAD THE VCF #####
##################################

VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf")





#open shortened file
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt", sep = "\t", header = TRUE)

#pop.data$onecluster<-rep("one")

head(pop.data)

#all(colnames(filtered.VCF@gt)[-1] == pop.data$seID)

### genlight from vcf #####
gl.toad <- vcfR2genlight(VCF)

## add region pop as pop to genlight object
pop(gl.toad) <- pop.data$fourclusters

#setting ploidy
ploidy(gl.toad) <-2
pop(gl.toad)


gl.toad

" /// GENLIGHT OBJECT /////////

 // 766 genotypes,  3,496 binary SNPs, size: 3.4 Mb
 349205 (13.04 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-287)
   @other: a list containing: elements without names "



####################################
################# MAF ###################
##############################
#loading the data file
VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf")

#minor allele freq
mafvcf_766INDIV_3496SNPS<-maf(VCF,element=2)



mafvcf_766INDIV_3496SNPS.df<-as.data.frame(mafvcf_766INDIV_3496SNPS)

head(mafvcf_766INDIV_3496SNPS.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))

# MAFplot_766INDIV_3496SNPS_2
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_766INDIV_3496SNPS_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf_766INDIV_3496SNPS.df$Frequency,breaks=seq(0,1,l=50), main="MAF whole population vcf file maf()",
     xlab="Minor allele frequency",
     ylab="Frequency",
     xlim=(0,3500))
dev.off()


### subset to show below 0.05
mafvcf_766INDIV_3496SNPS.df_below0.05<-mafvcf_766INDIV_3496SNPS.df[mafvcf_766INDIV_3496SNPS.df$Frequency<=0.05,]
dim(mafvcf_766INDIV_3496SNPS.df_below0.05) # 1526    4
dim(mafvcf_766INDIV_3496SNPS.df)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_766INDIV_3496SNPS_below0.05_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf_766INDIV_3496SNPS.df_below0.05$Frequency,breaks=seq(0,0.05,l=50), main="MAF metapop below 0.05",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()



### hg

VCF.HG<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/HG_766INDIV_3496SNPS_editlocinamesinR_manualVCFsubsetHG.vcf")

#minor allele freq
mafvcf_766INDIV_3496SNPS.HG<-maf(VCF.HG,element=2)



mafvcf_766INDIV_3496SNPS.HG.df<-as.data.frame(mafvcf_766INDIV_3496SNPS.HG)

hist(mafvcf_766INDIV_3496SNPS.HG.df$Frequency,breaks=seq(0,1,l=50), main="MAF Haida Gwaii",
     xlab="Minor allele frequency",
     ylab="Frequency")

head(mafvcf_766INDIV_3496SNPS.HG.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_766INDIV_3496SNPS_239INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf_766INDIV_3496SNPS.HG.df$Frequency,breaks=seq(0,1,l=50), main="MAF Haida Gwaii",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

### subset to show below 0.03
mafvcf_766INDIV_3496SNPS.HG.df_below0.03<-mafvcf_766INDIV_3496SNPS.HG.df[mafvcf_766INDIV_3496SNPS.HG.df$Frequency<=0.03,]
dim(mafvcf_766INDIV_3496SNPS.HG.df_below0.03) # 3280    4
dim(mafvcf_766INDIV_3496SNPS.HG.df) # 3496    4


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_766INDIV_3496SNPS_239INDIV_3496SNPS_below0.03_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(mafvcf_766INDIV_3496SNPS.HG.df_below0.03$Frequency,breaks=seq(0,0.03,l=50), main="MAF Haida Gwaii below 0.03",
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()


###############################
####### He vs Ho  ##########
#############################



########## H obs as four clusters #####

#### genind dartR from genligght ####
genind.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)
genind.toad
""

pop(genind.toad)<-pop.data$fourclusters
pop(genind.toad)

### check ###


#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genind.toad)

#hierf.toad<-genind2hierfstat(genind.toad,pop=TRUE)
#hierf.toad$pop<-pop.data$pop

hierf.toadbasic<-basic.stats(hierf.toad)

hierf.toadbasic$Ho

head(hierf.toadbasic$Ho)

Hodf<-as.data.frame(hierf.toadbasic$Ho)

Hodf$locus<-row.names(Hodf)

head(Hodf)

Hodf$locus<- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', Hodf$locus))))

head(Hodf)

## rename 
colnames(Hodf) <- c("VanIsland_Hobs", "LowerMain_Hobs","HaidaGwai_Hobs","Northwest_Hobs","locus")  
head(Hodf)


colnames(Hodf) <- c("VanIsland", "LowerMain","HaidaGwai","Northwest","locus")  
head(Hodf)




### turn wide to long
df_long_Hobs <- tidyr::gather(Hodf,
                              key = Region,
                              value = Obs_Het,
                              VanIsland ,  LowerMain ,  HaidaGwai,Northwest)

head(df_long_Hobs)
dim(df_long_Hobs)

str(df_long_Hobs)
summary(df_long_Hobs)



##################################
###########  Ext Het m- max obs het as one ####################
############################

# isn't working atm...


#### split matrix into 4 regions ###########

pop(gl.toad) <- pop.data$fourclusters
pop.data$fourclusters

gl.toadbypop<-seppop(gl.toad)

"$HaidaGwai
 /// GENLIGHT OBJECT /////////

 // 239 genotypes,  3,496 binary SNPs, size: 1.3 Mb
 114335 (13.68 %) missing data

 // Basic content
   @gen: list of 239 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  239 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 239-239)
   @other: a list containing: elements without names 


$LowerMain
 /// GENLIGHT OBJECT /////////

 // 284 genotypes,  3,496 binary SNPs, size: 1.4 Mb
 122283 (12.32 %) missing data

 // Basic content
   @gen: list of 284 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  284 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 284-284)
   @other: a list containing: elements without names 


$Northwest
 /// GENLIGHT OBJECT /////////

 // 48 genotypes,  3,496 binary SNPs, size: 479.1 Kb
 19867 (11.84 %) missing data

 // Basic content
   @gen: list of 48 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  48 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-48)
   @other: a list containing: elements without names 


$VanIsland
 /// GENLIGHT OBJECT /////////

 // 195 genotypes,  3,496 binary SNPs, size: 1.1 Mb
 92720 (13.6 %) missing data

 // Basic content
   @gen: list of 195 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  195 individual labels
   @loc.names:  3496 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-195)
   @other: a list containing: elements without names 

"





#pop(gl.toadbypop) <- pop.data$fourclusters


gl.toadVanIsland<-gl.toadbypop$VanIsland
gl.toadLowerMain<-gl.toadbypop$LowerMain
gl.toadHaidaGwai<-gl.toadbypop$HaidaGwai
gl.toadNorthwest<-gl.toadbypop$Northwest




### convert genlight to genind objects 
genind.toadVanIsland<-dartR::gl2gi(gl.toadVanIsland, probar = FALSE, verbose = NULL)


""


genind.toadLowerMain<-dartR::gl2gi(gl.toadLowerMain, probar = FALSE, verbose = NULL)

""


genind.toadHaidaGwai<-dartR::gl2gi(gl.toadHaidaGwai, probar = FALSE, verbose = NULL)

''
genind.toadNorthwest<-dartR::gl2gi(gl.toadNorthwest, probar = FALSE, verbose = NULL)

"Error in .nextMethod(x = x, value = value) : 
  Vector length does no match number of loci
In addition: Warning message:
In df2genind(xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names,  :
  Markers with no scored alleles have been removed"

############### summary on each one
genind.toadVanIsland_sum<-adegenet::summary(genind.toadVanIsland)

genind.toadLowerMain_sum<-adegenet::summary(genind.toadLowerMain)

genind.toadHaidaGwai_sum<-adegenet::summary(genind.toadHaidaGwai)

genind.toadNorthwest_sum<-adegenet::summary(genind.toadNorthwest) # not working


###### Hexp per region ############

## rename
divVanIsland<-genind.toadVanIsland_sum
divLowerMain<-genind.toadLowerMain_sum
divHaidaGwai<-genind.toadHaidaGwai_sum
divNorthwest<-genind.toadNorthwest_sum


divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp
divNorthwest$Northwest_Hexp<-divNorthwest$Hexp


divallfour<-cbind(divLowerMain$LowerMain_Hexp,divVanIsland$VanIsland_Hexp,divHaidaGwai$HaidaGwai_Hexp,divNorthwest$Northwest_Hexp)
head(divallfour)

## rename column names
names(divallfour)
colnames(divallfour) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp","Northwest_Hexp")    # Applying colnames
head(divallfour)
str(divallfour)

# turn to dataframe
dfdivallfour<-as.data.frame(divallfour)
head(dfdivallfour)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallfour$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallfour$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallfour$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallfour$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallfour$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallfour$HaidaGwai_Hexp) 


Northwest_mean<-mean(dfdivallfour$Northwest_Hexp) #
Northwest_sd<-sd(dfdivallfour$Northwest_Hexp) 

### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd,Northwest_mean,Northwest_sd)

#write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_601INDIV_173SNPS_maxhetobsasone_nohwe.txt")



dfdivallfour$locus<-row.names(dfdivallfour)
head(dfdivallfour)


colnames(dfdivallfour) <- c("LowerMain", "VanIsland","HaidaGwai","Northwest","locus")  
head(dfdivallfour)

### turn wide to long
df_long <- tidyr::gather(dfdivallfour,
                         key = Region,
                         value = Expect_Het,
                         LowerMain ,  VanIsland ,  HaidaGwai,Northwest)

head(df_long)
dim(df_long)

df_long_Hexp<-df_long
head(df_long_Hexp)



######### join Hobs and Hexp ####

Ho_He_df<-merge(x = df_long_Hobs, y = df_long_Hexp, by = c("locus","Region"), all = TRUE)

dim(Ho_He_df) # 360

head(Ho_He_df)
#View(Ho_He_df)

plot(Ho_He_df$Expect_Het,Ho_He_df$Obs_Het)

## rename regions

labels <- c(VanIsland = "Vancouver Island", LowerMain = "Lower mainland",HaidaGwai="Haida Gwaii",Northwest= "Northwest BC")

p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                     name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p



## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/He_vs_Ho_per_region_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Expect_Het, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_point()+
  facet_grid(Region ~ .,labeller=labeller(Region = labels))+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Expected heterozygosity") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
# scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
# ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()



## just Ho
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Ho_per_region_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Region, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

# just He
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/He_per_region_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(Ho_He_df, aes(x=Region, y=Expect_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()



## just Ho on hobs dataset alone
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Ho_per_region_766INDIV_3496SNPS_2.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long_Hobs, aes(x=Region, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()




##############################
####### CALL RATE #######
#################################

genid.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)

"/// GENIND OBJECT /////////

 // 766 individuals; 3,496 loci; 6,992 alleles; size: 22.5 Mb

 // Basic content
   @tab:  766 x 6992 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 6992 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 48-284)"

genid.toad.df<-genid.toad@tab
str(genid.toad.df)
genid.toad.df<-as.data.frame(genid.toad.df)
head(genid.toad.df)



genid.toad.df[is.na(genid.toad.df)] <- "00"


SNPscalled <- data.frame("Internal_ID"= row.names(genid.toad.df), "Called"=rowSums(genid.toad.df != "00"), "Call.Rate"=(rowSums(genid.toad.df != "00")/ncol(genid.toad.df)))
SNPscalled
hist(SNPscalled$Call.Rate)
hist(SNPscalled$Call.Rate,breaks=seq(0.60,1,l=100))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/CallRatePerPond_plot_766INDIV_3496SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(SNPscalled$Call.Rate,breaks=seq(0.60,1,l=100))
dev.off()







###########################################################
####################### LD Filtering ################
##############################################################
# haven't figured out how to remove snps


vcfF22 <- "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf"
snpgdsVCF2GDS(vcfF22, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.gds", method="biallelic.only")

snpgdsSummary("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.gds")
""

snpgdsClose("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.gds")

genofile <- SNPRelate::snpgdsOpen("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.gds")

# LD Pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, slide.max.bp = 100000, autosome.only = FALSE, ld.threshold = 0.2)

#3,390 markers are selected in total.

#### USE THIS ONE
#increase to 500,000 - composite and others!!
set.seed(1000)
snpset1 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, method="corr", autosome.only = FALSE, ld.threshold = 0.2)

# 3,213 markers are selected in total.


## try method = r 500k
set.seed(1000)
snpset2 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, autosome.only = FALSE, method="r",ld.threshold = 0.2)

# 3,270 markers are selected in total.


#dprime 500k
set.seed(1000)
snpset3 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, autosome.only = FALSE, method="dprime",ld.threshold = 0.2)

# 2,747 markers are selected in total.





#increase to 500,000 - composite and others!!
set.seed(1000)
snpset1 <- snpgdsLDpruning(genofile, slide.max.bp = 500000, method="corr", autosome.only = FALSE, ld.threshold = 0.2)

# 3,213 markers are selected in total. (from 3496) = removed 283 SNPs in high LD

str(snpset1)


"List of 65
 $ chrScaffold_1__2_contigs__length_508792519 : int [1:342] 1 2 3 5 6 7 8 9 10 11 ...
 $ chrScaffold_2__3_contigs__length_804548876 : int [1:519] 373 374 375 376 377 378 379 380 381 382 ...
 $ chrScaffold_3__2_contigs__length_655571293 : int [1:430] 934 935 936 937 938 939 940 941 942 943 ...
 $ chrScaffold_4__2_contigs__length_776011713 : int [1:539] 1398 1399 1400 1401 1402 1403 1404 1405 1407 1408 ...
 $ chrScaffold_5__2_contigs__length_593928783 : int [1:416] 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 ...
 $ chrScaffold_6__2_contigs__length_328854778 : int [1:230] 2437 2438 2439 2440 2441 2442 2443 2444 2445 2446 ...
 $ chrScaffold_7__1_contigs__length_206032987 : int [1:184] 2689 2690 2691 2694 2695 2696 2697 2698 2699 2700 ...
 $ chrScaffold_8__1_contigs__length_198584680 : int [1:159] 2897 2898 2899 2900 2901 2903 2904 2905 2906 2907 ...
 $ chrScaffold_9__1_contigs__length_196030596 : int [1:180] 3065 3066 3067 3068 3069 3071 3072 3073 3074 3075 ...
 $ chrScaffold_10__1_contigs__length_179080748: int [1:154] 3260 3262 3263 3264 3266 3267 3268 3269 3270 3271 ...
 $ chrScaffold_15__1_contigs__length_349370   : int 3435
 $ chrScaffold_20__1_contigs__length_332182   : int 3436
 $ chrScaffold_23__1_contigs__length_311584   : int 3437
 $ chrScaffold_33__1_contigs__length_253271   : int 3438
 $ chrScaffold_35__1_contigs__length_248379   : int 3439
 $ chrScaffold_37__1_contigs__length_251072   : int 3440
 $ chrScaffold_43__1_contigs__length_235388   : int 3441
 $ chrScaffold_53__1_contigs__length_217242   : int 3442
 $ chrScaffold_59__1_contigs__length_206510   : int 3443
 $ chrScaffold_61__1_contigs__length_203537   : int [1:2] 3445 3446
 $ chrScaffold_87__1_contigs__length_164672   : int 3447
 $ chrScaffold_93__1_contigs__length_161056   : int 3448
 $ chrScaffold_103__1_contigs__length_146291  : int 3449
 $ chrScaffold_116__1_contigs__length_134014  : int 3450
 $ chrScaffold_128__1_contigs__length_128476  : int [1:2] 3451 3452
 $ chrScaffold_131__1_contigs__length_127423  : int 3453
 $ chrScaffold_136__1_contigs__length_125746  : int 3454
 $ chrScaffold_141__1_contigs__length_124535  : int 3455
 $ chrScaffold_146__1_contigs__length_121756  : int 3457
 $ chrScaffold_156__1_contigs__length_117762  : int 3458
 $ chrScaffold_162__1_contigs__length_115834  : int 3459
 $ chrScaffold_164__1_contigs__length_115762  : int 3460
 $ chrScaffold_182__1_contigs__length_108084  : int 3461
 $ chrScaffold_211__1_contigs__length_100455  : int 3462
 $ chrScaffold_220__1_contigs__length_98461   : int 3463
 $ chrScaffold_257__1_contigs__length_89402   : int [1:3] 3464 3465 3466
 $ chrScaffold_263__1_contigs__length_87712   : int 3467
 $ chrScaffold_264__1_contigs__length_87568   : int 3468
 $ chrScaffold_286__1_contigs__length_84388   : int 3469
 $ chrScaffold_317__1_contigs__length_78170   : int 3470
 $ chrScaffold_397__1_contigs__length_67266   : int 3471
 $ chrScaffold_543__1_contigs__length_53344   : int 3472
 $ chrScaffold_559__1_contigs__length_52607   : int 3473
 $ chrScaffold_570__1_contigs__length_51877   : int 3474
 $ chrScaffold_638__1_contigs__length_47705   : int 3475
 $ chrScaffold_754__1_contigs__length_41843   : int [1:2] 3476 3477
 $ chrScaffold_835__1_contigs__length_38384   : int 3478
 $ chrScaffold_883__1_contigs__length_36667   : int 3479
 $ chrScaffold_916__1_contigs__length_35273   : int 3480
 $ chrScaffold_933__1_contigs__length_34700   : int 3481
 $ chrScaffold_967__1_contigs__length_33574   : int 3482
 $ chrScaffold_1008__1_contigs__length_32738  : int 3483
 $ chrScaffold_1011__1_contigs__length_32652  : int 3484
 $ chrScaffold_1138__1_contigs__length_29745  : int 3485
 $ chrScaffold_1166__1_contigs__length_29151  : int 3486
 $ chrScaffold_1194__1_contigs__length_28640  : int 3487
 $ chrScaffold_1271__1_contigs__length_27132  : int 3488
 $ chrScaffold_1447__1_contigs__length_24630  : int 3489
 $ chrScaffold_1914__1_contigs__length_19482  : int 3490
 $ chrScaffold_2127__1_contigs__length_17606  : int 3491
 $ chrScaffold_2129__1_contigs__length_17586  : int 3492
 $ chrScaffold_2462__1_contigs__length_15147  : int 3493
 $ chrScaffold_2540__1_contigs__length_14532  : int 3494
 $ chrScaffold_3333__1_contigs__length_9447   : int 3495
 $ chrScaffold_3438__1_contigs__length_8665   : int 3496"


## try to get names of loci to keep

#turn into dataframe... doesn't work
snpset1.df<-as.data.frame(snpset1)

"
Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : 
                      arguments imply differing number of rows: 342, 519, 430, 539, 416, 230, 184, 159, 180, 154, 1, 2, 3 
                      "

LD.keep.loci.names <- rownames(snpset1.df)

#subset
gl.toad.LD <- gl.toad[ , LD.keep.loci.names]
gl.toad.LD



#######################
## OBS HET FILTERING  ####
######################

#########################################
########## OBS HET PER REGION FILTER #####
####################################

#### genind dartR from genligght ####
genind.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)
genind.toad
""

pop(genind.toad)<-pop.data$fourclusters
pop(genind.toad)

### check ###


#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genind.toad)

#hierf.toad<-genind2hierfstat(genind.toad,pop=TRUE)
#hierf.toad$pop<-pop.data$pop

hierf.toadbasic<-basic.stats(hierf.toad)

hierf.toadbasic$Ho

head(hierf.toadbasic$Ho)

Hodf<-as.data.frame(hierf.toadbasic$Ho)

Hodf$locus<-row.names(Hodf)

head(Hodf)

Hodf$locus<- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', Hodf$locus))))

head(Hodf)

## rename 
colnames(Hodf) <- c("VanIsland_Hobs", "LowerMain_Hobs","HaidaGwai_Hobs","Northwest_Hobs","locus")  
head(Hodf)
dim(Hodf)


colnames(Hodf) <- c("VanIsland", "LowerMain","HaidaGwai","Northwest","locus")  
head #3496    5

str(Hodf) 



### bec 8th march 
ponds <- colnames(Hodf)[ 1 : (ncol(Hodf)-1) ]


# SNPS TO KEEP
Hodf.filt.OBSHET <- Hodf %>% filter( if_all(ponds, function(x) x < 0.6) )
dim(Hodf.filt.OBSHET) # 3356    5

Hodf.filt.OBSHET_loci<-Hodf.filt.OBSHET[,5]
str(Hodf.filt.OBSHET_loci) #  chr [1:3356]


#snps TO KEEP
gl.toad.obshetPerregion<- gl.toad[ , Hodf.filt.OBSHET_loci]
gl.toad.obshetPerregion

"  /// GENLIGHT OBJECT /////////

 // 766 genotypes,  3,356 binary SNPs, size: 3.3 Mb
 335947 (13.07 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  3356 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-284)
   @other: a list containing: elements without names  "



####### convert genlight to genind
genind.toad.obshetPerregion<-dartR::gl2gi(gl.toad.obshetPerregion, probar = FALSE, verbose = NULL)


###########################################
#### convert genind to structure ####
###########################################

########### genind2structureL ######
#https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R

source("F:/GBS_data_03_02_21/genind2structureL_function.R")



genind2structureL(genind.toad.obshetPerregion, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/structure_766INDIV_3356SNPS_genind_toad_obshetPerregion_regionID.str", pops=TRUE)



## convert structure file to vcf using pgdspider



## import filtered vcf

#VCF.obshet06perregionR<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf")


## convert vcf to genlight


### alternatively edit VCF in R!!! Then no need to go back to str and then to vcf


#############################################
#### SUBSET VCF SNP LOCI ID ######
#######################################
VCFfiletest

"***** Object of Class vcfR *****
766 samples
65 CHROMs
3,496 variants
Object size: 56.3 Mb
0 percent missing data
*****        *****         *****"

# LIST LOCI TO KEEP
str(Hodf.filt.OBSHET_loci) # chr [1:3356]

# SUBSET TO KEEP JUST LOCI IN THIS LIST
VCF.obshet06perregionR<-subset(VCFfiletest, VCFfiletest@fix[,'ID'] %in% Hodf.filt.OBSHET_loci)

VCF.obshet06perregionR

"***** Object of Class vcfR *****
766 samples
56 CHROMs
3,356 variants
Object size: 52.8 Mb
0 percent missing data
*****        *****         *****"

VCF.obshet06perregionR



#vcfR::write.vcf(VCF.obshet06perregionR, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf.gz")


## import filtered vcf

VCF.obshet06perregionR<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf")



## convert vcf to genlight



#open shortened file
pop.data <- read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.766_samples_header_region_rmINdiv_sorted.txt", sep = "\t", header = TRUE)

#pop.data$onecluster<-rep("one")

head(pop.data)

#all(colnames(filtered.VCF@gt)[-1] == pop.data$seID)

### genlight from vcf #####
gl.toad.obshet06perregionR <- vcfR2genlight(VCF.obshet06perregionR)

## add region pop as pop to genlight object
pop(gl.toad) <- pop.data$fourclusters

#setting ploidy
ploidy(gl.toad.obshet06perregionR) <-2
pop(gl.toad.obshet06perregionR)


gl.toad.obshet06perregionR




##### check obs het after filtering ######


########## H obs as four clusters #####

#### genind dartR from genligght ####
genind.toad.obshetPerregion<-dartR::gl2gi(gl.toad.obshetPerregion, probar = FALSE, verbose = NULL)
genind.toad.obshetPerregion
""

pop(genind.toad.obshetPerregion)<-pop.data$fourclusters
pop(genind.toad.obshetPerregion)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion<-genind2hierfstat(genind.toad.obshetPerregion)

#hierf.toad.obshetPerregion<-genind2hierfstat(genind.toad.obshetPerregion,pop=TRUE)
#hierf.toad.obshetPerregion$pop<-pop.data$pop

hierf.toad.obshetPerregionbasic<-basic.stats(hierf.toad.obshetPerregion)

hierf.toad.obshetPerregionbasic$Ho

head(hierf.toad.obshetPerregionbasic$Ho)

Hobsdf_afterfilter<-as.data.frame(hierf.toad.obshetPerregionbasic$Ho)

Hobsdf_afterfilter$locus<-row.names(Hobsdf_afterfilter)

head(Hobsdf_afterfilter)

Hobsdf_afterfilter$locus<- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', Hobsdf_afterfilter$locus))))

head(Hobsdf_afterfilter)

## rename 
colnames(Hobsdf_afterfilter) <- c("VanIsland_Hobs", "LowerMain_Hobs","HaidaGwai_Hobs","Northwest_Hobs","locus")  
head(Hobsdf_afterfilter)


colnames(Hobsdf_afterfilter) <- c("VanIsland", "LowerMain","HaidaGwai","Northwest","locus")  
head(Hobsdf_afterfilter)
dim(Hobsdf_afterfilter) # 3356



### turn wide to long
long_Hobsdf_afterfilter <- tidyr::gather(Hobsdf_afterfilter,
                                         key = Region,
                                         value = Obs_Het,
                                         VanIsland ,  LowerMain ,  HaidaGwai,Northwest)

head(long_Hobsdf_afterfilter)
dim(long_Hobsdf_afterfilter)

str(long_Hobsdf_afterfilter)
summary(long_Hobsdf_afterfilter)



## just Ho on hobs dataset alone
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Ho_per_region_766INDIV_3356SNPS_obshet06perRegion.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(long_Hobsdf_afterfilter, aes(x=Region, y=Obs_Het,fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()







##########################################
####### RM SIBS ###########
#####################################

# run plink on server using 766INDIV_3356SNPS_editlocinamesinR_obshet06Perregion_R.vcf file to generate list of sibs to remove

#list sibs to remove
listindiv_plink_sub_sib_rel_over_04_names_unique.char<-read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/766INDIV_3356SNPS_listindiv_plink_sub_sib_rel_over_04_names_unique.char.txt")



write.table(listindiv_plink_sub_sib_rel_over_04_names_unique.char, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/listindiv_plink_sub_sib_rel_over_04_names_unique.indv", row.names = FALSE, col.names = FALSE,quote=FALSE )



listindiv_plink_sub_sib_rel_over_04_names_unique.char$V1

listsibs<-listindiv_plink_sub_sib_rel_over_04_names_unique.char$V1


# function alba made to remove siblings
source("F:/GBS_data_03_02_21/Function_only_remove_siblings_genlight_Alba_nov_8th_2021.R")

# run function
gl.toad.obshetPerregion.nosibs<- remove_sibs_genlight(genlight_object = gl.toad.obshetPerregion,
                                                      names = listsibs)
gl.toad.obshetPerregion.nosibs


" 
/// GENLIGHT OBJECT /////////

 // 510 genotypes,  3,356 binary SNPs, size: 2.7 Mb
 259114 (15.14 %) missing data

 // Basic content
   @gen: list of 510 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  510 individual labels
   @loc.names:  3356 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 34-204)
   @other: a list containing: elements without names 
"

### edit pop file

pop.data<-pop.data[which(!pop.data$sample.id %in% listsibs),]
dim(pop.data) # 510   7

write.table(pop.data,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.510_samples_header_region_rmINdiv_rmsibs_sorted.txt",quote=FALSE,sep = "\t",row.names=FALSE)

# count pops and get list
str(pop.data$pop)

pop.data2<-pop.data

pop.data2$pop<-as.factor(pop.data2$pop)

levels(pop.data2$pop)

ponds_kept_510INDIV_3356SNPS<-levels(pop.data2$pop)
ponds_kept_510INDIV_3356SNPS



#######################
## HWE filtering per pond p=0.01 ####
######################
#https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html


## check new pop.data file
dim(pop.data) #510

pop(gl.toad.obshetPerregion.nosibs)

#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs, probar = FALSE, verbose = NULL)

"/// GENIND OBJECT /////////

 // 510 individuals; 3,356 loci; 6,712 alleles; size: 15 Mb

 // Basic content
   @tab:  510 x 6712 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 6712 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 34-204)"


library("pegas")
#(nanhwe.full <- hw.test(genid.toad.obshetPerregion.nosibs, B = 1000)) # performs 1000 permuatations

## hwe per pop

pop(genid.toad.obshetPerregion.nosibs)<-pop.data$pop

### check how many pops
pop(genid.toad.obshetPerregion.nosibs) # 44

#run hwe test for each pop seprately - and only focus on the analytical P value (B=0)
(nanhwe.obshetPerregion.pop <- seppop(genid.toad.obshetPerregion.nosibs) %>% lapply(hw.test, B = 0))

# Take the third column with all rows - p value column
(nanhwe.obshetPerregion.mat <- sapply(nanhwe.obshetPerregion.pop, "[", i = TRUE, j = 3)) 


#turn to dataframe
nanhwe.obshetperregion.mat.df<-as.data.frame(nanhwe.obshetPerregion.mat)

str(nanhwe.obshetperregion.mat.df) # 
dim(nanhwe.obshetperregion.mat.df) #3356   45

#View(nanhwe.obshetperregion.mat.df)

### subset loci to remove - ie anything under 0.05 #####
library(dplyr)


head(nanhwe.obshetperregion.mat.df)

summary(nanhwe.obshetperregion.mat.df)


ponds <- colnames(nanhwe.obshetperregion.mat.df)[ 1 : (ncol(nanhwe.obshetperregion.mat.df)-1) ]


nanhwe.obshetperregion.mat.df.filt.hwe <- nanhwe.obshetperregion.mat.df %>% filter( if_all(ponds, function(x) x >= 0.01|is.na(.)) )

dim(nanhwe.obshetperregion.mat.df.filt.hwe ) # 1434   44

str(nanhwe.obshetperregion.mat.df.filt.hwe)

summary(nanhwe.obshetperregion.mat.df.filt.hwe)

hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe <- rownames(nanhwe.obshetperregion.mat.df.filt.hwe)
str(hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe) # chr [1:1434]


gl.toad.obshetPerregion.nosibs.hwe <- gl.toad.obshetPerregion.nosibs[ , hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe]
gl.toad.obshetPerregion.nosibs.hwe

"/// GENLIGHT OBJECT /////////

 // 510 genotypes,  1,434 binary SNPs, size: 1.5 Mb
 106586 (14.57 %) missing data

 // Basic content
   @gen: list of 510 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  510 individual labels
   @loc.names:  1434 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 34-204)
   @other: a list containing: elements without names "

#################################
###        FIS filter  ###################
#########################



#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs.hwe<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs.hwe

"/// GENIND OBJECT /////////

 // 510 individuals; 1,434 loci; 2,868 alleles; size: 6.4 Mb

 // Basic content
   @tab:  510 x 2868 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2868 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 34-204)"

pop(genid.toad.obshetPerregion.nosibs.hwe)<-pop.data$pop
pop(genid.toad.obshetPerregion.nosibs.hwe)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe)

#hierf.toad.obshetPerregion.nosibs.hwe<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe,pop=TRUE)
#hierf.toad.obshetPerregion.nosibs.hwe$pop<-pop.data$pop

hierf.toad.obshetPerregion.nosibs.hwebasic<-basic.stats(hierf.toad.obshetPerregion.nosibs.hwe)

hierf.toad.obshetPerregion.nosibs.hwebasic$Fis

head(hierf.toad.obshetPerregion.nosibs.hwebasic$Fis)



### add new row names column
#hierf.toad.obshetPerregion.nosibs.hwebasic$SNP<-xxx

FISdataperpond<-hierf.toad.obshetPerregion.nosibs.hwebasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_510INDIV_1434SNPS_obshet06_rmsibs_hwe001_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()



######### filter

# snps I want to keep
nrow(FISdataperpondavg[which(FISdataperpondavg$Average<0.25 & FISdataperpondavg$Average>-0.31),]) 
# 1430

SNPkeepaveFIS<-FISdataperpondavg[which(FISdataperpondavg$Average<0.25 & FISdataperpondavg$Average>-0.31),]

SNPkeep<-SNPkeepaveFIS$SNP
str(SNPkeep) # chr [1:1430]

# edit names

SNPkeep2<- sub( '\\.', ':', SNPkeep)

str(SNPkeep2)

gl.toad.obshetPerregion.nosibs.hwe.FIS <- gl.toad.obshetPerregion.nosibs.hwe[ , SNPkeep2]
gl.toad.obshetPerregion.nosibs.hwe.FIS

"  /// GENLIGHT OBJECT /////////

 // 510 genotypes,  1,430 binary SNPs, size: 1.5 Mb
 106246 (14.57 %) missing data

 // Basic content
   @gen: list of 510 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  510 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 34-204)
   @other: a list containing: elements without names "


#######

#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs.hwe.FIS
"/// GENIND OBJECT /////////

 // 510 individuals; 1,430 loci; 2,860 alleles; size: 6.4 Mb

 // Basic content
   @tab:  510 x 2860 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2860 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 34-204)"

pop(genid.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data$pop
pop(genid.toad.obshetPerregion.nosibs.hwe.FIS)

### check after filter ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe.FIS)

#hierf.toad.obshetPerregion.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe.FIS,pop=TRUE)
#hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop<-pop.data$pop

hierf.toad.obshetPerregion.nosibs.hwe.FISbasic<-basic.stats(hierf.toad.obshetPerregion.nosibs.hwe.FIS)

hierf.toad.obshetPerregion.nosibs.hwe.FISbasic$Fis

head(hierf.toad.obshetPerregion.nosibs.hwe.FISbasic$Fis)



### add new row names column
#hierf.toad.obshetPerregion.nosibs.hwe.FISbasic$SNP<-xxx

FISdataperpond<-hierf.toad.obshetPerregion.nosibs.hwe.FISbasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_510INDIV_1430SNPS_obshet06_rmsibs_hwe001_FIS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_510INDIV_1434SNPS_obshet06_rmsibs_hwe001_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=50))
dev.off()




#######################################################
## SUBSET VCF TO LOCI AND INDIV FILTERED ########
#####################################################

## subset in vcftools and then use excel to subset regions, then upload each region

#VCFtest<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS_editlocinamesinR.vcf")
head(VCFtest@fix[,'ID'])

# REMOVE SNPS TO KEEP THAT HAVE BEEN THROUGH OBS HET, HWE AND FIS FILTER
VCF.FIS.SNPS<-subset(VCF, VCF@fix[,'ID'] %in% SNPkeep2)
VCF.FIS.SNPS

"***** Object of Class vcfR *****
766 samples
29 CHROMs
1,430 variants
Object size: 22.9 Mb
0 percent missing data
*****        *****         *****"

## save VCF

vcfR::write.vcf(VCF.FIS.SNPS, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_1430SNPS_editlocinamesinR.vcf.gz")


### VCFTOOLS ####

## now upload to server along with rm sibs list, and remove sibs using vcftools, then download the new file, upload to here.... 510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools.recode.vcf

# vcftools --vcf 766INDIV_1430SNPS_editlocinamesinR.vcf --remove listindiv_plink_sub_sib_rel_over_04_names_unique.indv --recode --recode-INFO-all --out 510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools




### SUBSET VCF REMOVE INDIV
str(pop.data)

pop.data$sample.id


indiv_samples_keep<-pop.data$sample.id

str(indiv_samples_keep) #  chr [1:510]

#should be 510, not 509
VCF.FIS.SNPS[,indiv_samples_keep]

"***** Object of Class vcfR *****
509 samples
29 CHROMs
1,430 variants
Object size: 15.2 Mb
0 percent missing data
*****        *****         *****"

test<-VCF.FIS.SNPS[,indiv_samples_keep]

"Warning message:
In .local(x, i, j, ..., drop = drop) :
  You have chosen to omit the FORMAT column, this is typically undesireable."


### IT'S CUTTING OUT FORMAT COLUMN SO ADD TO LIST

# add list at top called "FORMAT"
indiv_samples_keep

indiv_samples_keep_FORMAT<- append('FORMAT',indiv_samples_keep)

#WORKED!!!
VCF.FIS.SNPS_510<-VCF.FIS.SNPS[,indiv_samples_keep_FORMAT]

"***** Object of Class vcfR *****
510 samples
29 CHROMs
1,430 variants
Object size: 15.2 Mb
0 percent missing data
*****        *****         *****"




#######################
###### SFS SITE FREQUENCY SPECTRUM ##########
#########################


###########################################
#### SPLIT INTO four REGIONS ####
###########################################


pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.nosibs.hwe.FIS)

"$HaidaGwai
 /// GENLIGHT OBJECT /////////

 // 204 genotypes,  1,430 binary SNPs, size: 672.2 Kb
 43173 (14.8 %) missing data

 // Basic content
   @gen: list of 204 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  204 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 204-204)
   @other: a list containing: elements without names 


$LowerMain
 /// GENLIGHT OBJECT /////////

 // 160 genotypes,  1,430 binary SNPs, size: 540.8 Kb
 30523 (13.34 %) missing data

 // Basic content
   @gen: list of 160 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  160 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 160-160)
   @other: a list containing: elements without names 


$Northwest
 /// GENLIGHT OBJECT /////////

 // 34 genotypes,  1,430 binary SNPs, size: 208.2 Kb
 5394 (11.09 %) missing data

 // Basic content
   @gen: list of 34 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  34 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 34-34)
   @other: a list containing: elements without names 


$VanIsland
 /// GENLIGHT OBJECT /////////

 // 112 genotypes,  1,430 binary SNPs, size: 438.3 Kb
 27156 (16.96 %) missing data

 // Basic content
   @gen: list of 112 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  112 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size ran"



gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland

"   "



gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain
"  
"

gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai
"   "



gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$Northwest
gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest

" "

dartR::gl.report.maf(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, maf.limit = 0.5, ind.limit = 0, loc.limit = 0,
                     v = 2)

### convert genlight to genind objects 
genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.obshetPerregion.nosibs.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)


###########################################
#### convert to structure ####
###########################################

########### genind2structureL ######
#https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R

source("F:/GBS_data_03_02_21/genind2structureL_function.R")



genind2structureL(genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/HG_genind_structure_510NDIV_1348SNPS_regionID.str", pops=TRUE)


genind2structureL(genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/VI_VanIsland_genind_structure_510NDIV_1348SNPS_regionID.str", pops=TRUE)


genind2structureL(genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/ML_LowerMain_genind_structure_510NDIV_1348SNPS_regionID.str", pops=TRUE)


###########################################
#### plot folded MAF - as count of alleles, for each region ####
###########################################

# split VCF based on region



### HG #####
##############################


### SUBSET HAIDA GWAII

gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai@ind.names

HG_data_204_1430_indnames<-gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai@ind.names

HG_data_204_1430_indnames_FORMAT<-append('FORMAT',HG_data_204_1430_indnames)

VCF.FIS.SNPS.HG_204_1430<-VCF.FIS.SNPS[,HG_data_204_1430_indnames_FORMAT]
VCF.FIS.SNPS.HG_204_1430

"***** Object of Class vcfR *****
204 samples
29 CHROMs
1,430 variants
Object size: 6.3 Mb
0 percent missing data
*****        *****         *****"


#minor allele freq
VCF.HG_204_1430_mafvcf<-maf(VCF.FIS.SNPS.HG_204_1430,element=2)

VCF.HG_204_1430_mafvcf.df<-as.data.frame(VCF.HG_204_1430_mafvcf)
summary(VCF.HG_204_1430_mafvcf.df)

head(VCF.HG_204_1430_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_204INDIV_510INDIV_1430SNPS_freq_cutindvR_maf_xlim0_1.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_204_1430_mafvcf.df$Frequency,breaks=seq(0,1,l=50), main="MAF Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_204INDIV_510INDIV_1430SNPS_freq_cutindvR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_204_1430_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main="MAF Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_204INDIV_510INDIV_1430SNPS_COUNT_cutindvR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_204_1430_mafvcf.df$Count, main="MAC Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()



HGcounthist<-hist(VCF.HG_204_1430_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Frequency")






HG_mafvcf.df<-VCF.HG_204_1430_mafvcf.df
hist(HG_mafvcf.df$Count, main="MAC Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")

### VI #####
##############################

## SUBSET VanIsland
gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland

" /// GENLIGHT OBJECT /////////

 // 112 genotypes,  1,430 binary SNPs, size: 438.3 Kb
 27156 (16.96 %) missing data

 // Basic content
   @gen: list of 112 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  112 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 112-112)
   @other: a list containing: elements without names "

gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland@ind.names

VI_data_112_1430_indnames<-gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland@ind.names


VI_data_112_1430_indnames_FORMAT<-append('FORMAT',VI_data_112_1430_indnames)

VCF.FIS.SNPS.VI_112_1430<-VCF.FIS.SNPS[,VI_data_112_1430_indnames_FORMAT]
VCF.FIS.SNPS.VI_112_1430


VI_mafvcf<-maf(VCF.FIS.SNPS.VI_112_1430,element=2)


VI_mafvcf.df<-as.data.frame(VI_mafvcf)

head(VI_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_112INDIV_510NDIV_1430SNPS_freq_cutinR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VI_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" MAF Vancouver Island filtered cut in R",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_112INDIV_510NDIV_1430SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VI_mafvcf.df$Count, main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()


VIcounthist<-hist(VI_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Frequency")

### LM #####
##############################
#loading the data file

## SUBSET LowerMain
gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain

"/// GENLIGHT OBJECT /////////

 // 160 genotypes,  1,430 binary SNPs, size: 540.8 Kb
 30523 (13.34 %) missing data

 // Basic content
   @gen: list of 160 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  160 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each indiLMdual (group size range: 160-160)
   @other: a list containing: elements without names  "

gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain@ind.names

LM_data_160_1430_indnames<-gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain@ind.names


LM_data_160_1430_indnames_FORMAT<-append('FORMAT',LM_data_160_1430_indnames)

VCF.FIS.SNPS.LM_160_1430<-VCF.FIS.SNPS[,LM_data_160_1430_indnames_FORMAT]
VCF.FIS.SNPS.LM_160_1430

"***** Object of Class vcfR *****
160 samples
29 CHROMs
1,430 variants
Object size: 8.5 Mb
0 percent missing data
*****        *****         *****"


LM_mafvcf<-maf(VCF.FIS.SNPS.LM_160_1430,element=2)


LM_mafvcf.df<-as.data.frame(LM_mafvcf)

head(LM_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_160INDIV_510NDIV_1430SNPS_freq_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(LM_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_160INDIV_510NDIV_1430SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(LM_mafvcf.df$Count, main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()



LMcounthist<-hist(LM_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Number of SNPs")



### NW #####
##############################
#loading the data file
# SUBSET Northwest
gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest

" /// GENLIGHT OBJECT /////////

 // 34 genotypes,  1,430 binary SNPs, size: 208.2 Kb
 5394 (11.09 %) missing data

 // Basic content
   @gen: list of 34 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  34 individual labels
   @loc.names:  1430 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 34-34)
   @other: a list containing: elements without names "

gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest@ind.names

NW_data_34_1430_indnames<-gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest@ind.names


NW_data_34_1430_indnames_FORMAT<-append('FORMAT',NW_data_34_1430_indnames)

VCF.FIS.SNPS.NW_34_1430<-VCF.FIS.SNPS[,NW_data_34_1430_indnames_FORMAT]
VCF.FIS.SNPS.NW_34_1430

"***** Object of Class vcfR *****
34 samples
29 CHROMs
1,430 variants
Object size: 3.7 Mb
0 percent missing data
*****        *****         *****"


NW_mafvcf<-maf(VCF.FIS.SNPS.NW_34_1430,element=2)


NW_mafvcf.df<-as.data.frame(NW_mafvcf)

head(NW_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_NW_34INDIV_510NDIV_1430SNPS_freq_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(NW_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_NW_34INDIV_510NDIV_1430SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(NW_mafvcf.df$Count, main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()



NWcounthist<-hist(NW_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Number of SNPs")





######################################
############### make r facet plot of SFS ###########
####################################
HG_mafvcf.df$Region<-rep("HaidaGwai")
VI_mafvcf.df$Region<-rep("VanIsland")
LM_mafvcf.df$Region<-rep("LowerMain")
NW_mafvcf.df$Region<-rep("Northwest")


allfourregions<-rbind(HG_mafvcf.df,VI_mafvcf.df)
allfourregions<-rbind(allfourregions,LM_mafvcf.df)
allfourregions<-rbind(allfourregions,NW_mafvcf.df)


head(allfourregions)
dim(allfourregions) # 5720    5
## loci name vs freq

### maf - freq no bin
labels4clustersgenetics <- c(VanIsland = "Vancouver Island", LowerMain = "Lower mainland", HaidaGwai = "Haida Gwaii", Northwest =  "Northwest BC")


## maf
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_510_1430.png", width = 3.5, height = 6.5, units = 'in', res = 300)
ggplot(data=allfourregions, aes(x=Frequency,fill=Region)) +
  geom_histogram()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                    name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1,labeller = labeller(Region = labels4clustersgenetics))+
  ylab("Number of loci") + xlab("Minor allele frequency")+
  theme_bw()+
  xlim(0,0.5)+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()


### count
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_COUNT_vs_locicount_510INDIV_1430SNPS.png", width = 3.5, height = 6.5, units = 'in', res = 300)
ggplot(data=allfourregions, aes(x=Count,fill=Region)) +
  geom_histogram(bins = 25)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                    name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1,labeller = labeller(Region = labels4clustersgenetics))+
  ylab("Number of loci") + xlab("Minor allele count")+
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()


ggsave("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_510INDIV_1430SNPS.svg", width = 3.5, height = 6.5, units = 'in')
ggplot(data=allfourregions, aes(x=Frequency,fill=Region)) +
  geom_histogram()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                    name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  facet_wrap(~Region,ncol=1,labeller = labeller(Region = labels4clustersgenetics))+
  ylab("Number of loci") + xlab("Minor allele frequency")+
  theme_bw()+
  xlim(0,0.5)+
  theme(legend.position="none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()




#############################
######### PCA FIS ##########
###############################

##per region
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.obshetPerregion.nosibs.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)

#gl.toad.obshetPerregion.nosibs.hwe.FIS.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(1,2)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(1,2)]/sum(my_pca$eig)
#  32.958713  2.966368

##per region
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_510INDIV_1430SNPS_obshet06_rmsibs_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (3.0% explained variance)") + xlab("PC1 (33.0% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=2))
p
dev.off()



##per pop
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$pop
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_perpop_4axes_510INDIV_1430SNPS_obshet06_rmsibs_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
#                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='pop') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + guides(colour=guide_legend(nrow=5))
p
dev.off()



## NEW COLS 
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_MOEpop_4axes_510INDIV_1430SNPS_obshet06_rmsibs_hwe001_FIS_orangegrad_legright2.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()






######################################
########## FST FIS ##############
########################################

genid.toad.obshetPerregion.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs.hwe.FIS

" "

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe.FIS,pop=NULL)





#has pop - as fourclusters
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop

hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- pop.data$fourclusters
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- factor(pop.data$fourclusters)

FST_hierf.toad.obshetPerregion.nosibs.hwe.FIS<-pairwise.WCfst(hierf.toad.obshetPerregion.nosibs.hwe.FIS,diploid=T)
FST_hierf.toad.obshetPerregion.nosibs.hwe.FIS

"       HaidaGwai  LowerMain  Northwest  VanIsland
HaidaGwai        NA 0.49161051 0.66036243 0.53306052
LowerMain 0.4916105         NA 0.09216986 0.05699455
Northwest 0.6603624 0.09216986         NA 0.10754298
VanIsland 0.5330605 0.05699455 0.10754298         NA      "

write.table(FST_hierf.toad.obshetPerregion.nosibs.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_510NDIV_1430SNPS_obshet06_rmsibs_hwe001_FIS.txt",sep = "\t",row.names = TRUE,col.names = TRUE )






#has pop - as four areas
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop

hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- pop.data$fourclusters
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- factor(pop.data$fourclusters)

trial2<-pairwise.WCfst(hierf.toad.obshetPerregion.nosibs.hwe.FIS,diploid=T)
trial2

"                   "

write.table(trial2, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_3clusters_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.txt",sep = "\t",row.names = TRUE,col.names = TRUE )


# as pop
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop

hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- pop.data$pop
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- factor(pop.data$pop)

trial3<-pairwise.WCfst(hierf.toad.obshetPerregion.nosibs.hwe.FIS,diploid=T)
trial3

"               "

write.csv(trial3, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_POP_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_take2.csv")




##### p value pop
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$pop

(Fstp<-stamppFst(gl.toad.obshetPerregion.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_POP_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )



##### p value two areas
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$twoclusters

(Fstp<-stamppFst(gl.toad.obshetPerregion.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_2clusters_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )

##### p value four areas
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters

(Fstp<-stamppFst(gl.toad.obshetPerregion.nosibs.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_3regions_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.csv" )





### FST per loc overall pops ########

hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- pop.data$twoclusters
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop <- factor(pop.data$twoclusters)

basic.nosibs.hwe.FIS<-basic.stats(hierf.toad.obshetPerregion.nosibs.hwe.FIS)

perloc.basic.nosibs.hwe.FIS<-basic.nosibs.hwe.FIS$perloc

par(mar = c(3,3,3,3))

summary(perloc.basic.nosibs.hwe.FIS$Fst)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_2.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perloc.basic.nosibs.hwe.FIS$Fst,breaks=seq(-0.003100,1,l=50),
     main = " ",
     xlab="Pairwise Fst per locus between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()


ggsave("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1.svg", width = 8, height = 6.5, units = 'in')
hist(perloc.basic.nosibs.hwe.FIS$Fst,breaks=seq(-0.003100,1,l=50),
     xlab="Pairwise Fst between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()

######## wc

(FSTperloc<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS,diploid=TRUE))
perlocobject<-FSTperloc$per.loc

perlocobject$FST

summary(perlocobject$FST)

hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2regions_510NDIV_1367SNPS_hwe44pops0.01_FISperpond0.1_wcfunc.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50),
     main = " ",
     xlab="Pairwise Fst per locus between Haida Gwaii and southwest BC",
     ylab="Frequency")
dev.off()



####### FST per loc per region #####

#### genind dartR from genligght ####
genidtoad.nosibs.hwe<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe, probar = FALSE, verbose = NULL)
genidtoad.nosibs.hwe
""

pop(genidtoad.nosibs.hwe)<-pop.data$pop
pop(genidtoad.nosibs.hwe)

### check ###


#### heirfstat conversion from genind ########
hierftoad.nosibs.hwe<-genind2hierfstat(genidtoad.nosibs.hwe)









#################
######### fst per locus ######
#########################

#library(devtools)
#install_github('j-a-thia/genomalicious')

## run on vcftools

fstperloc_swbc_hg_vcftools<-read.csv("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_VCFTOOLS_766_1503_swbc_HG_perloc.csv")

hist(fstperloc_swbc_hg_vcftools$WEIR_AND_COCKERHAM_FST)

summary(fstperloc_swbc_hg_vcftools$WEIR_AND_COCKERHAM_FST)


################################
##### FST within region  -  FST overall per pop FIS ########
#################################

########### four areas #############


pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$HaidaGwai



### convert genlight to genind objects 
gind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


gind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


gind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''


### make new pop files for each dataset
pop.data_VanIsland <- filter(pop.data, 4clusters == "VanIsland")
pop.data_LowerMain <- filter(pop.data, 4clusters == "LowerMain")
pop.data_HaidaGwai <- filter(pop.data, 4clusters == "HaidaGwai")




##### VanIsland #############

# make pop breeding pond
pop(gind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland)<-pop.data_VanIsland$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-genind2hierfstat(gind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland$pop


(VanIslandFST<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland,diploid=TRUE))




##### LowerMain #############


# make pop breeding pond
pop(gind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain)<-pop.data_LowerMain$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-genind2hierfstat(gind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain$pop


(LowerMainFST<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain,diploid=TRUE))



##### HaidaGwai #############


# make pop breeding pond
pop(gind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai)<-pop.data_HaidaGwai$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-genind2hierfstat(gind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai$pop


(HaidaGwaiFST<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai,diploid=TRUE))


######## save values

VanIslandFST1<-VanIslandFST$FST
LowerMainFST1<-LowerMainFST$FST
HaidaGwaiFST1<-HaidaGwaiFST$FST


FSTperpop<-cbind(VanIslandFST1,LowerMainFST1,HaidaGwaiFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_4clusters_766_1503.txt")



###################################
########### two clusters #############
########################################

pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$twoclusters

gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


gl.toad.obshetPerregion.nosibs.hwe.FIS.HG<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$HG
gl.toad.obshetPerregion.nosibs.hwe.FIS.swBC<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$swBC



### convert genlight to genind objects 
gind.toad.obshetPerregion.nosibs.hwe.FIS.HG<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.HG, probar = FALSE, verbose = NULL)


""


gind.toad.obshetPerregion.nosibs.hwe.FIS.swBC<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.swBC, probar = FALSE, verbose = NULL)

""





### make new pop files for each dataset
pop.data_HG <- filter(pop.data, twoclusters == "HG")
pop.data_swBC <- filter(pop.data, twoclusters == "swBC")




##### HG #############

# make pop breeding pond
pop(gind.toad.obshetPerregion.nosibs.hwe.FIS.HG)<-pop.data_HG$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS.HG<-genind2hierfstat(gind.toad.obshetPerregion.nosibs.hwe.FIS.HG,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.nosibs.hwe.FIS.HG$pop


(HGFST<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS.HG,diploid=TRUE))




##### swBC #############


# make pop breeding pond
pop(gind.toad.obshetPerregion.nosibs.hwe.FIS.swBC)<-pop.data_swBC$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS.swBC<-genind2hierfstat(gind.toad.obshetPerregion.nosibs.hwe.FIS.swBC,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.nosibs.hwe.FIS.swBC$pop


(swBCFST<-hierfstat::wc(hierf.toad.obshetPerregion.nosibs.hwe.FIS.swBC,diploid=TRUE))


######## save values

HGFST1<-HGFST$FST
swBCFST1<-swBCFST$FST


FSTperpop<-cbind(HGFST1,swBCFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_twoclusters_766_1503.txt")






#############################
##### missingness #########
#######################################


#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.hwe, probar = FALSE, verbose = NULL)


"/// GENIND OBJECT /////////

 // 510 individuals; 1,430 loci; 2,860 alleles; size: 6.4 Mb

 // Basic content
   @tab:  510 x 2860 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2860 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 1-24)"

#set pop as region
pop(genid.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data$fourclusters


#has pop
genid.toad.obshetPerregion.nosibs.hwe.FIS@pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe.FIS,pop=NULL)



basicstat4clusters<-basic.stats(hierf.toad.obshetPerregion.nosibs.hwe.FIS,diploid=TRUE)
hierf.toad.obshetPerregion.ind<-basicstat4clusters$n.ind.samp


### count number indiv per region

numVanIsland<-sum(pop.data$fourclusters=="VanIsland") # 112
numLowerMain<-sum(pop.data$fourclusters=="LowerMain") # 160
numHaidaGwai<-sum(pop.data$fourclusters=="HaidaGwai") # 204
numNorthwest<-sum(pop.data$fourclusters=="Northwest") # 34


hierf.toad.obshetPerregion.ind




dfhierf.toad.obshetPerregion.ind<-as.data.frame(hierf.toad.obshetPerregion.ind)

dfhierf.toad.obshetPerregion.ind$VanIsland_pct<-dfhierf.toad.obshetPerregion.ind$VanIsland/numVanIsland
dfhierf.toad.obshetPerregion.ind$LowerMain_pct<-dfhierf.toad.obshetPerregion.ind$LowerMain/numLowerMain
dfhierf.toad.obshetPerregion.ind$HaidaGwai_pct<-dfhierf.toad.obshetPerregion.ind$HaidaGwai/numHaidaGwai
dfhierf.toad.obshetPerregion.ind$Northwest_pct<-dfhierf.toad.obshetPerregion.ind$Northwest/numNorthwest



dfhierf.toad.obshetPerregion.ind$VanIsland_pctmissing<-(1-(dfhierf.toad.obshetPerregion.ind$VanIsland/numVanIsland))*100
dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing<-(1-(dfhierf.toad.obshetPerregion.ind$LowerMain/numLowerMain))*100
dfhierf.toad.obshetPerregion.ind$HaidaGwai_pctmissing<-(1-(dfhierf.toad.obshetPerregion.ind$HaidaGwai/numHaidaGwai))*100
dfhierf.toad.obshetPerregion.ind$Northwest_pctmissing<-(1-(dfhierf.toad.obshetPerregion.ind$Northwest/numNorthwest))*100


summary(dfhierf.toad.obshetPerregion.ind)


"VanIsland_pctmissing LowerMain_pctmissing HaidaGwai_pctmissing Northwest_pctmissing 
 Min.   : 0.00        Min.   : 0.00        Min.   : 0.00        Min.   : 0.000             
 1st Qu.:12.50        1st Qu.:10.62        1st Qu.:11.27        1st Qu.: 5.882            
 Median :16.07        Median :13.12        Median :15.20        Median : 8.824            
 Mean   :16.96        Mean   :13.34        Mean   :14.80        Mean   :11.094          
 3rd Qu.:20.54        3rd Qu.:15.62        3rd Qu.:19.12        3rd Qu.:14.706            
 Max.   :42.86        Max.   :37.50        Max.   :34.80        Max.   :58.824           "

head(dfhierf.toad.obshetPerregion.ind)



hist(dfhierf.toad.obshetPerregion.ind$HaidaGwai_pctmissing)

hist(dfhierf.toad.obshetPerregion.ind$VanIsland_pctmissing)

hist(dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing)

hist(dfhierf.toad.obshetPerregion.ind$Northwest_pctmissing)

par(mfrow=c(2,2))

par(mfrow=c(1,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Missingness_VanIsland_LowerMain_HaidaGwai_Northwest_percentmissingloci_xlim_ylim_510INDIV_1430SNPS_long.png", width = 5, height = 10, units = 'in', res = 600)
par(mfrow=c(4,1))
hist(dfhierf.toad.obshetPerregion.ind$HaidaGwai_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$VanIsland_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$Northwest_pctmissing,xlim = c(0,100),ylim = c(0,1600))
dev.off()





png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Missingness_VanIsland_LowerMain_HaidaGwai_Northwest_percentmissingloci_xlim_ylim_510INDIV_1430SNPS.png", width = 20, height = 10, units = 'in', res = 600)
par(mfrow=c(1,4))
hist(dfhierf.toad.obshetPerregion.ind$HaidaGwai_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$VanIsland_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
hist(dfhierf.toad.obshetPerregion.ind$Northwest_pctmissing,xlim = c(0,100),ylim = c(0,1600))
dev.off()




png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/VanIsland_LowerMain_HaidaGwai_percentmissingloci_xlim_ylim_40pct_601INDIV_601SNPS.png", width = 10, height = 6.5, units = 'in', res = 600)
par(mfrow=c(1,4))
hist(dfhierf.toad.obshetPerregion.ind$HaidaGwai_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
hist(dfhierf.toad.obshetPerregion.ind$VanIsland_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
hist(dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
hist(dfhierf.toad.obshetPerregion.ind$LowerMain_pctmissing,xlim = c(0,100),ylim = c(0,1600))
abline(v=40,col="blue")
dev.off()





######################################
########## NJ TREE ##################
###############################################

library(ape)
tre <- nj(dist(as.matrix(gl.toad.obshetPerregion.nosibs.hwe.FIS)))
tre
plot(tre, typ="fan", cex=0.7)
#title("NJ tree of the 18 test samples filter 22")

#par(mar=c(5.1,4.1,4.1,2.1))



#colour 4clusters
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/NJTREE/NJTREE_col4clusters_510NDIV_1430SNPS.png", width = 12, height = 12, units = 'in', res = 300)
plot(tre, typ="fan", show.tip=TRUE, cex=0.5)
tiplabels(pch=20, col=brewer.pal(10, "Paired")[gl.toad.obshetPerregion.nosibs.hwe.FIS$pop], cex=2)
#title("NJ tree of Western toads for the 20 test samples/n colour from subregion/n denovo test 20 samples rm sibs gapped 0.9 M3 rm sibs fitlering may 11th")
dev.off()






###################################
######## DAPC ###############
#########################




# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf


# 1. find correct no. clusters with find.clusters()
# 2. run DAPC for each cluster that looks reasonable - using 200 PCs - using dapc()
# 3. find the optimal number of PCs using optim.a.score()
# 4. membership prob plots

#################################
######### find correct number clusters ##########
####################################

############################
####### 200 PCs, 3 groups

grp <- find.clusters(gl.toad.obshetPerregion.nosibs.hwe.FIS, max.n.clust=26)


# save variance explained PC plot
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_varianceexplainedbyPCA_findclusters_maxclusters26_510NDIV_1430SNPS.png")



# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# pick largest
#Choose the number PCs to retain (>=1): 
200

# save BIC plot

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_BIC_findclusters_maxclusters26_200PCs_510NDIV_1430SNPS.png")



#Choose the number of clusters (>=2): 
3

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_3clusters_510NDIV_1430SNPS.csv")



############################
####### 200 PCs, 2 groups #######

grp <- find.clusters(gl.toad.obshetPerregion.nosibs.hwe.FIS, max.n.clust=26)


#Choose the number PCs to retain (>=1): 
200


#Choose the number of clusters (>=2): 
2

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_2clusters_510NDIV_1430SNPS.csv")



############################
####### 200 PCs, 4 groups #######

grp <- find.clusters(gl.toad.obshetPerregion.nosibs.hwe.FIS, max.n.clust=44)


#Choose the number PCs to retain (>=1): 
350


#Choose the number of clusters (>=2): 
4

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters44_350PCs_4clusters_510NDIV_1430SNPS.csv")

DAPC_varianceexplainedbyPCA_findclusters_maxclusters44_510NDIV_1430SNPS

DAPC_BIC_findclusters_maxclusters44_350PCs_4clusters_510NDIV_1430SNPS

##############################
######## DAPC CLUSTER 4 ###########
###########################



clustass4clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters44_350PCs_4clusters_510NDIV_1430SNPS.csv")



clustass4clusters<-as.data.frame(clustass4clusters)

str(clustass4clusters)

clustass4clusters$sample.id<-clustass4clusters$X

head(clustass4clusters)

clustass4clusters$fourclusters350pca<-clustass4clusters$x

head(clustass4clusters)


#JUST SAVE 2 COLS
clustass4clusterssub<-clustass4clusters[c(3:4)]

head(clustass4clusterssub)




### join w pop data 
pop.data4clusters<-left_join(pop.data,clustass4clusterssub)

head(pop.data4clusters)


# work out which pop corresponds to which cluster number
#View(pop.data4clusters)


pop.data4clusters$fourclusters350pcaNAMES<-pop.data4clusters$fourclusters350pca

pop.data4clusters$fourclusters350pcaNAMES[which(pop.data4clusters$fourclusters350pca=="1")] <- "LowerMain"

pop.data4clusters$fourclusters350pcaNAMES[which(pop.data4clusters$fourclusters350pca=="4")] <- "VanIsland"

pop.data4clusters$fourclusters350pcaNAMES[which(pop.data4clusters$fourclusters350pca=="2")] <- "HaidaGwai"

pop.data4clusters$fourclusters350pcaNAMES[which(pop.data4clusters$fourclusters350pca=="3")] <- "Northwest"

head(pop.data4clusters)



pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data4clusters$fourclusters350pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 350, n.da = 10, parallel = require("parallel"))




#### check optimal number of PCs = 23
temp <- optim.a.score(pnw.dapc)

ascoreopt_PC_28_maxclusters44_350PCs_4clusters_510NDIV_1430SNPS


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 18, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_510NDIV_1430SNPS_optimal_18PCaxes_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()



# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 18, n.da = 4, parallel = require("parallel"))

scatter(pnw.dapc3, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



#p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
#                       name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c#("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))


labs <- c("HaidaGwai", "LowerMain","NorthWest","VanIsland")
cols=c("#009E73","#56B4E9", "pink","#E69F00" )


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_510NDIV_1430SNPS_optimal_18PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_510NDIV_1430SNPS_350pcaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



### DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_510NDIV_1430SNPS_80pctPCs_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_510NDIV_1430SNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()





#Loci contributing to observed differences, threshold set arbitrarily?
set.seed(4)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_loadings_col_4_clusters_510NDIV_1430SNPS_optimal_23PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)
dev.off()




#### MEMBERSHIP PROP ######
#probability of population belonging to the pop it's assigned

set.seed(999)
pramx <- xvalDapc(tab(gl.toad.obshetPerregion.nosibs.hwe.FIS), pop(gl.toad.obshetPerregion.nosibs.hwe.FIS), parallel = "snow")
###--->40



compoplot(pnw.dapc2,col = brewer.pal(4, "Paired"), posi = 'top')


dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")


labs <- c("HaidaGwai", "LowerMain","NorthWest","VanIsland")
cols=c("#009E73","#56B4E9", "pink","#E69F00" )


###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_510NDIV_1430SNPS_optimalPC_18PCaxes_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c( "#009E73","#56B4E9", "pink","#E69F00"),
                           name="Assigned Pop K=4 cluster", breaks = c("HaidaGwai", "LowerMain","NorthWest","VanIsland"), labels=c("HaidaGwai", "LowerMain","NorthWest","VanIsland"))
t
dev.off()





# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_510NDIV_1430SNPS_optimalPC_18PCaxes_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 2,vjust=0.2),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c( "#009E73","#56B4E9", "pink","#E69F00"),
                           name="Assigned Pop K=4 cluster", breaks = c("HaidaGwai", "LowerMain","NorthWest","VanIsland"), labels=c("HaidaGwai", "LowerMain","NorthWest","VanIsland"))
t
dev.off()




t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()



######## as four original pops 

pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data4clusters$fourclusters350pcaNAMES



#make labs
supp.labs <- c(  HaidaGwai=" HaidaGwai (original pop)",LowerMain="LowerMain (original pops)",NorthWest="NorthWest (original pop)",VanIsland ="VanIsland (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_510NDIV_1430SNPS_optimalPC_18PCaxes_4originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <-t + scale_fill_manual(values=c( "#009E73","#56B4E9", "pink","#E69F00"),
                          name="Assigned Pop K=4 cluster", breaks = c("HaidaGwai", "LowerMain","NorthWest","VanIsland"), labels=c("HaidaGwai", "LowerMain","NorthWest","VanIsland"))
t
dev.off()



## as two original pops....
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data4clusters$two_areas






#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Northern_HaidaGwaii="Northern Haida Gwaii (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_510NDIV_1430SNPS_optimalPC_23PCaxes_2originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#58d19f","#D55E00","pink"),
                           name="Assigned Pop K4 cluster", breaks = c("Gwaii_Haanas","Chikundal","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Chikundal","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()






##############################
######## DAPC CLUSTER 3 ###########
###########################



clustass3clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_3clusters_510NDIV_1430SNPS.csv")
dim(clustass3clusters)
dim(pop.data)

clustass3clusters<-as.data.frame(clustass3clusters)

str(clustass3clusters)

clustass3clusters$sample.id<-clustass3clusters$X

head(clustass3clusters)

clustass3clusters$fourclusters350pca<-clustass3clusters$x

#JUST SAVE 2 COLS
clustass3clusterssub<-clustass3clusters[c(3:4)]

head(clustass3clusterssub)




### join w pop data 
pop.data3clusters<-left_join(pop.data,clustass3clusterssub)

head(pop.data3clusters)


# work out which pop corresponds to which cluster number
#View(pop.data3clusters)


pop.data3clusters$fourclusters350pcaNAMES<-pop.data3clusters$fourclusters350pca

pop.data3clusters$fourclusters350pcaNAMES[which(pop.data3clusters$fourclusters350pca=="2")] <- "Haida_Gwaii"

pop.data3clusters$fourclusters350pcaNAMES[which(pop.data3clusters$fourclusters350pca=="1")] <- "SwBC1"

pop.data3clusters$fourclusters350pcaNAMES[which(pop.data3clusters$fourclusters350pca=="3")] <- "SwBC_Sea2Sky"


head(pop.data3clusters)



pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data3clusters$fourclusters350pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 200, n.da = 10, parallel = require("parallel"))




#### check optimal number of PCs = 24
temp <- optim.a.score(pnw.dapc)


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 24, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_510NDIV_1430SNPS_optimal_24PCaxes_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()




# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 16, n.da = 4, parallel = require("parallel"))

scatter(pnw.dapc3, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)






labs <- c( "Gudal Lake","Gwaii Haanas","Northern Haida Gwaii")
cols=c("pink","#0072B2", "#D55E00" )



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_510NDIV_1430SNPS_optimal_24PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_510NDIV_1430SNPS_350pcaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



## #DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_3_clusters_510NDIV_1430SNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()





#Loci contributing to observed differences, threshold set arbitrarily?
set.seed(4)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_loadings_col_3_clusters_510NDIV_1430SNPS_optimal_24PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)
dev.off()



#### MEMBERSHIP PROP ######
#probability of population belonging to the pop it's assigned


pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data3clusters$fourclusters350pcaNAMES

set.seed(999)
pramx <- xvalDapc(tab(gl.toad.obshetPerregion.nosibs.hwe.FIS), pop(gl.toad.obshetPerregion.nosibs.hwe.FIS), parallel = "snow")
###--->40



compoplot(pnw.dapc2,col = brewer.pal(4, "Paired"), posi = 'top')



dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")

###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_510NDIV_1430SNPS_optimalPC_24PCaxes_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()


# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_510NDIV_1430SNPS_optimalPC_24PCaxes_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 3),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c("#009E73", "#9999CC",  "#7949a6"),
                           name="Assigned Pop K=3 cluster", breaks = c("Haida_Gwaii","SwBC1", "SwBC_Sea2Sky"), labels=c("Haida_Gwaii", "SwBC1", "SwBC_Sea2Sky"))
t
dev.off()


######## as four original pops 

pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data4clusters$fourclusters350pcaNAMES



#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Chikundal="Chikundal (original pops)",Northern_Haida_Gwaii="Northern Haida Gwaii (original pop)",Gudal_Lake ="Gudal Lake (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_510NDIV_1430SNPS_optimalPC_24PCaxes_4originalareas.png", width = 12, height = 6.5, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#58d19f","#D55E00","pink"),
                           name="Assigned Pop K4 cluster", breaks = c("Gwaii_Haanas","Chikundal","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Chikundal","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()



## as two original pops....
pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data3clusters$two_areas





#make labs
supp.labs <- c( Gwaii_Haanas=" Gwaii Haanas (original pop)",Northern_HaidaGwaii="Northern Haida Gwaii (original pop)")


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_3_clusters_510NDIV_1430SNPS_optimalPC_24PCaxes_2originalareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + facet_grid(~Original_Pop, scales = "free",  labeller = labeller(Original_Pop = supp.labs))
t <- t + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c( "#0072B2","#D55E00","pink"),
                           name="Assigned Pop K=3 cluster", breaks = c("Gwaii_Haanas","Northern_Haida_Gwaii","Gudal_Lake"), labels=c("Gwaii Haanas","Northern Haida Gwaii","Gudal Lake"))
t
dev.off()









##############################
######## DAPC CLUSTER 2 ###########
###########################



clustass2clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_2clusters_510NDIV_1430SNPS.csv")



clustass2clusters<-as.data.frame(clustass2clusters)

str(clustass2clusters)

clustass2clusters$sample.id<-clustass2clusters$X

head(clustass2clusters)

clustass2clusters$twoclusters350pca<-clustass2clusters$x

#JUST SAVE 2 COLS
clustass2clusterssub<-clustass2clusters[c(3:4)]

head(clustass2clusterssub)




### join w pop data 
pop.data2clusters<-left_join(pop.data,clustass2clusterssub)

head(pop.data2clusters)


# work out which pop corresponds to which cluster number
#View(pop.data2clusters)


pop.data2clusters$twoclusters350pcaNAMES<-pop.data2clusters$twoclusters350pca


pop.data2clusters$twoclusters350pcaNAMES[which(pop.data2clusters$twoclusters350pca=="1")] <- "Southwest_BC"

pop.data2clusters$twoclusters350pcaNAMES[which(pop.data2clusters$twoclusters350pca=="2")] <- "Haida_Gwaii"


head(pop.data2clusters)



pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data2clusters$twoclusters350pcaNAMES




## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 200, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)

#### check optimal number of PCs = 1
temp <- optim.a.score(pnw.dapc)


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 1, n.da = 10, parallel = require("parallel"))

scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)

# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, n.pca = 1, n.da = 4, parallel = require("parallel"))


labs <- c("Haida Gwaii","SouthWest BC")
cols=c("#009E73","#E69F00" )



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_510NDIV_1430SNPSSNPS_optimal_1PCaxis.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_510NDIV_1430SNPSSNPS_350pcaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



## #DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.obshetPerregion.nosibs.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_2_clusters_510NDIV_1430SNPSSNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc,col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



########## membership prob

dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.obshetPerregion.nosibs.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")



###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_2_clusters_510NDIV_1430SNPS_optimalPC_1PCaxis_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c("#009E73","#E69F00" ),
                           name="Assigned Pop K=2 cluster", breaks = c("Haida_Gwaii","Southwest_BC"), labels=c("Haida_Gwaii","Southwest_BC"))
t
dev.off()


# with sample id
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_2_clusters_510NDIV_1430SNPS_optimalPC_1PCaxis_nosepareas_sampleid.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_text(angle = 90, size = 3),axis.title.y = element_text( size = 8))
t <- t + scale_fill_manual(values=c("#009E73","#E69F00" ),
                           name="Assigned Pop K=2 cluster", breaks = c("Haida_Gwaii","Southwest_BC"), labels=c("Haida_Gwaii","Southwest_BC"))
t
dev.off()













########################## BASIC STATS FIS ###################

#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toad.obshetPerregion.nosibs.hwe.FIS)<-pop.data$fourclusters


#has pop
genid.toad.obshetPerregion.nosibs.hwe.FIS@pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs.hwe.FIS,pop=NULL)

#has pop - as region
hierf.toad.obshetPerregion.nosibs.hwe.FIS$pop



dataIN<-hierf.toad.obshetPerregion.nosibs.hwe.FIS
pop.stats <- function(dataIN = dataIN, boots = 10000){
  
  # dataIN = population genomic data in hierfstat format
  
  #- calculate stats and fis bootstrapping
  tmp.stats <- hierfstat::basic.stats( dataIN )
  tmp.boot <- hierfstat::boot.ppfis( dataIN, nboot = boots )
  
  #- create output df
  
  tmp.summary <- matrix( nrow = 8, ncol = length( levels(dataIN$pop)),
                         dimnames = list( c( "Hs_mean", "Hs_SD", "Ho_mean", "Ho_SD", "P", "Fis", "Fis_ll", "Fis_hl" ),
                                          levels(dataIN$pop) )
  )
  
  #- calculate mean and sd for Hs and Ho from individual loci estimates
  
  tmp.summary[ "Hs_mean", ] <- apply( tmp.stats$Hs, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Hs_SD", ] <- apply( tmp.stats$Hs, 2, function(x) sd( x, na.rm = T) )
  
  tmp.summary[ "Ho_mean", ] <- apply( tmp.stats$Ho, 2, function(x) mean( x, na.rm = T ) )
  tmp.summary[ "Ho_SD", ] <- apply( tmp.stats$Ho, 2, function(x) sd( x, na.rm = T) )
  
  #- calculate percentage polymorphic loci from 'pop.freq' data
  
  tmp.freq <- do.call( rbind, lapply( tmp.stats$pop.freq, function(x) x[1,] ))
  tmp.P <- apply( tmp.freq, 2, function(x) 1 - ( length( x[ x == 0 | x == 1 ] ) / length(x) ) )
  tmp.summary[ "P", ] <- round( tmp.P*100, 1)
  
  #- calculate per population Fis as per hierfstat's formula: Fis = 1 - Ho/Hs.
  #- get 95% CI for Fist from bootstrap calculations
  
  tmp.summary[ "Fis", ] <- 1 - ( tmp.summary[ 'Ho_mean', ] / tmp.summary[ 'Hs_mean', ] )
  tmp.summary[ "Fis_ll", ] <- tmp.boot$fis.ci$ll
  tmp.summary[ "Fis_hl", ] <- tmp.boot$fis.ci$hl
  
  #- add overall stats
  
  tmp.overall <- c( tmp.stats$overall[[ 'Hs' ]],
                    sd( tmp.stats$Hs, na.rm = T ),
                    tmp.stats$overall[[ 'Ho' ]],
                    sd( tmp.stats$Ho, na.rm = T),
                    NA, # polymorphic loci
                    tmp.stats$overall[[ 'Fis' ]],
                    NA, NA) # overall Fis 95% CI
  
  tmp.summary <- cbind( tmp.summary, "OVERALL" = tmp.overall )
  
  return( as.data.frame( t( round( tmp.summary, 4 ) )))
  
}

forb2.pop.stats.hwe <- pop.stats(hierf.toad.obshetPerregion.nosibs.hwe.FIS)
forb2.pop.stats.hwe

### filter max het!! no negative FIS

"         Hs_mean  Hs_SD Ho_mean  Ho_SD    P     Fis  Fis_ll  Fis_hl
HaidaGwai  0.0080 0.0373  0.0086 0.0417 21.5 -0.0656 -0.1026 -0.0250
LowerMain  0.1333 0.0902  0.1343 0.0957 97.5 -0.0073 -0.0148  0.0002
Northwest  0.0980 0.1402  0.1032 0.1509 49.2 -0.0533 -0.0708 -0.0355
VanIsland  0.1326 0.1031  0.1329 0.1061 91.0 -0.0022 -0.0117  0.0073
OVERALL    0.0929 0.1120  0.0947 0.1177   NA -0.0195      NA      NA "

### pop stats output #####

write.table(forb2.pop.stats.hwe, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/hierfstat_popstats_region_510INDIV_1430SNPS.txt",sep = "\t",row.names = TRUE,col.names = TRUE)


########################################
############ Fis per region CI ################
#########################################
str(forb2.pop.stats.hwe)

forb2.pop.stats.hwe$Fis

Fis_95CI_region<-cbind(forb2.pop.stats.hwe$Fis,forb2.pop.stats.hwe$Fis_ll,forb2.pop.stats.hwe$Fis_hl)
dfFis_95CI_region<-as.data.frame(Fis_95CI_region)

#delete last row
dfFis_95CI_region<-dfFis_95CI_region[-c(4), ] 


## add region column
dfFis_95CI_region$Region<-c("Vancouver Island","Lower Mainland", "Haida Gwaii", "Northwest BC")

#rename column
colnames(dfFis_95CI_region) <- c("Fis","Fis_ll","Fis_hl","Region")    # Applying colnames

## reorder
dfFis_95CI_region$Region<- factor(dfFis_95CI_region$Region, levels = c("Vancouver Island", "Lower Mainland", "Haida Gwaii", "Northwest BC"))
dfFis_95CI_region<- droplevels(dfFis_95CI_region)

## plot - no legend
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p


## formatted correctly
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_510INDIV_1430SNPS_symbolsforamtted.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii","Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  ( ",italic(F)[is],paste(" )  ±  95% CIs"))))) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()


### save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_510INDIV_1430SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii","Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()





##  OLD with FIS all italic
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  ( ",italic(F[i][s]),paste("  ±  95% CIs)"))))) +
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p





## new COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_510INDIV_1430SNPS_neworange.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#ab4e03", "#fca45d","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding  (",italic(F)[plain(IS)],paste(")")))))+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()



####################
#### MAF VS Ho ######
#######################


pop(genid.toad.nolas)<-pop.data$onecluster

#### heirfstat conversion from genind ########
hierf.toad<-genind2hierfstat(genid.toad.nolas,pop=NULL)

#has pop
hierf.toad$pop


hierf.toad.basic<-basic.stats(hierf.toad)
hierf.toad.basic.perloc<-as.data.frame(hierf.toad.basic$perloc)

hierf.toad.basic.perloc$loci <- rownames(hierf.toad.basic.perloc)
head(hierf.toad.basic.perloc)
#rownames(hierf.toad.basic.perloc) <- NULL

#edit loci names to match how formatted in maf
hierf.toad.basic.perloc$loci <- sub( 'X', '', sub( '\\.', '\\/', sub( '\\.', '-', sub( '\\.', ':', hierf.toad.basic.perloc$loci))))
head(hierf.toad.basic.perloc)



#minor allele freq
mafvcf<-maf(filtered.VCF,element=2)



mafvcf.df<-as.data.frame(mafvcf)

head(mafvcf.df)
dim()
mafvcf.df$loci <- rownames(mafvcf.df)

## cind

maf_Ho<-cbind(hierf.toad.basic.perloc,mafvcf.df)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAF_vs_Ho_plot_371INDIV_2392SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
plot(maf_Ho$Frequency,maf_Ho$Ho, xlab="Minor allele frequency",
     ylab="Observed heterozygosity (Ho)")
dev.off()




####################################
######### FIS PER LOCUS ############
#########################

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fis_perlocus_plot_371INDIV_2392SNPS_title.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(maf_Ho$Fis,breaks=seq(-1,1,l=50), xlab="Fis",
     ylab="Frequency",main ="Fis per locus no max obs het filter")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fis_perlocus_plot_371INDIV_2392SNPS.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(maf_Ho$Fis,breaks=seq(-1,1,l=50), xlab="Fis",
     ylab="Frequency",main ="")
dev.off()

#count nymber loci that have FIS under 0
sum(maf_Ho$Fis<0) #177

177/602
[1] 0.8428571

# 84% of loci have Fis below zero




###### what if I split into regions before calc across loci? 




#############################
##### obs het vs missingness #########
#######################################
hierf.toadobs<-basic.stats(hierf.toad)


hierf.toadobs$Ho


















##################################
###########  Ext Het & Obs het FIS ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$VanIsland


" /// GENLIGHT OBJECT /////////

 // 33 genotypes,  732 binary SNPs, size: 286.8 Kb
 5026 (20.81 %) missing data

 // Basic content
   @gen: list of 33 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  33 individual labels
   @loc.names:  732 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 33-33)
   @other: a list containing: elements without names "



gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$LowerMain

"/// GENLIGHT OBJECT /////////

 // 71 genotypes,  732 binary SNPs, size: 374.8 Kb
 10519 (20.24 %) missing data

 // Basic content
   @gen: list of 71 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  71 individual labels
   @loc.names:  732 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 71-71)
   @other: a list containing: elements without names 
"

gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$HaidaGwai

" /// GENLIGHT OBJECT /////////

 // 267 genotypes,  732 binary SNPs, size: 860.7 Kb
 46991 (24.04 %) missing data

 // Basic content
   @gen: list of 267 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  267 individual labels
   @loc.names:  732 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 267-267)
   @other: a list containing: elements without names "

### convert genlight to genind objects 
gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''


############### summary on each one
gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland_sum<-adegenet::summary(gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland)

gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain_sum<-adegenet::summary(gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain)

gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai_sum<-adegenet::summary(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai)


## rename
divVanIsland<-gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland_sum
divLowerMain<-gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain_sum
divHaidaGwai<-gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai_sum

###############################################
############  Obs Het per region ############
##########################################



hist(divVanIsland$Hobs)

VanIslandHobs<-divVanIsland$Hobs
LowerMainHobs<-divLowerMain$Hobs
HaidaGwaiHobs<-divHaidaGwai$Hobs

str(HaidaGwaiHobs)

lsNames=c("VanIslandHobs","LowerMainHobs","HaidaGwaiHobs")


## creates list for each locus
do.call(mapply, c(FUN=c, sapply(lsNames, as.symbol), SIMPLIFY=FALSE))

""


### just join as columns - cbind
HaidaGwaiHobs<-as.data.frame(HaidaGwaiHobs)
LowerMainHobs<-as.data.frame(LowerMainHobs)
VanIslandHobs<-as.data.frame(VanIslandHobs)

Hobsdf <- data.frame(cbind(LowerMainHobs,VanIslandHobs,HaidaGwaiHobs))

head(Hobsdf)

"  "


### turn wide to long
df_longHobs <- tidyr::gather(Hobsdf,
                             key = Region,
                             value = Obs_Het,
                             LowerMainHobs ,  VanIslandHobs ,  HaidaGwaiHobs)

head(df_longHobs)




### plot ###

ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  ylim(0,1)

## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_obs_Het_per_region_510INDIV_1430SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_obs_Het_per_region_510INDIV_1430SNPS_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_longHobs, aes(x=Region, y=Obs_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  xlab("Region") + 
  ylab("Observed heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))
#ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()




##################################
###########  Ext Het & Obs het FIS ####################
############################


#### split matrix into 3 regions ###########

pop(gl.toad.obshetPerregion.nosibs.hwe.FIS) <- pop.data$moeregion

gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.nosibs.hwe.FIS)


gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland

"   "



gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain
"  
"

gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai
"   "



gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest<-gl.toad.obshetPerregion.nosibs.hwe.FIS.bypop$Northwest
gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest

" "

dartR::gl.report.maf(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, maf.limit = 0.5, ind.limit = 0, loc.limit = 0,
                     v = 2)

### convert genlight to genind objects 
genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.obshetPerregion.nosibs.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)

############### summary on each one
genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland_sum<-adegenet::summary(genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland)

genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain_sum<-adegenet::summary(genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain)

genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai_sum<-adegenet::summary(genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai)


genind.toad.obshetPerregion.nosibs.hwe.FIS.Northwest_sum<-adegenet::summary(genind.toad.obshetPerregion.nosibs.hwe.FIS.Northwest)


## rename
divVanIsland<-genind.toad.obshetPerregion.nosibs.hwe.FIS.VanIsland_sum
divLowerMain<-genind.toad.obshetPerregion.nosibs.hwe.FIS.LowerMain_sum
divHaidaGwai<-genind.toad.obshetPerregion.nosibs.hwe.FIS.HaidaGwai_sum
divNorthwest<-genind.toad.obshetPerregion.nosibs.hwe.FIS.Northwest_sum


##################################
###########  Ext Het FIS filter ####################
############################



###### Hexp per region ############

divLowerMain$LowerMain_Hexp<-divLowerMain$Hexp
divVanIsland$VanIsland_Hexp<-divVanIsland$Hexp
divHaidaGwai$HaidaGwai_Hexp<-divHaidaGwai$Hexp
divNorthwest$Northwest_Hexp<-divNorthwest$Hexp


LowerMain_Hexp<-divLowerMain$Hexp
LowerMain_Hexpdf<-as.data.frame(LowerMain_Hexp)

VanIsland_Hexp<-divVanIsland$Hexp
VanIsland_Hexpdf<-as.data.frame(VanIsland_Hexp)

HaidaGwai_Hexp<-divHaidaGwai$Hexp
HaidaGwai_Hexpdf<-as.data.frame(HaidaGwai_Hexp)

Northwest_Hexp<-divNorthwest$Hexp
Northwest_Hexpdf<-as.data.frame(Northwest_Hexp)


divallfour<-cbind(LowerMain_Hexp,VanIsland_Hexp,HaidaGwai_Hexpdf, Northwest_Hexpdf)
head(divallfour)

## rename column names
names(divallfour)
colnames(divallfour) <- c("LowerMain_Hexp", "VanIsland_Hexp","HaidaGwai_Hexp","Northwest_Hexp")    # Applying colnames
head(divallfour)
str(divallfour)

# turn to dataframe
dfdivallfour<-as.data.frame(divallfour)
head(dfdivallfour)


### Hexp means + stdev ######
VanIsland_mean<-mean(dfdivallfour$VanIsland_Hexp) #
VanIsland_sd<-sd(dfdivallfour$VanIsland_Hexp) #

LowerMain_mean<-mean(dfdivallfour$LowerMain_Hexp) #
LowerMain_sd<-sd(dfdivallfour$LowerMain_Hexp) #

HaidaGwai_mean<-mean(dfdivallfour$HaidaGwai_Hexp) #
HaidaGwai_sd<-sd(dfdivallfour$HaidaGwai_Hexp) 

Northwest_mean<-mean(dfdivallfour$Northwest_Hexp) #
Northwest_sd<-sd(dfdivallfour$Northwest_Hexp) 



### save output
Hetsum<-cbind(VanIsland_mean, VanIsland_sd, LowerMain_mean,LowerMain_sd,HaidaGwai_mean,HaidaGwai_sd)

write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_3regions_510INDIV_1430SNPS.txt")





### turn wide to long
df_long <- tidyr::gather(dfdivallfour,
                         key = Region,
                         value = Expect_Het,
                         LowerMain_Hexp ,  VanIsland_Hexp ,  HaidaGwai_Hexp,Northwest_Hexp)

head(df_long)




p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 


### plot ###

p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p


## save
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_510INDIV_1430SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()

## save no ylim
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_510INDIV_1430SNPS_hwe26pops_FISperpond0.1_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(df_long, aes(x=Region, y=Expect_Het, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))+
  xlab("Region") + 
  ylab("Expected heterozygosity")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii"))
#ylim(0,1)
p<-p+theme(legend.position = "none")
p
dev.off()





# with SE and 95% CI ###########


exphet_SE<-group_by(df_long, Region) %>%
  summarise_each(funs(mean=mean(Expect_Het),n=n(),sd=sd(Expect_Het),se=sd(.)/sqrt(n())),Expect_Het)


exphet_SEdf<-as.data.frame(exphet_SE)

exphet_SEdf$Exphet_highCI<-exphet_SEdf$mean+1.96*(exphet_SEdf$sd/sqrt(exphet_SEdf$n))
exphet_SEdf$Exphet_lowCI<-exphet_SEdf$mean-1.96*(exphet_SEdf$sd/sqrt(exphet_SEdf$n))

exphet_SEdf$Exphet_highSE<-exphet_SEdf$mean+exphet_SEdf$se

exphet_SEdf$Exphet_lowSE<-exphet_SEdf$mean-exphet_SEdf$se

exphet_SEdf



## reorder
exphet_SEdf$Region<- factor(exphet_SEdf$Region, levels = c("VanIsland_Hexp", "LowerMain_Hexp", "HaidaGwai_Hexp","Northwest_Hexp"))
exphet_SEdf<- droplevels(exphet_SEdf)

write.csv(exphet_SEdf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_mean_sd_SE_CI_MOE_510INDIV_1430SNPS.csv")


ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_95pctCI_4clusters_510INDIV_1430SNPS_correctCI.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(exphet_SEdf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Exphet_highCI, ymin = Exphet_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Expected Heterozygosity")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()


### NEW COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_95pctCI_4clusters_510INDIV_1430SNPS_correctCI.pdf",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(exphet_SEdf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Exphet_highCI, ymin = Exphet_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Expected Heterozygosity")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()


##########################
#### MAF PLOTS IANS CODE ######
##################################


VCF.FIS.SNPS.sibs<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools.recode.vcf")


# EXPORT GENOTYPE TABLE
GENOTYPE_VCF.FIS.SNPS.sibs<-extract.gt(VCF.FIS.SNPS.sibs)
View(GENOTYPE_VCF.FIS.SNPS.sibs)



VCF

GENOTYPE_VCF<-extract.gt(VCF)
View(GENOTYPE_VCF)


VCF_original

GENOTYPE_VCF_original<-extract.gt(VCF_original)
View(GENOTYPE_VCF_original)

#### plot MAF

#I have altered this so it will just calculate MAF, and following this code, output a dataframe of MAF per SNP
#I don't know if this is the most efficient way, but it should work!
#Genotable is your genotype table formatted in the same way described above
MAFCalc = function(x){ # x is a vector or genotypes
  y = x[x != '00'] # Extract a vector of genotypes and remove missing genotypes
  # Split the genotypes to a list, unlist them to characters, summarise as table with two elements,
  #  sort and extract the first element
  z = sort(table(unlist(strsplit(y, split = ''))))[1]
  #set each column equal to the frequency
  #here the frequency is calculated using only non-missing individuals (i.e. length of y rather than length of x)
  #this is how Ian calculated and appears to be how genind calculates it as well 
  x=rep((z/(2*length(y))), length(x))}



VCF.FIS.SNPS.sibs

"***** Object of Class vcfR *****
510 samples
29 CHROMs
1,430 variants
Object size: 15.2 Mb
0 percent missing data
*****        *****         *****"

GENOTYPE_VCF.FIS.SNPS.sibs<-extract.gt(VCF.FIS.SNPS.sibs)
GENOTYPE_VCF.FIS.SNPS.sibs
dim(GENOTYPE_VCF.FIS.SNPS.sibs) # [1] 1430  510

# snps as rows, indiv as cols data is whether minor or major allele (0 vs 1 ?
Genotable<-GENOTYPE_VCF.FIS.SNPS.sibs

View(Genotable)

# snps as rows, indv as cols, data is MAF
GenotableMAF <- apply(Genotable, MARGIN=2, MAFCalc)
View(GenotableMAF)
dim(GenotableMAF)

row.names(GenotableMAF)<-row.names(Genotable)

View(GenotableMAF)

# puts snps as columns, and indiv as rows - we dont' want that
Genotablet_MAF <- as.data.frame(t(GenotableMAF))
View(Genotablet_MAF)


#just need one column to get the frequency (because the value should be the same across columns)
MAF <- data.frame("SNP"=row.names(GenotableMAF), "MAF"=GenotableMAF[,1])
View(MAF)
dim(MAF) # 1430    2


#now you can run a histogram on the MAF


hist(MAF$MAF,breaks=seq(0,1,l=50))

hist(MAF$MAF,breaks=seq(0,1,l=50), main="MAF whole dataset from vcf genotype file",
     xlab="Minor allele frequency",
     ylab="Frequency")



# try converting genotype table to genind
df <- data.frame(locusA=c("11","11","12","32"),
                 locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
row.names(df) <- .genlab("genotype",4)
df



genind2df(obj)
genind2df(obj, sep="/")


#test it against the adegenet calculation 
## COMPUTE ALLELE FREQUENCIES

genid.toad.obshetPerregion.nosibs.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs.hwe.FIS
"/// GENIND OBJECT /////////

 // 510 individuals; 1,430 loci; 2,860 alleles; size: 6.4 Mb

 // Basic content
   @tab:  510 x 2860 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2860 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 34-204)"

x <-genid.toad.obshetPerregion.nosibs.hwe.FIS
apply(tab(x, freq=TRUE),2,mean, na.rm=TRUE)


## GET MINOR ALLELE FREQUENCY
m.freq <- as.data.frame(minorAllele(x))

head(m.freq)

hist(m.freq$`minorAllele(x)`)


hist(m.freq$`minorAllele(x)`,breaks=seq(0,1,l=50), main="MAF whole dataset from genind file",
     xlab="Minor allele frequency",
     ylab="Frequency")


# don't take out NAs
x<-genid.toad.obshetPerregion.nosibs.hwe.FIS

## GET MINOR ALLELE FREQUENCY
m.freq <- as.data.frame(minorAllele(x))

head(m.freq)
dim(m.freq)

hist(m.freq$`minorAllele(x)`)


hist(m.freq$`minorAllele(x)`,breaks=seq(0,1,l=50), main="MAF whole dataset from genind file no NA filter",
     xlab="Minor allele frequency",
     ylab="Frequency")


#######################################
####################################
###### Rm sibs 0.35 threshold ########
########################################

### Rm sibs

### try 0.35 threshold

sub_sibrelover035_names_unique.char

listsibs035<-sub_sibrelover035_names_unique.char
str(sub_sibrelover035_names_unique.char)

# function alba made to remove siblings
source("F:/GBS_data_03_02_21/Function_only_remove_siblings_genlight_Alba_nov_8th_2021.R")

# run function
gl.toad.obshetPerregion.nosibs035<- remove_sibs_genlight(genlight_object = gl.toad.obshetPerregion,
                                                         names = listsibs035)
gl.toad.obshetPerregion.nosibs035

" /// GENLIGHT OBJECT /////////

 // 289 genotypes,  3,356 binary SNPs, size: 1.8 Mb
 162867 (16.79 %) missing data

 // Basic content
   @gen: list of 289 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  289 individual labels
   @loc.names:  3356 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 31-122)
   @other: a list containing: elements without names "

### edit pop file

pop.data_289<-pop.data.766[which(!pop.data.766$sample.id %in% listsibs035),]
dim(pop.data_289) # 510   7

write.table(pop.data_289,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/popmap.289_samples_header_region_rmINdiv_rmsibs_sorted.txt",quote=FALSE,sep = "\t",row.names=FALSE)

# count pops and get list
str(pop.data_289$pop)

pop.data_289_2<-pop.data_289

pop.data_289_2$pop<-as.factor(pop.data_289_2$pop)

levels(pop.data_289_2$pop)

ponds_kept_289INDIV_3356SNPS<-levels(pop.data_289_2$pop)
ponds_kept_289INDIV_3356SNPS



#######################
## HWE filtering per pond p=0.01 ####
######################
#https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html


## check new pop.data file
dim(pop.data_289) #289

pop(gl.toad.obshetPerregion.nosibs035)

#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs035<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs035, probar = FALSE, verbose = NULL)




library("pegas")
#(nanhwe.full <- hw.test(genid.toad.obshetPerregion.nosibs035, B = 1000)) # performs 1000 permuatations

## hwe per pop

pop(genid.toad.obshetPerregion.nosibs035)<-pop.data_289$pop

### check how many pops
pop(genid.toad.obshetPerregion.nosibs035) # 44

#run hwe test for each pop seprately - and only focus on the analytical P value (B=0)
(nanhwe.obshetPerregion.pop <- seppop(genid.toad.obshetPerregion.nosibs035) %>% lapply(hw.test, B = 0))

# Take the third column with all rows - p value column
(nanhwe.obshetPerregion.mat <- sapply(nanhwe.obshetPerregion.pop, "[", i = TRUE, j = 3)) 


#turn to dataframe
nanhwe.obshetperregion.mat.df<-as.data.frame(nanhwe.obshetPerregion.mat)

str(nanhwe.obshetperregion.mat.df) # 
dim(nanhwe.obshetperregion.mat.df) #3356   44

#View(nanhwe.obshetperregion.mat.df)

### subset loci to remove - ie anything under 0.05 #####
library(dplyr)


head(nanhwe.obshetperregion.mat.df)

summary(nanhwe.obshetperregion.mat.df)


ponds <- colnames(nanhwe.obshetperregion.mat.df)[ 1 : (ncol(nanhwe.obshetperregion.mat.df)-1) ]


nanhwe.obshetperregion.mat.df.filt.hwe <- nanhwe.obshetperregion.mat.df %>% filter( if_all(ponds, function(x) x >= 0.01|is.na(.)) )

dim(nanhwe.obshetperregion.mat.df.filt.hwe ) # 801   44

str(nanhwe.obshetperregion.mat.df.filt.hwe)

summary(nanhwe.obshetperregion.mat.df.filt.hwe)

hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe <- rownames(nanhwe.obshetperregion.mat.df.filt.hwe)
str(hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe) # chr [1:801]


gl.toad.obshetPerregion.nosibs035.hwe <- gl.toad.obshetPerregion.nosibs035[ , hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe]
gl.toad.obshetPerregion.nosibs035.hwe

" /// GENLIGHT OBJECT /////////

 // 289 genotypes,  801 binary SNPs, size: 723.5 Kb
 36559 (15.79 %) missing data

 // Basic content
   @gen: list of 289 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  289 individual labels
   @loc.names:  801 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 31-122)
   @other: a list containing: elements without names "

#################################
###        FIS filter  ###################
#########################



#### genind dartR from genligght ####
genid.toad.obshetPerregion.nosibs035.hwe<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs035.hwe, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs035.hwe


pop(genid.toad.obshetPerregion.nosibs035.hwe)<-pop.data_289$pop
pop(genid.toad.obshetPerregion.nosibs035.hwe)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs035.hwe<-genind2hierfstat(genid.toad.obshetPerregion.nosibs035.hwe)

#hierf.toad.obshetPerregion.nosibs035.hwe<-genind2hierfstat(genid.toad.obshetPerregion.nosibs035.hwe,pop=TRUE)
#hierf.toad.obshetPerregion.nosibs035.hwe$pop<-pop.data_289$pop

hierf.toad.obshetPerregion.nosibs035.hwebasic<-basic.stats(hierf.toad.obshetPerregion.nosibs035.hwe)

hierf.toad.obshetPerregion.nosibs035.hwebasic$Fis

head(hierf.toad.obshetPerregion.nosibs035.hwebasic$Fis)



### add new row names column
#hierf.toad.obshetPerregion.nosibs035.hwebasic$SNP<-xxx

FISdataperpond_035<-hierf.toad.obshetPerregion.nosibs035.hwebasic$Fis

FISdataperpond_035avg <- data.frame("SNP"=rownames(FISdataperpond_035), "Average"=rowMeans(FISdataperpond_035, na.rm=T))
hist(FISdataperpond_035avg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_289INDIV_801SNPS_obshet06_rmsibs_hwe001_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpond_035avg$Average,breaks=seq(-1,1,l=100))
dev.off()



######### filter

# snps I want to keep
nrow(FISdataperpond_035avg[which(FISdataperpond_035avg$Average<0.25 & FISdataperpond_035avg$Average>-0.31),]) 
# 789

SNPkeepaveFIS_035<-FISdataperpond_035avg[which(FISdataperpond_035avg$Average<0.25 & FISdataperpond_035avg$Average>-0.31),]

SNPkeep_035<-SNPkeepaveFIS_035$SNP
str(SNPkeep_035) # chr [1:789]

# edit names

SNPkeep2_035<- sub( '\\.', ':', SNPkeep_035)

str(SNPkeep2_035)

gl.toad.obshetPerregion.nosibs035.hwe.FIS <- gl.toad.obshetPerregion.nosibs035.hwe[ , SNPkeep2_035]
gl.toad.obshetPerregion.nosibs035.hwe.FIS

"/// GENLIGHT OBJECT /////////

 // 289 genotypes,  789 binary SNPs, size: 720.7 Kb
 36088 (15.83 %) missing data

 // Basic content
   @gen: list of 289 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  289 individual labels
   @loc.names:  789 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 31-122)
   @other: a list containing: elements without names "


#############################
######### PCA FIS ##########
###############################

##per region
pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS) <- pop.data_289$fourclusters


my_pca <- glPca(gl.toad.obshetPerregion.nosibs035.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS)

#gl.toad.obshetPerregion.nosibs035.hwe.FIS.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(1,2)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(1,2)]/sum(my_pca$eig)
#  12.459151  4.679972

##per region
pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS) <- pop.data_289$fourclusters
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_289INDIV_789SNPS_obshet06_rmsibs_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (4.7% explained variance)") + xlab("PC1 (12.5% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=2))
p
dev.off()



##per pop
pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS) <- pop.data_289$pop
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.nosibs035.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_perpop_4axes_289INDIV_789SNPS_obshet06_rmsibs_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
#                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='pop') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + guides(colour=guide_legend(nrow=5))
p
dev.off()



## NEW COLS 
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_MOEpop_4axes_289INDIV_789SNPS_obshet06_rmsibs_hwe001_FIS_orangegrad_legright2.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73"),
                            name="Region", breaks = c("VanIsland","LowerMainland", "HaidaGwaii","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='moe region') 
p<-p + ylab("PC2 (2.71% explained variance)") + xlab("PC1 (31.1% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()






######################################
########## FST FIS ##############
########################################

genid.toad.obshetPerregion.nosibs035.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.nosibs035.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.nosibs035.hwe.FIS

"/// GENIND OBJECT /////////

 // 289 individuals; 789 loci; 1,578 alleles; size: 2.2 Mb

 // Basic content
   @tab:  289 x 1578 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 1578 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 31-122) "

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.nosibs035.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.nosibs035.hwe.FIS,pop=NULL)





#has pop - as fourclusters
hierf.toad.obshetPerregion.nosibs035.hwe.FIS$pop

hierf.toad.obshetPerregion.nosibs035.hwe.FIS$pop <- pop.data_289$fourclusters
hierf.toad.obshetPerregion.nosibs035.hwe.FIS$pop <- factor(pop.data_289$fourclusters)

FST_hierf.toad.obshetPerregion.nosibs035.hwe.FIS<-pairwise.WCfst(hierf.toad.obshetPerregion.nosibs035.hwe.FIS,diploid=T)
FST_hierf.toad.obshetPerregion.nosibs035.hwe.FIS

"         HaidaGwai  LowerMain  Northwest VanIsland
HaidaGwai        NA 0.32766600 0.42188478 0.3246088
LowerMain 0.3276660         NA 0.08152651 0.0594244
Northwest 0.4218848 0.08152651         NA 0.0999390
VanIsland 0.3246088 0.05942440 0.09993900        NA       "

write.table(FST_hierf.toad.obshetPerregion.nosibs035.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_289NDIV_789SNPS_obshet06_rmsibs_hwe001_FIS.txt",sep = "\t",row.names = TRUE,col.names = TRUE )







