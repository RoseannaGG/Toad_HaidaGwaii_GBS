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


#vcfR::write.vcf(VCF_original_test, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_3496SNPS_editlocinamesinR.vcf.gz")

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






#######################
## OBS HET FILTERING  ####
######################

#########################################
########## OBS HET PER REGION FILTER #####
####################################

#### genind dartR from genligght ####
genind.toad<-dartR::gl2gi(gl.toad, probar = FALSE, verbose = NULL)
genind.toad
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





#######################
## HWE filtering per pond p=0.01 ####
######################
#https://grunwaldlab.github.io/Population_Genetics_in_R/Locus_Stats.html


## check new pop.data file
dim(pop.data) #766

pop(gl.toad.obshetPerregion)

#### genind dartR from genligght ####
genid.toad.obshetPerregion<-dartR::gl2gi(gl.toad.obshetPerregion, probar = FALSE, verbose = NULL)

"/// GENIND OBJECT /////////

 // 766 individuals; 3,356 loci; 6,712 alleles; size: 21.6 Mb

 // Basic content
   @tab:  766 x 6712 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 6712 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 48-284)"


library("pegas")
#(nanhwe.full <- hw.test(genid.toad.obshetPerregion, B = 1000)) # performs 1000 permuatations

# check pop data is right file....
str(pop.data)
pop.data<-pop.data.766
str(pop.data)


## hwe per pop

pop(genid.toad.obshetPerregion)<-pop.data$pop

### check how many pops
pop(genid.toad.obshetPerregion) # 44

#run hwe test for each pop seprately - and only focus on the analytical P value (B=0)
(nanhwe.obshetPerregion.pop <- seppop(genid.toad.obshetPerregion) %>% lapply(hw.test, B = 0))

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

dim(nanhwe.obshetperregion.mat.df.filt.hwe ) # 1371   44

str(nanhwe.obshetperregion.mat.df.filt.hwe)

summary(nanhwe.obshetperregion.mat.df.filt.hwe)

hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe <- rownames(nanhwe.obshetperregion.mat.df.filt.hwe)
str(hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe) #  chr [1:1371] "ANBO_4657:127" "ANBO_13095:16" "ANBO_13454:49" "ANBO_18575:14" ...



gl.toad.obshetPerregion.hwe <- gl.toad.obshetPerregion[ , hwe.keep.loci_nanhwe.obshetperregion.mat.df.filt.hwe]
gl.toad.obshetPerregion.hwe

"/// GENLIGHT OBJECT /////////

 // 766 genotypes,  1,371 binary SNPs, size: 2 Mb
 132241 (12.59 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  1371 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-284)
   @other: a list containing: elements without names  "

#################################
###        FIS filter  ###################
#########################



#### genind dartR from genligght ####
genid.toad.obshetPerregion.hwe<-dartR::gl2gi(gl.toad.obshetPerregion.hwe, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.hwe

"/// GENIND OBJECT /////////

 // 766 individuals; 1,371 loci; 2,742 alleles; size: 8.8 Mb

 // Basic content
   @tab:  766 x 2742 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2742 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 48-284)"

pop(genid.toad.obshetPerregion.hwe)<-pop.data$pop
pop(genid.toad.obshetPerregion.hwe)

### check ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe<-genind2hierfstat(genid.toad.obshetPerregion.hwe)

#hierf.toad.obshetPerregion.hwe<-genind2hierfstat(genid.toad.obshetPerregion.hwe,pop=TRUE)
#hierf.toad.obshetPerregion.hwe$pop<-pop.data$pop

hierf.toad.obshetPerregion.hwebasic<-basic.stats(hierf.toad.obshetPerregion.hwe)

hierf.toad.obshetPerregion.hwebasic$Fis

head(hierf.toad.obshetPerregion.hwebasic$Fis)



### add new row names column
#hierf.toad.obshetPerregion.hwebasic$SNP<-xxx

FISdataperpond<-hierf.toad.obshetPerregion.hwebasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_766INDIV_1371SNPS_obshet06_hwe001_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()



######### filter

# snps I want to keep
nrow(FISdataperpondavg[which(FISdataperpondavg$Average<0.25 & FISdataperpondavg$Average>-0.31),]) 
# 1368

SNPkeepaveFIS<-FISdataperpondavg[which(FISdataperpondavg$Average<0.25 & FISdataperpondavg$Average>-0.31),]

SNPkeep<-SNPkeepaveFIS$SNP
str(SNPkeep) # chr [1:1368]

# edit names

SNPkeep2<- sub( '\\.', ':', SNPkeep)

str(SNPkeep2)

gl.toad.obshetPerregion.hwe.FIS <- gl.toad.obshetPerregion.hwe[ , SNPkeep2]
gl.toad.obshetPerregion.hwe.FIS

" /// GENLIGHT OBJECT /////////

 // 766 genotypes,  1,368 binary SNPs, size: 2 Mb
 131923 (12.59 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-284)
   @other: a list containing: elements without names  "


#######

#### genind dartR from genligght ####
genid.toad.obshetPerregion.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.hwe.FIS
""

pop(genid.toad.obshetPerregion.hwe.FIS)<-pop.data$pop
pop(genid.toad.obshetPerregion.hwe.FIS)

### check after filter ###


#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.hwe.FIS)

#hierf.toad.obshetPerregion.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.hwe.FIS,pop=TRUE)
#hierf.toad.obshetPerregion.hwe.FIS$pop<-pop.data$pop

hierf.toad.obshetPerregion.hwe.FISbasic<-basic.stats(hierf.toad.obshetPerregion.hwe.FIS)

hierf.toad.obshetPerregion.hwe.FISbasic$Fis

head(hierf.toad.obshetPerregion.hwe.FISbasic$Fis)



### add new row names column
#hierf.toad.obshetPerregion.hwe.FISbasic$SNP<-xxx

FISdataperpond<-hierf.toad.obshetPerregion.hwe.FISbasic$Fis

FISdataperpondavg <- data.frame("SNP"=rownames(FISdataperpond), "Average"=rowMeans(FISdataperpond, na.rm=T))
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/FISperpond_plot_766INDIV_1368SNPS_obshet06_hwe001_FIS_100bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=100))
dev.off()

#png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/FISperpond_plot_766INDIV_1371SNPS_obshet06_hwe001_50bins.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(FISdataperpondavg$Average,breaks=seq(-1,1,l=50))
dev.off()


###########################################
#### convert genind to structure ####
###########################################
genid.toad.obshetPerregion.hwe.FIS

"/// GENIND OBJECT /////////

 // 766 individuals; 1,368 loci; 2,736 alleles; size: 8.8 Mb

 // Basic content
   @tab:  766 x 2736 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 2-2)
   @loc.fac: locus factor for the 2736 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 2-2)
   @type:  codom
   @call: df2genind(X = xx[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
    pop = gl@pop, NA.char = "-", ploidy = 2)

 // Optional content
   @pop: population of each individual (group size range: 48-284)"

########### genind2structureL ######
#https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R

source("F:/GBS_data_03_02_21/genind2structureL_function.R")

## pop as four clusters
pop(genid.toad.obshetPerregion.hwe.FIS)
pop(genid.toad.obshetPerregion.hwe.FIS)<-pop.data$fourclusters

genind2structureL(genid.toad.obshetPerregion.hwe.FIS, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/structure_766INDIV_1368SNPS_genind2strR_toad_obshetPerregion_HWE_FIS_4clusters.str", pops=TRUE)


# pop as site
pop(genid.toad.obshetPerregion.hwe.FIS)
pop(genid.toad.obshetPerregion.hwe.FIS)<-pop.data$pop
pop(genid.toad.obshetPerregion.hwe.FIS)


genind2structureL(genid.toad.obshetPerregion.hwe.FIS, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/structure_766INDIV_1368SNPS_SITEID.str", pops=TRUE)


### add numbers to pop column of str file... b/c str file has two rows for each indiv (2 alleles) - hence wy cna't do manually

strfile<-read.table("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/structure_766INDIV_1368SNPS_SITEID.str",sep = "\t",header=TRUE)

View(strfile)

head(strfile)
str(strfile)

strfile$ind


strfile$pop<-substr(strfile$ind,1,9)
View(strfile)



#number ponds

strfile$pop <- as.numeric(as.factor(strfile$pop ))

View(strfile)
str(strfile)





write.table(strfile, file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/structure_766INDIV_1368SNPS_genind2strR_toad_obshetPerregion_HWE_FIS_SITEID_numbered.str",sep = "\t")


## save just numbers

strfile_num<-strfile[,1:2]
str(strfile_num)
View(strfile_num)

strfile_num_noheader<-strfile_num

#rm header
names(strfile_num_noheader)<- NULL

head(strfile_num_noheader)
#save shortened file
write.table(strfile_num_noheader,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/popmap766_BOTHLOCI_samples_header_siteID_rmINdiv_sorted_numbered.txt",quote=FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)
#rm(pop.data)

#######################################################
## SUBSET VCF TO LOCI AND INDIV FILTERED ########
#####################################################

## subset in vcftools and then use excel to subset regions, then upload each region

#VCF<-read.vcfR("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS_editlocinamesinR.vcf")


head(VCF@fix[,'ID'])

# REMOVE SNPS TO KEEP THAT HAVE BEEN THROUGH OBS HET, HWE AND FIS FILTER
VCF.FIS.SNPS.1368<-subset(VCF, VCF@fix[,'ID'] %in% SNPkeep2)
VCF.FIS.SNPS.1368

"***** Object of Class vcfR *****
766 samples
34 CHROMs
1,368 variants
Object size: 22 Mb
0 percent missing data
*****        *****         *****"

## save VCF

vcfR::write.vcf(VCF.FIS.SNPS.1368, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_1368SNPS_editlocinamesinR.vcf.gz")


### VCFTOOLS ####

## now upload to server along with rm sibs list, and remove sibs using vcftools, then download the new file, upload to here.... 766INDIV_1368SNPS_editlocinamesinRvcftools.recode.vcf

# vcftools --vcf 766INDIV_1368SNPS_editlocinamesinR.vcf --remove listindiv_plink_sub_sib_rel_over_04_names_unique.indv --recode --recode-INFO-all --out 766INDIV_1368SNPS_editlocinamesinRvcftools




### SUBSET VCF REMOVE INDIV
str(pop.data)

pop.data$sample.id


indiv_samples_keep<-pop.data$sample.id

str(indiv_samples_keep) #  chr [1:766]

#should be 766, not 509
VCF.FIS.SNPS.1368[,indiv_samples_keep]

""

test<-VCF.FIS.SNPS.1368[,indiv_samples_keep]

""


### IT'S CUTTING OUT FORMAT COLUMN SO ADD TO LIST

# add list at top called "FORMAT"
indiv_samples_keep

indiv_samples_keep_FORMAT<- append('FORMAT',indiv_samples_keep)

#WORKED!!!
VCF.FIS.SNPS.1368_766<-VCF.FIS.SNPS.1368[,indiv_samples_keep_FORMAT]

""




#######################
###### SFS SITE FREQUENCY SPECTRUM ##########
#########################


###########################################
#### SPLIT INTO four REGIONS ####
###########################################


pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.hwe.FIS)
gl.toad.obshetPerregion.hwe.FIS.bypop

"$HaidaGwai
 /// GENLIGHT OBJECT /////////

 // 239 genotypes,  1,368 binary SNPs, size: 730.1 Kb
 43394 (13.27 %) missing data

 // Basic content
   @gen: list of 239 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  239 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 239-239)
   @other: a list containing: elements without names 


$LowerMain
 /// GENLIGHT OBJECT /////////

 // 284 genotypes,  1,368 binary SNPs, size: 825.1 Kb
 46630 (12 %) missing data

 // Basic content
   @gen: list of 284 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  284 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 284-284)
   @other: a list containing: elements without names 


$Northwest
 /// GENLIGHT OBJECT /////////

 // 48 genotypes,  1,368 binary SNPs, size: 231.9 Kb
 6206 (9.45 %) missing data

 // Basic content
   @gen: list of 48 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  48 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-48)
   @other: a list containing: elements without names 


$VanIsland
 /// GENLIGHT OBJECT /////////

 // 195 genotypes,  1,368 binary SNPs, size: 618.4 Kb
 35693 (13.38 %) missing data

 // Basic content
   @gen: list of 195 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  195 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-195)
   @other: a list containing: elements without names "



gl.toad.obshetPerregion.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.hwe.FIS.VanIsland

"  /// GENLIGHT OBJECT /////////

 // 195 genotypes,  1,368 binary SNPs, size: 618.4 Kb
 35693 (13.38 %) missing data

 // Basic content
   @gen: list of 195 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  195 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 195-195)
   @other: a list containing: elements without names   "



gl.toad.obshetPerregion.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.hwe.FIS.LowerMain
" /// GENLIGHT OBJECT /////////

 // 284 genotypes,  1,368 binary SNPs, size: 825.1 Kb
 46630 (12 %) missing data

 // Basic content
   @gen: list of 284 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  284 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 284-284)
   @other: a list containing: elements without names   
"

gl.toad.obshetPerregion.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.hwe.FIS.HaidaGwai
"  /// GENLIGHT OBJECT /////////

 // 239 genotypes,  1,368 binary SNPs, size: 730.1 Kb
 43394 (13.27 %) missing data

 // Basic content
   @gen: list of 239 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  239 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 239-239)
   @other: a list containing: elements without names  "



gl.toad.obshetPerregion.hwe.FIS.Northwest<-gl.toad.obshetPerregion.hwe.FIS.bypop$Northwest
gl.toad.obshetPerregion.hwe.FIS.Northwest

" /// GENLIGHT OBJECT /////////

 // 48 genotypes,  1,368 binary SNPs, size: 231.9 Kb
 6206 (9.45 %) missing data

 // Basic content
   @gen: list of 48 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  48 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-48)
   @other: a list containing: elements without names  "

dartR::gl.report.maf(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, maf.limit = 0.5, ind.limit = 0, loc.limit = 0,
                     v = 2)

### convert genlight to genind objects 
genind.toad.obshetPerregion.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.obshetPerregion.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.obshetPerregion.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)

###########################################
#### plot folded MAF - as count of alleles, for each region ####
###########################################

# split VCF based on region



### HG #####
##############################
VCF.FIS.SNPS.1368

"***** Object of Class vcfR *****
766 samples
34 CHROMs
1,368 variants
Object size: 22 Mb
0 percent missing data
*****        *****         *****"

### SUBSET HAIDA GWAII

gl.toad.obshetPerregion.hwe.FIS.HaidaGwai@ind.names

HG_data_239_1368_indnames<-gl.toad.obshetPerregion.hwe.FIS.HaidaGwai@ind.names

HG_data_239_1368_indnames_FORMAT<-append('FORMAT',HG_data_239_1368_indnames)

VCF.FIS.SNPS.1368.HG_239_1368<-VCF.FIS.SNPS.1368[,HG_data_239_1368_indnames_FORMAT]
VCF.FIS.SNPS.1368.HG_239_1368

"***** Object of Class vcfR *****
239 samples
34 CHROMs
1,368 variants
Object size: 6.9 Mb
0 percent missing data
*****        *****         *****"


#minor allele freq
VCF.HG_239_1368_mafvcf<-maf(VCF.FIS.SNPS.1368.HG_239_1368,element=2)

VCF.HG_239_1368_mafvcf.df<-as.data.frame(VCF.HG_239_1368_mafvcf)
summary(VCF.HG_239_1368_mafvcf.df)

head(VCF.HG_239_1368_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_239INDIV_766INDIV_1368SNPS_freq_cutindvR_maf_xlim0_1.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_239_1368_mafvcf.df$Frequency,breaks=seq(0,1,l=50), main="MAF Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_239INDIV_766INDIV_1368SNPS_freq_cutindvR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_239_1368_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main="MAF Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_HG_239INDIV_766INDIV_1368SNPS_COUNT_cutindvR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VCF.HG_239_1368_mafvcf.df$Count, main="MAC Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")
dev.off()



HGcounthist<-hist(VCF.HG_239_1368_mafvcf.df$Count, main=" ",
                  ylim = c(0, 2000),
                  xlab="Minor allele count",
                  ylab="Frequency")






HG_mafvcf.df<-VCF.HG_239_1368_mafvcf.df
hist(HG_mafvcf.df$Count, main="MAC Haida Gwaii filtered calc from maf() - cut indv in R ",
     ylim = c(0, 1500),
     xlab="Minor allele count",
     ylab="Frequency")

### VI #####
##############################

## SUBSET VanIsland
gl.toad.obshetPerregion.hwe.FIS.VanIsland

""

gl.toad.obshetPerregion.hwe.FIS.VanIsland@ind.names

VI_data_195_1368_indnames<-gl.toad.obshetPerregion.hwe.FIS.VanIsland@ind.names


VI_data_195_1368_indnames_FORMAT<-append('FORMAT',VI_data_195_1368_indnames)

VCF.FIS.SNPS.1368.VI_195_1368<-VCF.FIS.SNPS.1368[,VI_data_195_1368_indnames_FORMAT]
VCF.FIS.SNPS.1368.VI_195_1368


VI_mafvcf<-maf(VCF.FIS.SNPS.1368.VI_195_1368,element=2)


VI_mafvcf.df<-as.data.frame(VI_mafvcf)

head(VI_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_195INDIV_766NDIV_1368SNPS_freq_cutinR_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(VI_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" MAF Vancouver Island filtered cut in R",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_VI_195INDIV_766NDIV_1368SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
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
gl.toad.obshetPerregion.hwe.FIS.LowerMain

" "

gl.toad.obshetPerregion.hwe.FIS.LowerMain@ind.names

LM_data_284_1368_indnames<-gl.toad.obshetPerregion.hwe.FIS.LowerMain@ind.names


LM_data_284_1368_indnames_FORMAT<-append('FORMAT',LM_data_284_1368_indnames)

VCF.FIS.SNPS.1368.LM_284_1368<-VCF.FIS.SNPS.1368[,LM_data_284_1368_indnames_FORMAT]
VCF.FIS.SNPS.1368.LM_284_1368

""


LM_mafvcf<-maf(VCF.FIS.SNPS.1368.LM_284_1368,element=2)


LM_mafvcf.df<-as.data.frame(LM_mafvcf)

head(LM_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_284INDIV_766NDIV_1368SNPS_freq_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(LM_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_LM_284INDIV_766NDIV_1368SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
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
gl.toad.obshetPerregion.hwe.FIS.Northwest

""

gl.toad.obshetPerregion.hwe.FIS.Northwest@ind.names

NW_data_48_1368_indnames<-gl.toad.obshetPerregion.hwe.FIS.Northwest@ind.names


NW_data_48_1368_indnames_FORMAT<-append('FORMAT',NW_data_48_1368_indnames)

VCF.FIS.SNPS.1368.NW_48_1368<-VCF.FIS.SNPS.1368[,NW_data_48_1368_indnames_FORMAT]
VCF.FIS.SNPS.1368.NW_48_1368

""


NW_mafvcf<-maf(VCF.FIS.SNPS.1368.NW_48_1368,element=2)


NW_mafvcf.df<-as.data.frame(NW_mafvcf)

head(NW_mafvcf.df)

#mafvcf.df$allelname <- row.names(mafvcf.df) 

par(mar=c(4,4,4,4))


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_NW_34INDIV_766NDIV_1368SNPS_freq_maf.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(NW_mafvcf.df$Frequency,breaks=seq(0,0.5,l=50), main=" ",
     ylim = c(0, 1500),
     xlab="Minor allele frequency",
     ylab="Frequency")
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/MAFplot_NW_34INDIV_766NDIV_1368SNPS_count.png", width = 8, height = 6.5, units = 'in', res = 300)
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
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_766_1368.png", width = 3.5, height = 6.5, units = 'in', res = 300)
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
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_COUNT_vs_locicount_766INDIV_1368SNPS.png", width = 3.5, height = 6.5, units = 'in', res = 300)
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


ggsave("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/SFS_4clusters_MAF_vs_locicount_766INDIV_1368SNPS.svg", width = 3.5, height = 6.5, units = 'in')
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
pop(gl.toad.obshetPerregion.hwe.FIS)
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.obshetPerregion.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)

#gl.toad.obshetPerregion.hwe.FIS.subset
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
#  28.20326  3.92967

##per region
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_766INDIV_1368SNPS_obshet06_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
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
p<-p + ylab("PC2 (3.9% explained variance)") + xlab("PC1 (28.2% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p<-p + guides(colour=guide_legend(nrow=2))
p
dev.off()



##per pop
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$pop
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_perpop_4axes_766INDIV_1368SNPS_obshet06_hwe001_FIS.png", width = 8, height = 6.5, units = 'in', res = 300)
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



## NEW str
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_str_4clusters_4axes_766INDIV_1368SNPS_obshet06_hwe001_FIS_orangegrad_legright2.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (3.9% explained variance)") + xlab("PC1 (28.2% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()

# new col
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_766INDIV_1368SNPS_obshet06_hwe001_FIS_orangegrad_legright2.pdf",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#ab4e03", "#fca45d","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC2 (3.9% explained variance)") + xlab("PC1 (28.2% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()





##############################
####### PC3 and 4 #############
##per region
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters


my_pca <- glPca(gl.toad.obshetPerregion.hwe.FIS, nf=6, parallel = require("parallel"))


my_scores <- as.data.frame(my_pca$scores)
my_scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)

#storing eigenvalues on a file


#Checking how many eigenvalues to keep
#plotting PCA
toad.pca.scores <- as.data.frame(my_pca$scores)
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)

#gl.toad.obshetPerregion.hwe.FIS.subset
"
"

library(ggplot2)

# get % pca for 1st two axes ????? it created 402 PCA axes, why????



my_pca$eig[c(3,4)]

#not a percent!!!!!!!!!
sum(my_pca$eig)

#### percentage variance explained for PC axes
100*my_pca$eig/sum(my_pca$eig)


100*my_pca$eig[c(3,4)]/sum(my_pca$eig)
# 1.5398337 0.9851565

##per region
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters
toad.pca.scores$pop <- pop(gl.toad.obshetPerregion.hwe.FIS)


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_col_4clusters_4axes_766INDIV_3496SNPS_hweperpond_FISperpond_PC3_PC4.png", width = 8, height = 6.5, units = 'in', res = 300)
set.seed(9)
p <- ggplot(toad.pca.scores, aes(x=PC3, y=PC4, colour=pop))
p <- p + geom_point(size=4, alpha=.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC4 (1.0% explained variance)") + xlab("PC3 (1.5% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="top",legend.text=element_text(size=16))
p
dev.off()

## NEW COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/PCA/PCA_str_4clusters_766INDIV_3496SNPS_hweperpond_FISperpond_orangegrad_legright2_PC3_PC4.png",  bg = "transparent",width =300, height = 150, units = c("mm"))
p <- ggplot(toad.pca.scores, aes(x=PC3, y=PC4, colour=pop))
p <- p + geom_point(size=4, alpha=0.9)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
#p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + theme_bw()
p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
                            name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))
p<-p + labs(colour='Region') 
p<-p + ylab("PC4 (1.0% explained variance)") + xlab("PC3 (1.5% explained variance)") + theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 20), title = element_text(size = 20))
p<-p + theme(legend.position="right",legend.text=element_text(size=18),legend.key.height=unit(2,"line"),panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank())
p
dev.off()




######################################
########## FST FIS ##############
########################################

genid.toad.obshetPerregion.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS, probar = FALSE, verbose = NULL)
genid.toad.obshetPerregion.hwe.FIS

" "

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.hwe.FIS,pop=NULL)





#has pop - as fourclusters
hierf.toad.obshetPerregion.hwe.FIS$pop

hierf.toad.obshetPerregion.hwe.FIS$pop <- pop.data$fourclusters
hierf.toad.obshetPerregion.hwe.FIS$pop <- factor(pop.data$fourclusters)

FST_hierf.toad.obshetPerregion.hwe.FIS<-pairwise.WCfst(hierf.toad.obshetPerregion.hwe.FIS,diploid=T)
FST_hierf.toad.obshetPerregion.hwe.FIS

"             HaidaGwai  LowerMain  Northwest  VanIsland
HaidaGwai        NA 0.43919665 0.64063228 0.48002269
LowerMain 0.4391967         NA 0.09246912 0.06552682
Northwest 0.6406323 0.09246912         NA 0.11186282
VanIsland 0.4800227 0.06552682 0.11186282         NA    "

write.table(FST_hierf.toad.obshetPerregion.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.txt",sep = "\t",row.names = TRUE,col.names = TRUE )



write.csv(FST_hierf.toad.obshetPerregion.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_4clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.csv",sep = "\t",row.names = TRUE,col.names = TRUE )


## fst twoclusters

hierf.toad.obshetPerregion.hwe.FIS$pop

hierf.toad.obshetPerregion.hwe.FIS$pop <- pop.data$twoclusters
hierf.toad.obshetPerregion.hwe.FIS$pop <- factor(pop.data$twoclusters)

FST_hierf.toad.obshetPerregion.hwe.FIS<-pairwise.WCfst(hierf.toad.obshetPerregion.hwe.FIS,diploid=T)
FST_hierf.toad.obshetPerregion.hwe.FIS

"            CoastalBC HaidaGwai
CoastalBC        NA 0.3625071
HaidaGwai 0.3625071        NA       "

write.table(FST_hierf.toad.obshetPerregion.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_2clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.txt",sep = "\t",row.names = TRUE,col.names = TRUE )



write.csv(FST_hierf.toad.obshetPerregion.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_2clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.csv",sep = "\t",row.names = TRUE,col.names = TRUE )




# as pop
hierf.toad.obshetPerregion.hwe.FIS$pop

hierf.toad.obshetPerregion.hwe.FIS$pop <- pop.data$pop
hierf.toad.obshetPerregion.hwe.FIS$pop <- factor(pop.data$pop)

trial<-pairwise.WCfst(hierf.toad.obshetPerregion.hwe.FIS,diploid=T)
trial

"               "

write.csv(trial, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_POP_766INDIV_3496SNPS_hwe26pops0.01_FISperpond.csv")








##### p value two areas
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$twoclusters

(Fstp<-stamppFst(gl.toad.obshetPerregion.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_2clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.csv" )

##### p value fourclusters
pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters

(Fstp<-stamppFst(gl.toad.obshetPerregion.hwe.FIS, nboots = 100, percent = 95, nclusters = 1))


write.csv(Fstp$Pvalues, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Fst_pvalue_4clusters_766NDIV_1368SNPS_obshet06_hwe001_FIS.csv" )





###################################
########### FST WITHIN two clusters #############
########################################

pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$twoclusters

gl.toad.obshetPerregion.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.hwe.FIS)


gl.toad.obshetPerregion.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.hwe.FIS.CoastalBC<-gl.toad.obshetPerregion.hwe.FIS.bypop$CoastalBC



### convert genlight to genind objects 
genind.toad.obshetPerregion.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.hwe.FIS.CoastalBC<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.CoastalBC, probar = FALSE, verbose = NULL)

""





### make new pop files for each dataset
pop.data_HaidaGwai <- filter(pop.data, twoclusters == "HaidaGwai")
pop.data_CoastalBC <- filter(pop.data, twoclusters == "CoastalBC")




##### HaidaGwai #############

# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.HaidaGwai)<-pop.data_HaidaGwai$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.HaidaGwai,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai$pop


(HaidaGwaiFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai,diploid=TRUE))




##### CoastalBC #############


# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.CoastalBC)<-pop.data_CoastalBC$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.CoastalBC<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.CoastalBC,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.CoastalBC$pop


(CoastalBCFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.CoastalBC,diploid=TRUE))


######## save values

HaidaGwaiFST1<-HaidaGwaiFST$FST
CoastalBCFST1<-CoastalBCFST$FST


FSTperpop<-cbind(HaidaGwaiFST1,CoastalBCFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_twoclusters_766INDIV_3496SNPS.txt")




###################################
########### FST WITHIN four clusters #############
########################################

pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.hwe.FIS)


gl.toad.obshetPerregion.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.hwe.FIS.VanIsland

"   "



gl.toad.obshetPerregion.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.hwe.FIS.LowerMain
"  
"

gl.toad.obshetPerregion.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.hwe.FIS.HaidaGwai
"   "



gl.toad.obshetPerregion.hwe.FIS.Northwest<-gl.toad.obshetPerregion.hwe.FIS.bypop$Northwest
gl.toad.obshetPerregion.hwe.FIS.Northwest

### convert genlight to genind objects 
genind.toad.obshetPerregion.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.obshetPerregion.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.obshetPerregion.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)


""



""





### make new pop files for each dataset
pop.data_HaidaGwai <- filter(pop.data, fourclusters == "HaidaGwai")
pop.data_LowerMain <- filter(pop.data, fourclusters == "LowerMain")
pop.data_VanIsland <- filter(pop.data, fourclusters == "VanIsland")
pop.data_Northwest <- filter(pop.data, fourclusters == "Northwest")



##### HaidaGwai #############

# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.HaidaGwai)<-pop.data_HaidaGwai$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.HaidaGwai,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai$pop


(HaidaGwaiFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.HaidaGwai,diploid=TRUE))



##### LowerMain #############

# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.LowerMain)<-pop.data_LowerMain$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.LowerMain<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.LowerMain,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.LowerMain$pop


(LowerMainFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.LowerMain,diploid=TRUE))




##### VanIsland #############

# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.VanIsland)<-pop.data_VanIsland$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.VanIsland<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.VanIsland,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.VanIsland$pop


(VanIslandFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.VanIsland,diploid=TRUE))



##### Northwest #############

# make pop breeding pond
pop(genind.toad.obshetPerregion.hwe.FIS.Northwest)<-pop.data_Northwest$pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS.Northwest<-genind2hierfstat(genind.toad.obshetPerregion.hwe.FIS.Northwest,pop=NULL)

#has pop - as breeding pond
hierf.toad.obshetPerregion.hwe.FIS.Northwest$pop


(NorthwestFST<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS.Northwest,diploid=TRUE))







######## save values

HaidaGwaiFST1<-HaidaGwaiFST$FST
LowerMainFST1<-LowerMainFST$FST
VanIslandFST1<-VanIslandFST$FST
NorthwestFST1<-NorthwestFST$FST

FSTperpop<-cbind(HaidaGwaiFST1,LowerMainFST1,VanIslandFST1,NorthwestFST1)


write.table(FSTperpop[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/FST_within_fourclusters_766INDIV_3496SNPS.txt")








##########################################
### FST per loc overall pops ########

hierf.toad.obshetPerregion.hwe.FIS$pop <- pop.data$twoclusters
hierf.toad.obshetPerregion.hwe.FIS$pop <- factor(pop.data$twoclusters)

basic.nosibs.hwe.FIS<-basic.stats(hierf.toad.obshetPerregion.hwe.FIS)

perloc.basic.nosibs.hwe.FIS<-basic.nosibs.hwe.FIS$perloc



######## wc function

(FSTperloc<-hierfstat::wc(hierf.toad.obshetPerregion.hwe.FIS,diploid=TRUE))
perlocobject<-FSTperloc$per.loc

perlocobject$FST

summary(perlocobject$FST)

hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50))

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Fst/Hist_Fst_perlocus_2clusters_766INDIV_3496SNPS_hwe26pops0.01_FISperpond_wcfunc.png", width = 8, height = 6.5, units = 'in', res = 300)
hist(perlocobject$FST,breaks=seq(-0.006251 ,1,l=50),
     main = " ",
     xlab="Pairwise Fst per locus between Haida Gwaii and coastal BC",
     ylab="Frequency")
dev.off()





##################################
###########  Ext Het & Obs het FIS ####################
############################

gl.toad.obshetPerregion.hwe.FIS

"
 // 766 genotypes,  1,368 binary SNPs, size: 2 Mb
 131923 (12.59 %) missing data

 // Basic content
   @gen: list of 766 SNPbin
   @ploidy: ploidy of each individual  (range: 2-2)

 // Optional content
   @ind.names:  766 individual labels
   @loc.names:  1368 locus labels
   @chromosome: factor storing chromosomes of the SNPs
   @position: integer storing positions of the SNPs
   @pop: population of each individual (group size range: 48-284)
   @other: a list containing: elements without names "


#### split matrix into 3 regions ###########

pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data$fourclusters

gl.toad.obshetPerregion.hwe.FIS.bypop<-seppop(gl.toad.obshetPerregion.hwe.FIS)


gl.toad.obshetPerregion.hwe.FIS.VanIsland<-gl.toad.obshetPerregion.hwe.FIS.bypop$VanIsland
gl.toad.obshetPerregion.hwe.FIS.VanIsland

"   "



gl.toad.obshetPerregion.hwe.FIS.LowerMain<-gl.toad.obshetPerregion.hwe.FIS.bypop$LowerMain
gl.toad.obshetPerregion.hwe.FIS.LowerMain
"  
"

gl.toad.obshetPerregion.hwe.FIS.HaidaGwai<-gl.toad.obshetPerregion.hwe.FIS.bypop$HaidaGwai
gl.toad.obshetPerregion.hwe.FIS.HaidaGwai
"   "



gl.toad.obshetPerregion.hwe.FIS.Northwest<-gl.toad.obshetPerregion.hwe.FIS.bypop$Northwest
gl.toad.obshetPerregion.hwe.FIS.Northwest

" "

#dartR::gl.report.maf(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, maf.limit = 0.5, ind.limit = 0, loc.limit = 0,
                     v = 2)

### convert genlight to genind objects 
genind.toad.obshetPerregion.hwe.FIS.VanIsland<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.VanIsland, probar = FALSE, verbose = NULL)


""


genind.toad.obshetPerregion.hwe.FIS.LowerMain<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.LowerMain, probar = FALSE, verbose = NULL)

""


genind.toad.obshetPerregion.hwe.FIS.HaidaGwai<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.HaidaGwai, probar = FALSE, verbose = NULL)

''

genind.toad.obshetPerregion.hwe.FIS.Northwest<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS.Northwest, probar = FALSE, verbose = NULL)

############### summary on each one
genind.toad.obshetPerregion.hwe.FIS.VanIsland_sum<-adegenet::summary(genind.toad.obshetPerregion.hwe.FIS.VanIsland)

genind.toad.obshetPerregion.hwe.FIS.LowerMain_sum<-adegenet::summary(genind.toad.obshetPerregion.hwe.FIS.LowerMain)

genind.toad.obshetPerregion.hwe.FIS.HaidaGwai_sum<-adegenet::summary(genind.toad.obshetPerregion.hwe.FIS.HaidaGwai)


genind.toad.obshetPerregion.hwe.FIS.Northwest_sum<-adegenet::summary(genind.toad.obshetPerregion.hwe.FIS.Northwest)


## rename
divVanIsland<-genind.toad.obshetPerregion.hwe.FIS.VanIsland_sum
divLowerMain<-genind.toad.obshetPerregion.hwe.FIS.LowerMain_sum
divHaidaGwai<-genind.toad.obshetPerregion.hwe.FIS.HaidaGwai_sum
divNorthwest<-genind.toad.obshetPerregion.hwe.FIS.Northwest_sum


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

write.table(Hetsum[1,], "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exp_Het_summary_4regions_766INDIV_1368SNPS.txt")





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
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_766INDIV_1368SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
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
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Boxplot_expect_Het_per_region_766INDIV_1368SNPS_hwe26pops_FISperpond0.1_noylim.png", width = 8, height = 6.5, units = 'in', res = 600)
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

write.csv(exphet_SEdf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_mean_sd_SE_CI_MOE_766INDIV_1368SNPS.csv")


ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_95pctCI_4clusters_766INDIV_1368SNPS_correctCI.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
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
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Exphet_95pctCI_4clusters_766INDIV_1368SNPS_correctCI.pdf",  bg = "transparent",width =200, height = 180, units = c("mm"))
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




########################## BASIC STATS FIS ###################

#### genind dartR from genligght ####
genid.toad.obshetPerregion.hwe.FIS<-dartR::gl2gi(gl.toad.obshetPerregion.hwe.FIS, probar = FALSE, verbose = NULL)

#set pop as region
pop(genid.toad.obshetPerregion.hwe.FIS)<-pop.data$fourclusters


#has pop
genid.toad.obshetPerregion.hwe.FIS@pop

#### heirfstat conversion from genind ########
hierf.toad.obshetPerregion.hwe.FIS<-genind2hierfstat(genid.toad.obshetPerregion.hwe.FIS,pop=NULL)

#has pop - as region
hierf.toad.obshetPerregion.hwe.FIS$pop



dataIN<-hierf.toad.obshetPerregion.hwe.FIS
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

forb2.pop.stats.hwe.FIS <- pop.stats(hierf.toad.obshetPerregion.hwe.FIS)
forb2.pop.stats.hwe.FIS

### filter max het!! no negative FIS

"        Hs_mean  Hs_SD Ho_mean  Ho_SD    P     Fis  Fis_ll  Fis_hl
HaidaGwai  0.0079 0.0378  0.0084 0.0423 23.2 -0.0616 -0.1048 -0.0158
LowerMain  0.1399 0.0917  0.1383 0.0936 98.2  0.0114  0.0041  0.0185
Northwest  0.1007 0.1402  0.1061 0.1508 53.4 -0.0538 -0.0695 -0.0379
VanIsland  0.1381 0.1074  0.1361 0.1068 93.2  0.0142  0.0058  0.0225
OVERALL    0.0966 0.1146  0.0972 0.1182   NA -0.0062      NA      NA "

### pop stats output #####

write.table(forb2.pop.stats.hwe.FIS, "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/hierfstat_popstats_region_766INDIV_1368SNPS.txt",sep = "\t",row.names = TRUE,col.names = TRUE)


########################################
############ Fis per region CI ################
#########################################
str(forb2.pop.stats.hwe.FIS)

forb2.pop.stats.hwe.FIS$Fis

Fis_95CI_region<-cbind(forb2.pop.stats.hwe.FIS$Fis,forb2.pop.stats.hwe.FIS$Fis_ll,forb2.pop.stats.hwe.FIS$Fis_hl)
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

#export
write.csv(dfFis_95CI_region,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_766INDIV_1368SNPS.csv")

## plot - no legend
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Inbreeding (Fis ± 95% CIs)")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p


## formatted correctly
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_766INDIV_1368SNPS_symbolsforamtted.png", width = 8, height = 6.5, units = 'in', res = 600)
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
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_766INDIV_1368SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
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


## new str

ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_str_95CI_4clusters_766INDIV_1368SNPS.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(dfFis_95CI_region, aes(x = Region, color = Region, y=Fis)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii","Northwest BC"))+
  geom_errorbar(aes(ymax = Fis_hl, ymin = Fis_ll),
                position = "dodge",
                width=0.2,
                size=0.5)+
  labs(x="Region",y=expression(bold(paste("Inbreeding coefficient  (",italic(F)[plain(IS)],paste(")")))))+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
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
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Diversity/Fis_95CI_per_region_766INDIV_1368SNPS_neworange.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
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

grp <- find.clusters(gl.toad.obshetPerregion.hwe.FIS, max.n.clust=26)


# save variance explained PC plot
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_varianceexplainedbyPCA_findclusters_maxclusters26_766INDIV_1368SNPS.png")



# https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# pick largest
#Choose the number PCs to retain (>=1): 
600




# save BIC plot

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_BIC_findclusters_maxclusters26_600PCs_766INDIV_1368SNPS.png")



#Choose the number of clusters (>=2): 
2

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_600PCs_2clusters_766INDIV_1368SNPS.csv")



############################
####### 200 PCs, 2 groups #######

grp <- find.clusters(gl.toad.obshetPerregion.hwe.FIS, max.n.clust=26)


#Choose the number PCs to retain (>=1): 
200


#Choose the number of clusters (>=2): 
2

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters26_200PCs_2clusters_766INDIV_1368SNPS.csv")



############################
####### 200 PCs, 4 groups #######

grp <- find.clusters(gl.toad.obshetPerregion.hwe.FIS, max.n.clust=44)


# save variance explained PC plot
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_varianceexplainedbyPCA_findclusters_maxclusters44_766INDIV_1368SNPS.png")



#Choose the number PCs to retain (>=1): 

600

# save BIC plot

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_BIC_findclusters_maxclusters44_600PCs_766INDIV_1368SNPS.png")



#Choose the number of clusters (>=2): 
4

grp$grp

# saving cluster allocation
write.csv(grp$grp, file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters44_600PCs_4clusters_766INDIV_1368SNPS.csv")

DAPC_varianceexplainedbyPCA_findclusters_maxclusters44_766INDIV_1368SNPS

DAPC_BIC_findclusters_maxclusters44_600PCs_4clusters_766INDIV_1368SNPS

##############################
######## DAPC CLUSTER 4 ###########
###########################



clustass4clusters<-read.csv( file = "F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/clustRumAssign_findclusters_maxclusters44_600PCs_4clusters_766INDIV_1368SNPS.csv")



clustass4clusters<-as.data.frame(clustass4clusters)

str(clustass4clusters)

clustass4clusters$sample.id<-clustass4clusters$X

head(clustass4clusters)

clustass4clusters$fourclusters600pca<-clustass4clusters$x

head(clustass4clusters)


#JUST SAVE 2 COLS
clustass4clusterssub<-clustass4clusters[c(3:4)]

head(clustass4clusterssub)




### join w pop data 
pop.data4clusters<-left_join(pop.data,clustass4clusterssub)

head(pop.data4clusters)


View(pop.data4clusters)

# work out which pop corresponds to which cluster number
#View(pop.data4clusters)


pop.data4clusters$fourclusters600pcaNAMES<-pop.data4clusters$fourclusters600pca

pop.data4clusters$fourclusters600pcaNAMES[which(pop.data4clusters$fourclusters600pca=="1")] <- "LowerMain"

pop.data4clusters$fourclusters600pcaNAMES[which(pop.data4clusters$fourclusters600pca=="4")] <- "VanIsland"

pop.data4clusters$fourclusters600pcaNAMES[which(pop.data4clusters$fourclusters600pca=="2")] <- "HaidaGwai"

pop.data4clusters$fourclusters600pcaNAMES[which(pop.data4clusters$fourclusters600pca=="3")] <- "Northwest"

head(pop.data4clusters)



pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data4clusters$fourclusters600pcaNAMES


pop(gl.toad.obshetPerregion.hwe.FIS) <- pop.data4clusters$fourclusters


## dapc with max pcs and das
pnw.dapc <- dapc(gl.toad.obshetPerregion.hwe.FIS, n.pca = 600, n.da = 10, parallel = require("parallel"))




#### check optimal number of PCs = 48
temp <- optim.a.score(pnw.dapc)
temp
#ascoreopt_PC_28_maxclusters44_600PCs_4clusters_766INDIV_1368SNPS


# run dapc with optimal no. pcs
pnw.dapc2 <- dapc(gl.toad.obshetPerregion.hwe.FIS, n.pca = 48, n.da = 10, parallel = require("parallel"))


scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_766INDIV_1368SNPS_optimal_48PCaxes_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()



# reduce da axes from 10 to 4 to test - makes no difference
pnw.dapc3 <- dapc(gl.toad.obshetPerregion.hwe.FIS, n.pca = 18, n.da = 4, parallel = require("parallel"))

scatter(pnw.dapc3, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)



#p <- p + scale_color_manual(values=c("#E69F00", "#56B4E9","#009E73","pink"),
#                       name="Region", breaks = c("VanIsland","LowerMain", "HaidaGwai","Northwest"), labels=c#("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))


labs <- c("HaidaGwai", "LowerMain","NorthWest","VanIsland")
cols=c("#009E73","#56B4E9", "pink","#E69F00" )


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_766INDIV_1368SNPS_optimal_48PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc2, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_766INDIV_1368SNPS_600pcaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()



### DAPC with 80% variation explained####

highdapc <- dapc(gl.toad.obshetPerregion.hwe.FIS, parallel = require("parallel"), n.da = 10, pca.select = "percVar", perc.pca=80)



png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_766INDIV_1368SNPS_80pctPCs_defaultcols.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75)
dev.off()


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_col_4_clusters_766INDIV_1368SNPS_80pctPCs.png", width = 12, height = 6.5, units = 'in', res = 300)
scatter(highdapc, col = cols, cex = 2, legend = TRUE, clabel = F, 
        posi.leg = "bottomleft", scree.pca = T,posi.da = "bottomright",
        posi.pca = "topright", cleg = 0.75,txt.leg =labs)
dev.off()





#Loci contributing to observed differences, threshold set arbitrarily?
set.seed(4)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)

png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/DAPC_loadings_col_4_clusters_766INDIV_1368SNPS_optimal_23PCaxes.png", width = 12, height = 6.5, units = 'in', res = 300)
contrib <- loadingplot(pnw.dapc2$var.contr, axis = 2, thres = 0.02, lab.jitter = 1)
dev.off()




#### MEMBERSHIP PROP ######
#probability of population belonging to the pop it's assigned

set.seed(999)
pramx <- xvalDapc(tab(gl.toad.obshetPerregion.hwe.FIS), pop(gl.toad.obshetPerregion.hwe.FIS), parallel = "snow")
###--->40



compoplot(pnw.dapc2,col = brewer.pal(4, "Paired"), posi = 'top')


dapc.results <- as.data.frame(pnw.dapc2$posterior)
dapc.results$oldpop <- pop(gl.toad.obshetPerregion.hwe.FIS)
dapc.results$indNames <- rownames(dapc.results)


dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample",
                            "Assigned_Pop","Posterior_membership_probability")


labs <- c("HaidaGwai", "LowerMain","NorthWest","VanIsland")
cols=c("#009E73","#56B4E9", "pink","#E69F00" )


###### no sep areas 
png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/DAPC/Membershipprob_col_4_clusters_766INDIV_1368SNPS_optimalPC_18PCaxes_nosepareas.png", width = 12, height = 3, units = 'in', res = 300)
t <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
t <- t + geom_bar(stat='identity') 
t <- t + theme(axis.text.x = element_blank(),axis.title.y = element_text( size = 10),axis.ticks.x = element_blank(),panel.background = element_blank())
t <- t + scale_fill_manual(values=c( "#009E73","#56B4E9", "pink","#E69F00"),
                           name="Assigned Pop K=4 cluster", breaks = c("HaidaGwai", "LowerMain","NorthWest","VanIsland"), labels=c("HaidaGwai", "LowerMain","NorthWest","VanIsland"))
t
dev.off()





