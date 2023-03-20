

library(dplyr)
### edit populations.sumstats_summary file - all sites - copy and paste all sites part of file into excel sheet and save as CSV and move into Pi directory #####



############ 27 ponds ##################


# open file

popsumstats<-read.csv("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/populations.sumstats_summary_allsites_regioneditman.csv", sep = ",")

str(popsumstats)


# amke subregion
#popsumstats$subregion<-popsumstats$ï..Pop_ID

#popsumstats$subregion[substr(popsumstats$subregion,1,3),]

#subset = ibd[substr(ibd$name1,1,9)==substr(ibd$name2,1,9),]

# make region
#popsumstats$Region<-popsumstats$ï..Pop_ID
#popsumstats$Region[(1:14),] <- "VanIsland"

#pop.data_rm$twoclusters[pop.data_rm$twoclusters=="swBC"] <- "CoastalBC"
#pop.data_rm$twoclusters[pop.data_rm$twoclusters=="Northwest"] <- "CoastalBC"

View(popsumstats)




head(popsumstats)


str(popsumstats$Region)

popsumstats$Region<-as.factor(popsumstats$Region)



## reorder
popsumstats$Region<- factor(popsumstats$Region, levels = c("VanIsland", "LowerMain", "HaidaGwai","Northwest"))
popsumstats<- droplevels(popsumstats)

## box plot
ggplot(popsumstats, aes(x=Region, y=Pi, fill=Region)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi)")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))





### plots

# get mean and sd
sumdata<-group_by(popsumstats, Region) %>%
  summarise(mean=mean(Pi), sd=sd(Pi))


sumdatadf<-as.data.frame(sumdata)

sumdatadf$PihighSD<-sumdatadf$mean+sumdatadf$sd

sumdatadf$PilowSD<-sumdatadf$mean-sumdatadf$sd





## reorder
sumdatadf$Region<- factor(sumdatadf$Region, levels = c("VanIsland", "LowerMain", "HaidaGwai", "Northwest"))
sumdatadf<- droplevels(sumdatadf)

write.csv(sumdatadf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_4clusters_1370INDIV_165731SNPS.csv")

## plot - no legend
p<-ggplot(sumdatadf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Pihigh, ymin = Pilow),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi) ± stdev")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p


png("F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_4clusters_1370INDIV_165731SNPS.png", width = 8, height = 6.5, units = 'in', res = 600)
p<-ggplot(sumdatadf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Pihigh, ymin = Pilow),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi) ± stdev")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()



ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_4clusters_1370INDIV_165731SNPS_ggsave.png",  bg = "transparent",width =177.8, height = 160, units = c("mm"))
p<-ggplot(sumdatadf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 4)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  
  geom_errorbar(aes(ymax = Pihigh, ymin = Pilow),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi) ± stdev")+theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))
p<-p+theme(legend.position = "none")
p
dev.off()





# with SE + 95% CI ###########
sumdatat_SE<-group_by(popsumstats, Region) %>%
  summarise_each(funs(mean=mean(Pi),n=n(),sd=sd(Pi),se=sd(.)/sqrt(n())),Pi)



sumdatat_SEdf<-as.data.frame(sumdatat_SE)

sumdatat_SEdf$PihighSE<-sumdatat_SEdf$mean+sumdatat_SEdf$se

sumdatat_SEdf$PilowSE<-sumdatat_SEdf$mean-sumdatat_SEdf$se


sumdatat_SEdf$Pi_highCI<-sumdatat_SEdf$mean+1.96*(sumdatat_SEdf$sd/sqrt(sumdatat_SEdf$n))
sumdatat_SEdf$Pi_lowCI<-sumdatat_SEdf$mean-1.96*(sumdatat_SEdf$sd/sqrt(sumdatat_SEdf$n))


## reorder
sumdatat_SEdf$Region<- factor(sumdatat_SEdf$Region, levels = c("VanIsland", "LowerMain", "HaidaGwai", "Northwest"))
sumdatat_SEdf<- droplevels(sumdatat_SEdf)


write.csv(sumdatat_SEdf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_SE_4clusters_CI_1370INDIV_165731SNPS.csv")




## NEw COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS//Pi/Pi_mean_95pctCI_4clusters_1370INDIV_165731SNPS_neworange_correctCI.pdf",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(sumdatat_SEdf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#E69F00","#56B4E9","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  #ylim(0,0.0008)+
  geom_errorbar(aes(ymax = Pi_highCI, ymin = Pi_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi)")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()


ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS//Pi/Pi_mean_95pctCI_4clusters_1370INDIV_165731SNPS_neworange_correctCI_xlim.png",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(sumdatat_SEdf, aes(x = Region, color = Region, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#ab4e03", "#fca45d","#009E73","pink"),name = "Region", labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  scale_x_discrete(labels=c("Vancouver Island", "Lower mainland", "Haida Gwaii", "Northwest BC"))+
  ylim(0,0.0008)+
  geom_errorbar(aes(ymax = Pi_highCI, ymin = Pi_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi)")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()




########### two clusters #######


str(popsumstats)
str(popsumstats$twoclusters)

popsumstats$twoclusters<-as.factor(popsumstats$twoclusters)



## reorder
popsumstats$twoclusters<- factor(popsumstats$twoclusters, levels = c("CoastalBC", "HaidaGwai"))
popsumstats<- droplevels(popsumstats)

## box plot
ggplot(popsumstats, aes(x=twoclusters, y=Pi, fill=twoclusters)) +
  #stat_boxplot(geom ='errorbar') +  # horizontal whiskers
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c( "#FC8C3C","#009E73"),name = "Region", labels=c("CoastalBC", "HaidaGwai"))+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi)")+
  theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 16), title = element_text(size = 16))+
  scale_x_discrete(labels=c("Southwest BC", "Haida Gwaii"))





### plots

# get mean and sd
sumdata<-group_by(popsumstats, twoclusters) %>%
  summarise(mean=mean(Pi), sd=sd(Pi))


sumdatadf<-as.data.frame(sumdata)

sumdatadf$PihighSD<-sumdatadf$mean+sumdatadf$sd

sumdatadf$PilowSD<-sumdatadf$mean-sumdatadf$sd





## reorder
sumdatadf$twoclusters<- factor(sumdatadf$twoclusters, levels = c("CoastalBC", "HaidaGwai"))
sumdatadf<- droplevels(sumdatadf)

write.csv(sumdatadf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_4clusters_1370INDIV_165731SNPS_twoclusters.csv")


# with SE + 95% CI ###########
sumdatat_SE<-group_by(popsumstats, twoclusters) %>%
  summarise_each(funs(mean=mean(Pi),n=n(),sd=sd(Pi),se=sd(.)/sqrt(n())),Pi)



sumdatat_SEdf<-as.data.frame(sumdatat_SE)

sumdatat_SEdf$PihighSE<-sumdatat_SEdf$mean+sumdatat_SEdf$se

sumdatat_SEdf$PilowSE<-sumdatat_SEdf$mean-sumdatat_SEdf$se


sumdatat_SEdf$Pi_highCI<-sumdatat_SEdf$mean+1.96*(sumdatat_SEdf$sd/sqrt(sumdatat_SEdf$n))
sumdatat_SEdf$Pi_lowCI<-sumdatat_SEdf$mean-1.96*(sumdatat_SEdf$sd/sqrt(sumdatat_SEdf$n))


## reorder
sumdatat_SEdf$twoclusters<- factor(sumdatat_SEdf$twoclusters, levels = c("CoastalBC", "HaidaGwai"))
sumdatat_SEdf<- droplevels(sumdatat_SEdf)


write.csv(sumdatat_SEdf,"F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Pi/Pi_mean_sd_SE_4clusters_CI_1370INDIV_165731SNPS_twoclusters.csv")




## NEw COL
ggsave(file="F:/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS//Pi/Pi_mean_95pctCI_4clusters_1370INDIV_165731SNPS_neworange_correctCI_twoclusters.pdf",  bg = "transparent",width =200, height = 180, units = c("mm"))
p<-ggplot(sumdatat_SEdf, aes(x = twoclusters, color = twoclusters, y=mean)) +
  geom_point(size = 8)+
  theme_classic() +
  scale_color_manual(values=c( "#FC8C3C","#009E73"),name = "Region", labels=c("Coastal BC", "Haida Gwaii"))+
  scale_x_discrete(labels=c("Coastal BC", "Haida Gwaii"))+
  
  geom_errorbar(aes(ymax = Pi_highCI, ymin = Pi_lowCI),
                position = "dodge",
                width=0.2,
                size=0.5)+
  xlab("Region") + 
  ylab("Nucleotide diversity (Pi)")+theme(axis.title = element_text(face = "bold",size=20), axis.text = element_text(size = 18))
p<-p+theme(legend.position = "none",axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
p
dev.off()



