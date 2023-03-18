



## loci names I want to keep
gl.toad.nosibs.hwe.FIS@loc.names

locinamnes_371_732<-gl.toad.nosibs.hwe.FIS@loc.names






### object I need to remove them from "ID column"
filtered.VCF@fix[,'ID']
filtered.VCF@fix[,3]


#### the way to do it!! either way works
################

subset(filtered.VCF, filtered.VCF@fix[,3]%in%locinamnes_371_732)


filtered.VCF_371_732<-filtered.VCF[filtered.VCF@fix[,3]%in%locinamnes_371_732,]
###############




fix_filtered.VCF<-filtered.VCF@fix
fix_filtered.VCF@ID




for(i in 1:length(locinamnes_371_732)){
  
  filtered.VCF_test <- gl.toad_test[indNames(filtered.VCF) != listsibs[i]]
  
}


filtered.VCF_test@fix#ID




filtered.VCF_test2@fix<-subset(filtered.VCF@fix, filtered.VCF@fix[,3] == "17:17")


subset(filtered.VCF, filtered.VCF@fix[,3] == "17:17")




subset(filtered.VCF, filtered.VCF@fix[,3] == "17:17")



for(i in 1:length(locinamnes_371_732)){
  
  filtered.VCF_test <-subset(filtered.VCF@fix, filtered.VCF@fix[,3] = locinamnes_371_732[i])
  
}
