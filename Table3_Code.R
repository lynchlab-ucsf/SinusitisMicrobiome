# Beta Diversity Analysis

#Determine factors associated with community composition. Need to employ a few different methods in order to determine what is associated with community composition in this longitudinal dataset.

#Maybe some of these things that are coming up as non-significant might only be significant in the context of the greater group? I.e., I should find a better way than this bootstrap (maybe PC1?) to figure this out?

library(vegan)
library(phyloseq)

uwuf.dm <- read.table("/data/Users/kmccauley/WaldStudy/MakeOTUtable/bdiv_mr_drop/unweighted_unifrac_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
wuf.dm <- read.table("/data/Users/kmccauley/WaldStudy/MakeOTUtable/bdiv_mr_drop/weighted_unifrac_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
bray.dm <- read.table("/data/Users/kmccauley/WaldStudy/MakeOTUtable/bdiv_mr_drop/bray_curtis_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
can.dm <- read.table("/data/Users/kmccauley/WaldStudy/MakeOTUtable/bdiv_mr_drop/canberra_dm.txt", header=T, sep="\t",row.names=1,check.names=F)

map <- sample_data(readRDS("/data/Users/kmccauley/WaldStudy/AnalysisPhyloseq.rds"))
#map$Ethnicity <- map$Ethniciy
#map$Ethniciy <- NULL
map$log.SpnQPCR <- log(map$SpnQPCR.CFUe.ml+1)
map$log.HfluQPCR <- log(map$HfluQPCR.CFUe.ml+1)
map$log.McQPCR <- log(map$McQPCR.CFUe.ml+1)

set.seed(123)

bootstrap.adonis <- function(dm,resp,rep.tot=100) {
  R2 <- NULL
  p <- NULL
  pb <- txtProgressBar(min=0,max=rep.tot,style=3)
  for(i in 1:rep.tot) {
    #print(i)
    set.seed(123)
    rand.subset <- tapply(rownames(map),map$PatID,sample,size=1)
    r.subset.map <- map[row.names(map) %in% rand.subset & !is.na(map[,resp]),]
    r.subset.dm <- dm[row.names(r.subset.map),row.names(r.subset.map)]
    my.adonis <- adonis(as.dist(r.subset.dm) ~ r.subset.map[[resp]])
    R2[i] <- my.adonis$aov.tab$R2[1]
    p[i] <- my.adonis$aov.tab[1,6]
    setTxtProgressBar(pb, i)
  }
  return(cbind(resp,mean(R2),mean(p)))
  close(pb)
}

map$any.virus <- map$tot.virus.in.samp > 0

vars_of_int <- c("any.virus","season","StudySubGroup","Visit.Type_Lab","Gender","Race","Ethnicity","Mother.Age","Mother.Education","Number.of.Surveillance.visits","Size.of.Household","Other.Children_Initial","Share.room","Tobacco.Exposure_Initial","Dog.at.Home","Animal.Exposure_Dog.at.School_Initial","Cat.at.Home","Cat.at.School","Live.on.Farm","Daycare_Initial","age_init","dom_fam","dom_gen","healthy_vs_sick","sinusitis_vs_nonsin","sinusitis_vs_uri","Surv_Sick_Recov","PatID","season","site","totURI","tot.virus.in.samp","age_init","hrv_type","log.SpnQPCR","log.HfluQPCR","log.McQPCR")

fin <- NULL
for(i in vars_of_int) {
print(i)
wuf.run <- data.frame(bootstrap.adonis(wuf.dm, i))
names(wuf.run) <- c("Variable","Weighted R2","P-value")
uwuf.run <- data.frame(bootstrap.adonis(uwuf.dm, i))
names(uwuf.run) <- c("Variable","Unweighted R2","P-value")
can.run <- data.frame(bootstrap.adonis(can.dm, i))
names(can.run) <- c("Variable","Canberra R2","P-value")
bray.run <- data.frame(bootstrap.adonis(bray.dm, i))
names(bray.run) <- c("Variable","Bray R2","P-value")
fin <- rbind(fin, cbind(wuf.run, uwuf.run, can.run, bray.run))
}
fin
write.table(fin,"/data/Users/kmccauley/WaldStudy/Analysis/BetaDiversityAnalysis/BDiv_Bootstrap_Results.txt",sep="\t",quote=F, row.names=F)

