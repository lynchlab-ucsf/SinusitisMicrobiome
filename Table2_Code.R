# Table of Variables Associated with Community Composition

## May choose to focus solely on one distance matrix, so supplemental tables may be generated for the other DMs

## Sue denotes Table 1 as being the LME table. Perhaps the bootstrapped analyses will all go into supplemental?

library(geepack)
rm(list=ls())

myDat <- sample_data(readRDS("/data/Users/kmccauley/WaldStudy/AnalysisPhyloseq.rds"))
myDat$row.names <- NULL
myDat$any.virus <- !myDat$hrv_type %in% "None"
myDat$log.SpnQPCR <- log(myDat$SpnQPCR.CFUe.ml+1)
myDat$log.McQPCR <- log(myDat$McQPCR.CFUe.ml+1)
myDat$log.HfluQPCR <- log(myDat$HfluQPCR.CFUe.ml+1)
myDat$site <- factor(myDat$site, labels=c("Site 1","Site 2"))
myDat$healthy_sinusitis_uri <- as.character(myDat$sinusitis_vs_uri)
myDat$healthy_sinusitis_uri[is.na(myDat$healthy_sinusitis_uri)] <- "Healthy"
myDat$healthy_sinusitis_uri <- factor(myDat$healthy_sinusitis_uri)
myDat$Surv_Sick_Recov <- factor(myDat$Surv_Sick_Recov, levels=c("Surveillance","Sick","Recovery"))

levels(myDat$Mother.Education)[levels(myDat$Mother.Education) %in% c("Grade School","High School","Vo/Tech")] <- "Other"


tnames <- NULL
names <- strsplit(as.character(myDat$dom_gen),";")
tnames <- t(matrix(unlist(names),nrow=6))
tnames <- substr(tnames,4,30)
tnames[,6][tnames[,6] == ""] <- tnames[,5][tnames[,6] == ""]
tnames[,6][tnames[,6] == ""] <- tnames[,4][tnames[,6] == ""]
colnames(tnames) <- c("domgen_kingdom","domgen_phylum","domgen_class","domgen_order","domgen_family","domgen_genus")
myDat <- cbind(myDat, tnames[,c(5,6)])
tnames <- NULL
names <- strsplit(as.character(myDat$dom_fam),";")
tnames <- t(matrix(unlist(names),nrow=5))
tnames <- substr(tnames,4,30)
tnames[,5][tnames[,5] == ""] <- tnames[,4][tnames[,5] == ""]
tnames[,5][tnames[,5] == ""] <- tnames[,4][tnames[,5] == ""]
colnames(tnames) <- c("domfam_kingdom","domfam_phylum","domfam_class","domfam_order","domfam_family")
myDat <- cbind(myDat, tnames[,c(4,5)])

levels(myDat$domgen_genus)[levels(myDat$domgen_genus) %in% names(table(myDat$domgen_genus)[table(myDat$domgen_genus) < 2])] <- "Other"
#levels(myDat_test$domfam_family)[levels(myDat_test$domfam_family) %in% names(table(myDat_test$domfam_family)[table(myDat_test$domfam_family) < 2])] <- "Other"

myDat$domgen_genus <- relevel(myDat$domgen_genus, ref="Moraxella")

uwuf.dm <- read.table("bdiv_mr_drop/unweighted_unifrac_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
wuf.dm <- read.table("bdiv_mr_drop/weighted_unifrac_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
bray.dm <- read.table("bdiv_mr_drop/bray_curtis_dm.txt", header=T, sep="\t",row.names=1,check.names=F)
can.dm <- read.table("bdiv_mr_drop/canberra_dm.txt", header=T, sep="\t",row.names=1,check.names=F)

#varlink <- read.csv("/Volumes/data/Users/kmccauley/WaldStudy/Analysis/VarLinkageTable.csv")

pc1.lme <- function(dm, map, varlist, idvar="PatID", sampid=0) {
  pacman::p_load(ape, geepack)
  full.results <- NULL
  for(i in 1:length(varlist)) {
    resp <- varlist[i]
    print(resp)
    #If the values are all the same across all samples or all different, don't run LME on that variable.
    if(length(unique(map[,resp][!is.na(map[,resp])])) > 1 & length(unique(map[,resp])) < length(map[,resp])) {
    map.sub <- map[!is.na(map[,resp]),]
    if(sampid == 0) { 
      dm.sub <- dm[rownames(map.sub), rownames(map.sub)] 
    } else { 
      dm.sub <- dm[as.character(map.sub[,sampid]), as.character(map.sub[,sampid])]
      }
    pcs <- pcoa(dm.sub)
    PC1 <- pcs$vectors[,1:3]
    print(pcs$values$Relative_eig[1:3])
    if(sampid == 0) {
      merge.mapA <- merge(map.sub, PC1, by=0)
      if(dim(merge.mapA)[1] ==0 | dim(merge.mapA)[2] == 0) {
        merge.mapA <- merge(map.sub, PC1, by.x=sampid, by.y=0)
      } 
    } else {
      merge.mapA <- merge(map.sub, PC1, by.x=sampid, by.y=0)
    }
    map.fin <- merge.mapA[order(merge.mapA[,idvar]),]
    gee.mod1 <- geeglm(Axis.1 ~ map.fin[,resp], id=map.fin[,idvar], corstr="exchangeable", data=map.fin)
    gee.mod2 <- geeglm(Axis.2 ~ map.fin[,resp], id=map.fin[,idvar], corstr="exchangeable", data=map.fin)
    gee.mod3 <- geeglm(Axis.3 ~ map.fin[,resp], id=map.fin[,idvar], corstr="exchangeable", data=map.fin)
    
    p.pc1 <- summary(gee.mod1)$coef[2,4]
    p.pc2 <- summary(gee.mod2)$coef[2,4]
    p.pc3 <- summary(gee.mod3)$coef[2,4]
    
    #I thought Sample Size in groups would be useful information to have, so I've added it here...
    #if((is.factor(map.fin[,resp]) == TRUE | is.integer(map.fin[,resp]) == TRUE) & length(levels(as.factor(map[,resp]))) < 15) {
    #  gettable <- table(map.fin[,resp])
    #  sample.size <- paste(gettable, collapse ="/")
    #} else {
    #  sample.size <- length(map.fin[,resp])
    #}
    #results <- cbind(resp, sample.size, p.value)
    results <- cbind(resp, p.pc1, p.pc2, p.pc3)
    full.results <- rbind(full.results, results)
    }}
  full.resultsB <- as.data.frame(full.results)
  full.resultsB[,2:ncol(full.resultsB)] <- lapply(full.resultsB[,2:ncol(full.resultsB)], function(x) as.numeric(as.character(x)))
  full.resultsB$fdr.p1 <- p.adjust(full.resultsB$p.pc1, method="fdr")
  full.resultsB$fdr.p2 <- p.adjust(full.resultsB$p.pc2, method="fdr")
  full.resultsB$fdr.p3 <- p.adjust(full.resultsB$p.pc3, method="fdr")
  full.resultsB$fdr.p1 <- ifelse(full.resultsB$fdr.p1 < 0.001, "<0.001", as.character(round(full.resultsB$fdr.p1,3)))
  full.resultsB$fdr.p2 <- ifelse(full.resultsB$fdr.p2 < 0.001, "<0.001", as.character(round(full.resultsB$fdr.p2,3)))
  full.resultsB$fdr.p3 <- ifelse(full.resultsB$fdr.p3 < 0.001, "<0.001", as.character(round(full.resultsB$fdr.p3,3)))
  
  #full.resultsB$p.value <- ifelse(full.resultsB$p.value < 0.001, "<0.001", as.character(round(full.resultsB$p.value,3)))
  
  #rownames(full.results) <- full.results$resp
  #full.results$resp <- NULL
  return(full.resultsB)
}

mynewlist <- names(myDat)[!names(myDat) %in% c("PatID","#SampleID", "Date","Parity","Symp_Init_Date","Symp_Start_Date","Symp_Day7_Date","Symp_Day10_Date","Symp_Day15_Date","Initial.Visit.Date", "Race_Other")]

vars_of_int <- c("chao1","PD_whole_tree","equitability","any.virus","StudySubGroup","Visit.Type_Lab","Gender","Race","Ethnicity","Mother.Age","Mother.Education","Number.of.Surveillance.visits","Size.of.Household","Other.Children_Initial","Share.room","Tobacco.Exposure_Initial","Dog.at.Home","Animal.Exposure_Dog.at.School_Initial","Cat.at.Home","Cat.at.School","Live.on.Farm","Daycare_Initial","age_init","healthy_vs_sick","sinusitis_vs_nonsin","sinusitis_vs_uri","Surv_Sick_Recov","PatID","season","site","totURI","tot.virus.in.samp","age_init","hrv_type","domgen_genus", "healthy_sinusitis_uri", "log.SpnQPCR","log.HfluQPCR","log.McQPCR")

wuf.pc1 <- pc1.lme(dm=wuf.dm, map=myDat, varlist=vars_of_int, idvar="PatID", sampid=0)

uwuf.pc1 <- pc1.lme(dm=uwuf.dm, map=myDat, varlist=vars_of_int, idvar="PatID", sampid=0)

bray.pc1 <- pc1.lme(dm=bray.dm, map=myDat, varlist=vars_of_int, idvar="PatID", sampid=0)

can.pc1 <- pc1.lme(dm=can.dm, map=myDat, varlist=vars_of_int, idvar="PatID", sampid=0)

#sav <- cbind(wuf.pc1, uwuf.pc1, bray.pc1, can.pc1)


pcs <- pcoa(can.dm)
PC1 <- pcs$vectors[,1:3]
merge.mapA <- merge(myDat, PC1, by.x=sampid, by.y=0)
map.fin <- merge.mapA[order(merge.mapA[,idvar]),]
gee.mod1 <- geeglm(Axis.1 ~ healthy_sinusitis_uri + Tobacco.Exposure_Initial + Daycare_Initial + Race + season + Cat.at.Home + Dog.at.Home + Mother.Education, id=PatID, corstr="exchangeable", data=map.fin)
summary(gee.mod1)$coef

gee.mod2 <- geeglm(Axis.2 ~ healthy_sinusitis_uri + Tobacco.Exposure_Initial + Daycare_Initial + Race + season, id=PatID, corstr="exchangeable", data=map.fin)
summary(gee.mod2)$coef

gee.mod3 <- geeglm(Axis.3 ~ healthy_sinusitis_uri + Tobacco.Exposure_Initial + Daycare_Initial + Race + season, id=PatID, corstr="exchangeable", data=map.fin)
summary(gee.mod3)$coef

gee.mod1 <- geeglm(!healthy_sinusitis_uri %in% "Healthy" ~ season, id=PatID, corstr="exchangeable", data=map.fin)
summary(gee.mod1)
# To be useful later...
final.pc1 <- merge(varlink, wuf.pc1, by.x="Old.Name", by.y="resp", all=TRUE)
final.pc1$Old.Name <- NULL
row.names(final.pc1) <- NULL
final.pc1_clean <- final.pc1[!is.na(final.pc1$p.value),]

write.csv(x=final.pc1_clean, file="/Volumes/data/Users/kmccauley/WaldStudy/Analysis/PC1_Table_05Apr.csv")
