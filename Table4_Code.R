#Association of Clinical Variables with States (using geeglm)

library(geepack)
rm(list=ls())
myDat <- sample_data(readRDS("/data/Users/kmccauley/WaldStudy/AnalysisPhyloseq.rds"))
myDat$any.virus <- !myDat$hrv_type %in% "None"
myDat$log.SpnQPCR <- log(myDat$SpnQPCR.CFUe.ml+1)
myDat$log.McQPCR <- log(myDat$McQPCR.CFUe.ml+1)
myDat$log.HfluQPCR <- log(myDat$HfluQPCR.CFUe.ml+1)

levels(myDat$Mother.Education)[levels(myDat$Mother.Education) %in% c("Grade School","High School","Vo/Tech")] <- "Other"

myDat.Clust <- myDat
names(myDat.Clust)[names(myDat.Clust) %in% "Row.names"] <- "#SampleID"
myDat.Clust <- myDat.Clust[order(myDat.Clust$PatID),]

#Use clusters as outcomes:
myDat.Clust$C1 <- myDat.Clust$CanberraClusters == "I"
myDat.Clust$C2 <- myDat.Clust$CanberraClusters == "II"
myDat.Clust$C3 <- myDat.Clust$CanberraClusters == "III"
myDat.Clust$C4 <- myDat.Clust$CanberraClusters == "IV"

vars_of_int <- c("healthy_sinusitis_uri","Visit.Type_Lab","Surv_Sick_Recov","season", "viral_pos","tot.virus.in.samp","hrv_type")



#Change reference groups here:
myDat.Clust$hrv_type <- relevel(myDat.Clust$hrv_type, ref="None")
myDat.Clust$viral_pos <- myDat.Clust$tot.virus.in.samp > 0
myDat.Clust$StudySubGroup <- relevel(myDat.Clust$StudySubGroup, ref="URI0")
myDat.Clust$Visit.Type_Lab <- relevel(myDat.Clust$Visit.Type_Lab, ref="Surveillance/Entry")
myDat.Clust$Race <- relevel(myDat.Clust$Race, ref="White or Caucasian")
myDat.Clust$Ethnicity <- relevel(myDat.Clust$Ethnicity, ref="Non-Hispanic")
myDat.Clust$site <- factor(myDat.Clust$site, labels=c("Site 1","Site 2"))
myDat.Clust$healthy_sinusitis_uri <- as.character(myDat.Clust$sinusitis_vs_uri)
myDat.Clust$healthy_sinusitis_uri[is.na(myDat.Clust$healthy_sinusitis_uri)] <- "Healthy"
myDat.Clust$healthy_sinusitis_uri <- factor(myDat.Clust$healthy_sinusitis_uri)
myDat.Clust$Surv_Sick_Recov <- factor(myDat.Clust$Surv_Sick_Recov, levels=c("Surveillance","Sick","Recovery"))
#I want to make a nice dominant family and dominant genus variable
tnames <- NULL
names <- strsplit(as.character(myDat.Clust$dom_gen),";")
tnames <- t(matrix(unlist(names),nrow=6))
tnames <- substr(tnames,4,30)
tnames[,6][tnames[,6] == ""] <- tnames[,5][tnames[,6] == ""]
tnames[,6][tnames[,6] == ""] <- tnames[,4][tnames[,6] == ""]
colnames(tnames) <- c("domgen_kingdom","domgen_phylum","domgen_class","domgen_order","domgen_family","domgen_genus")
myDat.Clust <- cbind(myDat.Clust, tnames[,c(5,6)])
tnames <- NULL
names <- strsplit(as.character(myDat.Clust$dom_fam),";")
tnames <- t(matrix(unlist(names),nrow=5))
tnames <- substr(tnames,4,30)
tnames[,5][tnames[,5] == ""] <- tnames[,4][tnames[,5] == ""]
tnames[,5][tnames[,5] == ""] <- tnames[,4][tnames[,5] == ""]
colnames(tnames) <- c("domfam_kingdom","domfam_phylum","domfam_class","domfam_order","domfam_family")
myDat.Clust <- cbind(myDat.Clust, tnames[,c(4,5)])

levels(myDat.Clust$domgen_genus)[levels(myDat.Clust$domgen_genus) %in% names(table(myDat.Clust$domgen_genus)[table(myDat.Clust$domgen_genus) < 2])] <- "Other"
#levels(myDat.Clust_test$domfam_family)[levels(myDat.Clust_test$domfam_family) %in% names(table(myDat.Clust_test$domfam_family)[table(myDat.Clust_test$domfam_family) < 2])] <- "Other"

myDat.Clust$domgen_genus <- relevel(myDat.Clust$domgen_genus, ref="Moraxella")

full_results <- NULL

clusts <- c("C1","C2","C3","C4")
for(i in vars_of_int) {
  oneresult <- NULL
  clust.res <- NULL
  for(j in clusts) {
    sub.dat <- myDat.Clust[!is.na(myDat.Clust[,i]),]
    clust.formula <- as.formula(paste(j, "~", i))
    run.model <- geeglm(clust.formula, data=sub.dat, id=PatID, corstr="exchangeable")
    clust.est <- c(1.000, round(exp(summary(run.model)$coef[-(1),1]), 3))
    clust.p <-  c("--", round(summary(run.model)$coef[-1,4], 3))
    oneresult <- cbind(clust.est, clust.p)
    colnames(oneresult) <- paste0(colnames(oneresult), j)
    clust.res <- cbind(clust.res, oneresult)
  }
  level.names <- c(levels(factor(sub.dat[,i]))[1], rownames(summary(run.model)$coef)[-1])
  clustres.l <- cbind(level.names, clust.res)
  full_results <- rbind(full_results, clustres.l)
}

write.csv(full_results, "/data/Users/kmccauley/WaldStudy/PLoS_Revision/Table4_Results.csv")
