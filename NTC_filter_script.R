# Removing the effect of NTCs

# New method as discussed on Feb 3rd lab meeting.

setwd("/data/Users/kmccauley/WaldStudy/MakeOTUtable")
old.OTU <- read.table("otu_table_nochimera_filtered_0_00001.txt",header=TRUE, check.names=F, comment="", skip=1,sep="\t")
#First, identify the columns that are negative controls, and keep them in a new data frame? Or just define the column numbers?

tax.data <- old.OTU[,c("#OTU ID","taxonomy")]

row.names(old.OTU) <- old.OTU[,c("#OTU ID")]
old.OTU[,c("#OTU ID")] <- NULL

#Rename some negatives first (to comply with some logic that all negatives start with a letter (not a number).
names(old.OTU)[names(old.OTU) %in% c("24.35.")] <- "neg.24.35"
names(old.OTU)[names(old.OTU) %in% c("1.22.")] <- "neg.1.22"

#Negatives were not named in a standardized way, so I had to do this logic, which was to take the first letter and ask whether or not it was anything that indicated a negative (ie, n=neg or w=water)
#Only change the last part of this code where the letters are identified. You can specify any number of letters, but they are case-sensitive and should be specific to ONLY NTCs
neg.names <- names(old.OTU)[substr(names(old.OTU),1,1) %in% c("n","N","w")]

neg.OTUtable <- old.OTU[,names(old.OTU) %in% c(neg.names,"taxonomy")]
samp.OTUtable <- old.OTU[,!names(old.OTU) %in% c(neg.names,"taxonomy")]
samp.OTUtable$taxonomy <- NULL

# I first believe that a number of my negative controls aren't "real" negative controls in my data, and therefore can't be trusted as such, so I will identify them and remove them from the "negative OTU" table. Wouldn't always be the case for all studies, but here, the mean works quite nicely to remove the negatives that are not "true" negatives.

names(neg.OTUtable)[colSums(neg.OTUtable[,1:(ncol(neg.OTUtable)-1)]) > mean(colSums(neg.OTUtable[,1:(ncol(neg.OTUtable)-1)]))]

untrue.negs <- names(neg.OTUtable)[colSums(neg.OTUtable[,1:(ncol(neg.OTUtable)-1)]) > mean(colSums(neg.OTUtable[,1:(ncol(neg.OTUtable)-1)]))]
untrue.negs <- untrue.negs[untrue.negs != "taxonomy"]

neg.OTUtable2 <- neg.OTUtable[,!names(neg.OTUtable) %in% untrue.negs]

#The first step of the current plan is to remove the OTUs that are in greater than 50% of the NTCs
#To do this, I need to specify the number of NTCs that each OTU is present in.

num.NTCs.perOTU <- rowSums(neg.OTUtable2[,1:(ncol(neg.OTUtable2)-1)] > 0) > max(rowSums(neg.OTUtable2[,1:(ncol(neg.OTUtable2)-1)] > 0))/2
neg.OTUtable2[num.NTCs.perOTU,]

#Move to relabund temporarily?
#neg.relabund <- apply(neg.OTUtable2[,1:(ncol(neg.OTUtable2)-1)], 2, function(x) x/sum(x))
#rowSums(neg.relabund)
#Some samples have 0,1, or 2 reads in them
#samp.OTUtable.filt <- samp.OTUtable[,colSums(samp.OTUtable) >10]
#samp.relabund <- apply(samp.OTUtable.filt[,1:(ncol(samp.OTUtable.filt)-1)],2, function(x) x/sum(x))
#rowSums(samp.relabund)
#new.samp <- samp.OTUtable[!rownames(samp.OTUtable) %in% rownames(neg.OTUtable2[rowSums(neg.relabund) > rowSums(samp.relabund),]),]
#final.OTUtable.tax <- merge(new.samp,tax.data, by.x=0,by.y=1)
#names(final.OTUtable.tax)[names(final.OTUtable.tax) %in% "Row.names"] <- "#OTU ID"

#write.table(final.OTUtable.tax, "proportion_based_neg_filter.txt", row.names=FALSE, sep="\t", quote=F)

#Initial filtering, now that I've determined which rows to keep in a LOGICAL vector.
samp.OTUtable2 <- samp.OTUtable[!num.NTCs.perOTU,]

#The second step is to do subtraction 
#get a vector of the maximum in the remaining negatives (neg.OTUtable2)
neg.OTUtable2$taxonomy <- NULL
max.in.NTC <- apply(neg.OTUtable2,1, max)
summary(max.in.NTC)

#Now remove the count from the actual samples.
#test <- t(apply(samp.OTUtable2, 1, function(x) x - max.in.NTC))

#Clunky for-loop
samp.OTUtable3 <- samp.OTUtable2
for(i in 1:nrow(samp.OTUtable2)) {
   print(i)
   name <- names(max.in.NTC[i])
   samp.OTUtable3[name,] <- samp.OTUtable2[name,] - max.in.NTC[name]
   samp.OTUtable3[name,][samp.OTUtable3[name,] < 0] <- 0
}

final.OTUtable <- samp.OTUtable3[rowSums(samp.OTUtable3) != 0,]

final.OTUtable.tax <- merge(final.OTUtable,tax.data, by.x=0,by.y=1)
names(final.OTUtable.tax)[names(final.OTUtable.tax) %in% "Row.names"] <- "#OTU ID"
write.table(final.OTUtable.tax, "otu_table_nochimera_filtered_0_00001_noneg.txt", row.names=FALSE, sep="\t",quote=F)
