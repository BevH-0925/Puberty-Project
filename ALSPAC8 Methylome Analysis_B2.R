setwd("~/Documents/ALSPAC Methylome")

library(ChAMP)
library(randomForest)
library(Boruta)
library(pROC)
library(ROCR)
library(smotefamily)
library(caret)
library(mltools)
library(ggplot2)
library(tidyverse)

#Data loading and set up for chAMP analysis----
load("~/Documents/ALSPAC Methylome/Data/samplesheet/samplesheet.Robj")
load("~/Documents/ALSPAC Methylome/Data/detection_p_values/data.Robj")
load("~/Documents/ALSPAC Methylome/Data/betas/data.Robj")

#Remove duplicate rows as per ALSPAC methylome user guide
Dups<-which(samplesheet$duplicate.rm=="Remove")
samplesheet<-samplesheet[-Dups,]

#Change row name to participant ID plus rank (combine columns 2 & 3, no spaces), as for Clin Data file
rownames(samplesheet)<-paste(samplesheet[,2],samplesheet[,3],sep="")

B2QMGp <- read.csv("~/Documents/ALSPAC Methylome/Data/B2QMGp.csv")
rownames(B2QMGp) <- B2QMGp$X
B2QMGp <- B2QMGp[,-1]

#Closest to B2 vs rest
ClosesttB2<-B2QMGp[,c(1,2)]
ClosesttB2$Sample_Group[ClosesttB2$DmB2Quartile=="1st Quartile"] <- "FirstQuartile"
ClosesttB2$Sample_Group[ClosesttB2$DmB2Quartile=="2nd Quartile"] <- "SecondtoFourthQuartiles"
ClosesttB2$Sample_Group[ClosesttB2$DmB2Quartile=="3rd Quartile"] <- "SecondtoFourthQuartiles"
ClosesttB2$Sample_Group[ClosesttB2$DmB2Quartile=="4th Quartile"] <- "SecondtoFourthQuartiles"

#Furthest from B2 vs rest
FurthestfB2<- B2QMGp[,c(1,2)]
FurthestfB2$Sample_Group[FurthestfB2$DmB2Quartile=="4th Quartile"] <- "FourthQuartile"
FurthestfB2$Sample_Group[FurthestfB2$DmB2Quartile=="1st Quartile"] <- "FirsttoThirdQuartiles"
FurthestfB2$Sample_Group[FurthestfB2$DmB2Quartile=="2nd Quartile"] <- "FirsttoThirdQuartiles"
FurthestfB2$Sample_Group[FurthestfB2$DmB2Quartile=="3rd Quartile"] <- "FirsttoThirdQuartiles"

#Subset of samplesheet for participants with time to B2
samplesheetB2 <- samplesheet[match(rownames(B2QMGp),rownames(samplesheet)),]

#Subset of methylome data to include only those participants with a time to B2
betasB2 <- betas[,samplesheetB2$Sample_Name]

#Subset of the detection.p file to include participants with a time to B2
detectionB2 <- detection.p[,samplesheetB2$Sample_Name]
                          
#As rownames need to match exactly between betas & detection.p, subset
detectionB2<- detectionB2[(rownames(betasB2)),]

#Add puberty groups to sample sheet
samplesheetclosest <- samplesheetB2
samplesheetclosest$Sample_Group <- ClosesttB2[,3]

samplesheetfurthest <- samplesheetB2
samplesheetfurthest$Sample_Group <- FurthestfB2[,3]

rm(betas)
rm(detection.p)
rm(samplesheet)

#chAMP analysis of closest to B2 versus rest----
myLoadc <- champ.filter(beta=betasB2,pd=samplesheetclosest,detP=detectionB2,arraytype="450K")
mynormc <- champ.norm(beta=myLoadc$beta,arraytype="450K",cores = 8)
myDMPc <- champ.DMP(beta=mynormc,pheno=myLoadc$pd$Sample_Group,adjPVal = 1)
myDMRc <- champ.DMR(beta=mynormc,pheno=myLoadc$pd$Sample_Group, method = "Bumphunter",adjPvalDmr=1)

#Transpose normalised data so individuals in rows and probes in columns
B2normc <- (t(mynormc))
B2normc <- as.data.frame(B2normc)
rownames(B2normc) <- paste(rownames(samplesheetB2))
save(B2normc, file= "B2normc.robj")

#Make new data object (data frame) from first part of list
DMPc <- myDMPc[[1]]
DMRc <- myDMRc[[1]]

#Subset of DMP list for top 500 (list is in order of p-value)
DMPctop500 <- DMPc[c(1:500),]

write.csv(DMPc,"DMPcB2.csv")
write.csv(DMPctop500,"DMPcB2top500.csv")
write.csv(DMRc,"DMRcB2.csv")

#Subset of mynorm file for top 500 DMPs
B2normc_top500 <- B2normc[,(rownames(DMPctop500))]
write.csv(B2normc_top500,"B2normc_top500.csv")

#chAMP analysis of furthest from B2 versus rest----
myLoadf <- champ.filter(beta=betasB2,pd=samplesheetfurthest,detP=detectionB2,arraytype="450K")
mynormf <- champ.norm(beta=myLoadf$beta,arraytype="450K",cores = 8)
myDMPf <- champ.DMP(beta=mynormf,pheno=myLoadf$pd$Sample_Group,adjPVal = 1)
myDMRf <- champ.DMR(beta=mynormf,pheno=myLoadf$pd$Sample_Group, method = "Bumphunter",adjPvalDmr=1)

#Transpose normalised data so individuals in rows and probes in columns
B2normf <- (t(mynormf))
B2normf <- as.data.frame(B2normf)
rownames(B2normf) <- paste(rownames(samplesheetB2))
save(B2normf, file= "B2normf.robj")

#Make new data object (data frame) from first part of list
DMPf <- myDMPf[[1]]
DMRf <- myDMRf[[1]]

#Subset of DMP list for top 500
DMPftop500 <- DMPf[c(1:500),]

write.csv(DMPf,"DMPfB2.csv")
write.csv(DMPftop500,"DMPfB2top500.csv")
write.csv(DMRf,"DMRfB2.csv")

#Subset of mynorm file for top 500 DMPs
B2normf_top500 <- B2normf[,(rownames(DMPftop500))]
write.csv(B2normf_top500,"B2normf_top500.csv")

