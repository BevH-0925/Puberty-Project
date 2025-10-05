library(ChAMP)
#MF----
myLoad <- champ.load(directory = "~/Documents/Code/RMCH_MethArrayBatch1/Arrays1to46", arraytype="EPIC")
champ.QC()
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPIC")
QC.GUI(beta=myNorm, pheno = myLoad$pd$Sample_Group, arraytype = "EPIC")
champ.SVD(myNorm, pd = myLoad$pd)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC",adjPVal = 1)
DMP.GUI(DMP=myDMP[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Group)
myDMR <- champ.DMR(beta=myNorm, 
                   pheno=myLoad$pd$Sample_Group,
                   method="Bumphunter", 
                   arraytype = "EPIC", 
                   cores=10, adjPvalDmr=1)
DMR.GUI(DMR=myDMR, beta=myNorm, pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")

#Transpose normalised data so individuals in rows and probes in columns
MyNormt <- (t(myNorm))
MyNormt <- as.data.frame(MyNormt)
save(MyNormt, file= "MyNormt.robj")

#Make new data object (data frame) from first part of list
DMP <- myDMP[[1]]
DMR <- myDMR[[1]]

#Subset of DMP list for top 500 (list is in order of p-value)
DMPtop500 <- DMP[c(1:500),]

write.csv(DMP,"DMP.csv")
write.csv(DMPtop500,"DMPtop500.csv")
write.csv(DMR,"DMR.csv")

#F only----
Sample_sheetF <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/Arrays1to46Fonly/Sample_sheet.csv", colClasses = c("Slide" = "character"))
Sample_sheetF <- subset(Sample_sheetF, !(Sample_Name %in% c("B01642_12", "B01642_15", "B01642_18","B01642_22","B01642_23","B01642_28","B01642_30","B01642_35",
                                                            "B01642_39","B01642_42","B01642_45","B01642_48","B01642_50","B01642_53")))
write.csv(Sample_sheetF,file="Sample_sheet.csv")

myLoadF <- champ.load(directory = "~/Documents/Code/RMCH_MethArrayBatch1/Arrays1to46Fonly", arraytype="EPIC")

myNormF <- champ.norm(myLoadF$beta, method = "BMIQ", arraytype = "EPIC")
QC.GUI(beta=myNormF, pheno = myLoadF$pd$Sample_Group, arraytype = "EPIC")
champ.SVD(myNormF, pd = myLoadF$pd)
myDMPF <- champ.DMP(beta = myNormF,pheno=myLoadF$pd$Sample_Group, arraytype = "EPIC",adjPVal = 1)
DMP.GUI(DMP=myDMPF[[1]], beta=myNormF, pheno=myLoadF$pd$Sample_Group)
myDMRF <- champ.DMR(beta=myNormF, 
                   pheno=myLoadF$pd$Sample_Group,
                   method="Bumphunter", 
                   arraytype = "EPIC", 
                   cores=10, adjPvalDmr=1)
DMR.GUI(DMR=myDMRF, beta=myNormF, pheno=myLoadF$pd$Sample_Group, arraytype = "EPIC")

#Transpose normalised data so individuals in rows and probes in columns
MyNormFt <- (t(myNormF))
MyNormFt <- as.data.frame(MyNormFt)
save(MyNormFt, file= "MyNormFt.robj")

#Make new data object (data frame) from first part of list
DMPF <- myDMPF[[1]]
DMRF <- myDMRF[[1]]

#Subset of DMP list for top 500 (list is in order of p-value)
DMPtop500F <- DMPF[c(1:500),]

write.csv(DMPF,"DMPF.csv")
write.csv(DMPtop500F,"DMPtop500F.csv")
write.csv(DMRF,"DMRF.csv")