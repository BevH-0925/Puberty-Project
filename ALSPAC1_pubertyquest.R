setwd("~/Documents/Code")
library(ggplot2)

load("~/Documents/Code/ALSPACClinData.rdata")
Clindata<-ALSPACClinData

# In transcriptomic data files, child participants are identified as ID then birth rank e.g.2A so need 
# to name each row as ID then birth rank (combine columns 1 & 2, no spaces).
# sep="" means no spaces between the two values.
rownames(Clindata)<-paste(Clindata[,1],Clindata[,2],sep="")

#Markers of puberty: 
#Pubic hair stage in males (MPS)
#Pubic hair stage in females (FPS)
#Genitalia stage in males (GS)
#Breast stage in females (BS)
#Age at first period in females (MENS)

#Males
#Male Pubic Hair

#Subset of male pubic hair stage - all questionnaires to see spread of data
Male.Pubic.Hair <- Clindata[,c("pub155","pub255","pub355","pub455","pub555","pub655","pub755","pub855","pub955")]
apply(Male.Pubic.Hair,2,table)

#Replacing male pubic hair stage with NA if answer was "not sure" on questionnaire ("6")
Clindata$pub155[Clindata$pub155==6]<-NA
Clindata$pub255[Clindata$pub255==6]<-NA
Clindata$pub355[Clindata$pub355==6]<-NA
Clindata$pub455[Clindata$pub455==6]<-NA
Clindata$pub555[Clindata$pub555==6]<-NA
Clindata$pub655[Clindata$pub655==6]<-NA
Clindata$pub755[Clindata$pub755==6]<-NA
Clindata$pub855[Clindata$pub855==6]<-NA
Clindata$pub955[Clindata$pub955==6]<-NA

#Check no values of 6 left for male pubic hair
Male.Pubic.Hair <- Clindata[,c("pub155","pub255","pub355","pub455","pub555","pub655","pub755","pub855","pub955")]
apply(Male.Pubic.Hair,2,table)

#Creating a new variable called MPS1 (NA for all rows)
Clindata$MPS1<-NA

#If male pubic hair stage from 1st questionnaire is 2 populate MPS1 with age in months at questionnaire completion otherwise use NA
Clindata$MPS1<-ifelse(Clindata$pub155==2,Clindata$pub195,NA)

#Repeated for next questionnaire:
Clindata$MPS2<-NA
Clindata$MPS2<-ifelse(Clindata$pub255==2,Clindata$pub295,NA)

#New variable called matched which gives the row numbers where MPS1 is not NA
matched<-which(!is.na(Clindata$MPS1))

#Calls MPS2 NA if already a value in for that row
Clindata$MPS2[matched]<-NA

#Next questionnaire:
Clindata$MPS3<-NA
Clindata$MPS3<-ifelse(Clindata$pub355==2,Clindata$pub397a,NA)

#New variable called combi which gives the row numbers where MPS2 is not NA
combi<-which(!is.na(Clindata$MPS2))

#matched combined with combi, so gives rows where value is earliest non-NA value
matched<-c(matched,combi)

#Calls MPS3 NA if already value in for that row
Clindata$MPS3[matched]<-NA

Clindata$MPS4<-NA
Clindata$MPS4<-ifelse(Clindata$pub455==2,Clindata$pub497a,NA)
combi<-which(!is.na(Clindata$MPS3))
matched<-c(matched,combi)
Clindata$MPS4[matched]<-NA

Clindata$MPS5<-NA
Clindata$MPS5<-ifelse(Clindata$pub555==2,Clindata$pub597a,NA)
combi<-which(!is.na(Clindata$MPS4))
matched<-c(matched,combi)
Clindata$MPS5[matched]<-NA

Clindata$MPS6<-NA
Clindata$MPS6<-ifelse(Clindata$pub655==2,Clindata$pub697a,NA)
combi<-which(!is.na(Clindata$MPS5))
matched<-c(matched,combi)
Clindata$MPS6[matched]<-NA

Clindata$MPS7<-NA
Clindata$MPS7<-ifelse(Clindata$pub755==2,Clindata$pub797a,NA)
combi<-which(!is.na(Clindata$MPS6))
matched<-c(matched,combi)
Clindata$MPS7[matched]<-NA

Clindata$MPS8<-NA
Clindata$MPS8<-ifelse(Clindata$pub855==2,Clindata$pub897a,NA)
combi<-which(!is.na(Clindata$MPS7))
matched<-c(matched,combi)
Clindata$MPS8[matched]<-NA

Clindata$MPS9<-NA
Clindata$MPS9<-ifelse(Clindata$pub955==2,Clindata$pub997a,NA)
combi<-which(!is.na(Clindata$MPS8))
matched<-c(matched,combi)
Clindata$MPS9[matched]<-NA

#Create a new column 'EAMPS2'. Adds up columns MPS1-MPS9, excluding missing values. 
#As all should be NA except the earliest one, EAMPS2=earliest age in months at male pubic hair stage 2
Clindata$EAMPS2<-rowSums(Clindata[,2231:2239],na.rm=TRUE)
Clindata$EAMPS2[Clindata$EAMPS2==0]<-NA
summary(Clindata$EAMPS2)

#Male Genitalia Stage

#Subset of male genitalia stage - all questionnaires to see spread of Clindata
Male.Genitalia.Stage <- Clindata[,c("pub150","pub250","pub350","pub450","pub550","pub650","pub750","pub850","pub950")]
apply(Male.Genitalia.Stage,2,table)

#Replacing male genitalia stage with NA if answer was "not sure" on questionnaire ("6")
Clindata$pub150[Clindata$pub150==6]<-NA
Clindata$pub250[Clindata$pub250==6]<-NA
Clindata$pub350[Clindata$pub350==6]<-NA
Clindata$pub450[Clindata$pub450==6]<-NA
Clindata$pub550[Clindata$pub550==6]<-NA
Clindata$pub650[Clindata$pub650==6]<-NA
Clindata$pub750[Clindata$pub750==6]<-NA
Clindata$pub850[Clindata$pub850==6]<-NA
Clindata$pub950[Clindata$pub950==6]<-NA

#Check no values of 6 left for male genitalia stage
Male.Genitalia.Stage <- Clindata[,c("pub150","pub250","pub350","pub450","pub550","pub650","pub750","pub850","pub950")]
apply(Male.Genitalia.Stage,2,table)

#Creating a new variable called GS1 (NA for all rows)
Clindata$GS1<-NA

#If male genitalia stage from 1st questionnaire is 2 populate GS1 with age in months at questionnaire completion otherwise use NA
Clindata$GS1<-ifelse(Clindata$pub150==2,Clindata$pub195,NA)

#Repeated for next questionnaire:
Clindata$GS2<-NA
Clindata$GS2<-ifelse(Clindata$pub250==2,Clindata$pub295,NA)

#New variable called matched which gives the row numbers where GS1 is not NA
matched<-which(!is.na(Clindata$GS1))

#Calls GS2 NA if already a value in for that row
Clindata$GS2[matched]<-NA

#Next questionnaire:
Clindata$GS3<-NA
Clindata$GS3<-ifelse(Clindata$pub350==2,Clindata$pub397a,NA)

#New variable called combi which gives the row numbers where GS2 is not NA
combi<-which(!is.na(Clindata$GS2))

#matched combined with combi, so gives rows where value is earliest non-NA value
matched<-c(matched,combi)

#Calls GS3 NA if already value in for that row
Clindata$GS3[matched]<-NA

Clindata$GS4<-NA
Clindata$GS4<-ifelse(Clindata$pub450==2,Clindata$pub497a,NA)
combi<-which(!is.na(Clindata$GS3))
matched<-c(matched,combi)
Clindata$GS4[matched]<-NA

Clindata$GS5<-NA
Clindata$GS5<-ifelse(Clindata$pub550==2,Clindata$pub597a,NA)
combi<-which(!is.na(Clindata$GS4))
matched<-c(matched,combi)
Clindata$GS5[matched]<-NA

Clindata$GS6<-NA
Clindata$GS6<-ifelse(Clindata$pub650==2,Clindata$pub697a,NA)
combi<-which(!is.na(Clindata$GS5))
matched<-c(matched,combi)
Clindata$GS6[matched]<-NA

Clindata$GS7<-NA
Clindata$GS7<-ifelse(Clindata$pub750==2,Clindata$pub797a,NA)
combi<-which(!is.na(Clindata$GS6))
matched<-c(matched,combi)
Clindata$GS7[matched]<-NA

Clindata$GS8<-NA
Clindata$GS8<-ifelse(Clindata$pub850==2,Clindata$pub897a,NA)
combi<-which(!is.na(Clindata$GS7))
matched<-c(matched,combi)
Clindata$GS8[matched]<-NA

Clindata$GS9<-NA
Clindata$GS9<-ifelse(Clindata$pub950==2,Clindata$pub997a,NA)
combi<-which(!is.na(Clindata$GS8))
matched<-c(matched,combi)
Clindata$GS9[matched]<-NA

#Create a new column 'EAGS2'. Adds up columns GS1-GS9, excluding missing values. 
#As all should be NA except the earliest one, EAGS2=earliest age in months at male genitalia stage 2.
Clindata$EAGS2<-rowSums(Clindata[,2241:2249],na.rm=TRUE)
Clindata$EAGS2[Clindata$EAGS2==0]<-NA
summary(Clindata$EAGS2)

#Female

#Breast stage

#Subset of female breast stage - all questionnaires to see spread of Clindata
Female.Breast.Stage <- Clindata[,c("pub130","pub230","pub330","pub430","pub530","pub630","pub730","pub830","pub930")]
apply(Female.Breast.Stage,2,table)

#Replacing female breast stage with NA if answer was "not sure" on questionnaire ("6")
Clindata$pub130[Clindata$pub130==6]<-NA
Clindata$pub230[Clindata$pub230==6]<-NA
Clindata$pub330[Clindata$pub330==6]<-NA
Clindata$pub430[Clindata$pub430==6]<-NA
Clindata$pub530[Clindata$pub530==6]<-NA
Clindata$pub630[Clindata$pub630==6]<-NA
Clindata$pub730[Clindata$pub730==6]<-NA
Clindata$pub830[Clindata$pub830==6]<-NA
Clindata$pub930[Clindata$pub930==6]<-NA

#Check no values of 6 left for female breast stage
Female.Breast.Stage <- Clindata[,c("pub130","pub230","pub330","pub430","pub530","pub630","pub730","pub830","pub930")]
apply(Female.Breast.Stage,2,table)

Clindata$BS1<-NA
Clindata$BS1<-ifelse(Clindata$pub130==2,Clindata$pub195,NA)

Clindata$BS2<-NA
Clindata$BS2<-ifelse(Clindata$pub230==2,Clindata$pub295,NA)
matched<-which(!is.na(Clindata$BS1))
Clindata$BS2[matched]<-NA

Clindata$BS3<-NA
Clindata$BS3<-ifelse(Clindata$pub330==2,Clindata$pub397a,NA)
combi<-which(!is.na(Clindata$BS2))
matched<-c(matched,combi)
Clindata$BS3[matched]<-NA

Clindata$BS4<-NA
Clindata$BS4<-ifelse(Clindata$pub430==2,Clindata$pub497a,NA)
combi<-which(!is.na(Clindata$BS3))
matched<-c(matched,combi)
Clindata$BS4[matched]<-NA

Clindata$BS5<-NA
Clindata$BS5<-ifelse(Clindata$pub530==2,Clindata$pub597a,NA)
combi<-which(!is.na(Clindata$BS4))
matched<-c(matched,combi)
Clindata$BS5[matched]<-NA

Clindata$BS6<-NA
Clindata$BS6<-ifelse(Clindata$pub630==2,Clindata$pub697a,NA)
combi<-which(!is.na(Clindata$BS5))
matched<-c(matched,combi)
Clindata$BS6[matched]<-NA

Clindata$BS7<-NA
Clindata$BS7<-ifelse(Clindata$pub730==2,Clindata$pub797a,NA)
combi<-which(!is.na(Clindata$BS6))
matched<-c(matched,combi)
Clindata$BS7[matched]<-NA

Clindata$BS8<-NA
Clindata$BS8<-ifelse(Clindata$pub830==2,Clindata$pub897a,NA)
combi<-which(!is.na(Clindata$BS7))
matched<-c(matched,combi)
Clindata$BS8[matched]<-NA

Clindata$BS9<-NA
Clindata$BS9<-ifelse(Clindata$pub930==2,Clindata$pub997a,NA)
combi<-which(!is.na(Clindata$BS8))
matched<-c(matched,combi)
Clindata$BS9[matched]<-NA

#Create a new column 'EABS2'. Adds up columns BS1-BS9, excluding missing values. 
#As all should be NA except the earliest one, EABS2=earliest age in months at breast stage 2.
Clindata$EABS2<-rowSums(Clindata[,2251:2259],na.rm=TRUE)
Clindata$EABS2[Clindata$EABS2==0]<-NA
summary(Clindata$EABS2)

#Female pubic hair

#Subset of female pubic hair stage - all questionnaires to see spread of Clindata
Female.Pubic.Hair <- Clindata[,c("pub135","pub235","pub335","pub435","pub535","pub635","pub735","pub835","pub935")]
apply(Female.Pubic.Hair,2,table)

#Replacing female public hair stage with NA if answer was "not sure" on questionnaire ("6")
Clindata$pub135[Clindata$pub135==6]<-NA
Clindata$pub235[Clindata$pub235==6]<-NA
Clindata$pub335[Clindata$pub335==6]<-NA
Clindata$pub435[Clindata$pub435==6]<-NA
Clindata$pub535[Clindata$pub535==6]<-NA
Clindata$pub635[Clindata$pub635==6]<-NA
Clindata$pub735[Clindata$pub735==6]<-NA
Clindata$pub835[Clindata$pub835==6]<-NA
Clindata$pub935[Clindata$pub935==6]<-NA

#To check no values of 6
Female.Pubic.Hair <- Clindata[,c("pub135","pub235","pub335","pub435","pub535","pub635","pub735","pub835","pub935")]
apply(Female.Pubic.Hair,2,table)

Clindata$FPS1<-NA
Clindata$FPS1<-ifelse(Clindata$pub135==2,Clindata$pub195,NA)

Clindata$FPS2<-NA
Clindata$FPS2<-ifelse(Clindata$pub235==2,Clindata$pub295,NA)
matched<-which(!is.na(Clindata$FPS1))
Clindata$FPS2[matched]<-NA

Clindata$FPS3<-NA
Clindata$FPS3<-ifelse(Clindata$pub335==2,Clindata$pub397a,NA)
combi<-which(!is.na(Clindata$FPS2))
matched<-c(matched,combi)
Clindata$FPS3[matched]<-NA

Clindata$FPS4<-NA
Clindata$FPS4<-ifelse(Clindata$pub435==2,Clindata$pub497a,NA)
combi<-which(!is.na(Clindata$FPS3))
matched<-c(matched,combi)
Clindata$FPS4[matched]<-NA

Clindata$FPS5<-NA
Clindata$FPS5<-ifelse(Clindata$pub535==2,Clindata$pub597a,NA)
combi<-which(!is.na(Clindata$FPS4))
matched<-c(matched,combi)
Clindata$FPS5[matched]<-NA

Clindata$FPS6<-NA
Clindata$FPS6<-ifelse(Clindata$pub635==2,Clindata$pub697a,NA)
combi<-which(!is.na(Clindata$FPS5))
matched<-c(matched,combi)
Clindata$FPS6[matched]<-NA

Clindata$FPS7<-NA
Clindata$FPS7<-ifelse(Clindata$pub735==2,Clindata$pub797a,NA)
combi<-which(!is.na(Clindata$FPS6))
matched<-c(matched,combi)
Clindata$FPS7[matched]<-NA

Clindata$FPS8<-NA
Clindata$FPS8<-ifelse(Clindata$pub835==2,Clindata$pub897a,NA)
combi<-which(!is.na(Clindata$FPS7))
matched<-c(matched,combi)
Clindata$FPS8[matched]<-NA

Clindata$FPS9<-NA
Clindata$FPS9<-ifelse(Clindata$pub935==2,Clindata$pub997a,NA)
combi<-which(!is.na(Clindata$FPS8))
matched<-c(matched,combi)
Clindata$FPS9[matched]<-NA

#Create a new column 'EAFPS2'. Adds up columns FPS1-FPS9, excluding missing values. 
#As all should be NA except the earliest one, EAFPS2=earliest age in months at female pubic hair stage 2
Clindata$EAFPS2<-rowSums(Clindata[,2261:2269],na.rm=TRUE)
Clindata$EAFPS2[Clindata$EAFPS2==0]<-NA
summary(Clindata$EAFPS2)

#Age at 1st period

#Subset of age at 1st period 
Age.1st.Period <- Clindata[,c("pub112","pub212","pub312","pub412","pub512","pub612","pub712","pub812","pub912")]
apply(Age.1st.Period,2,table)

#Take value of earliest age given for 1st period, rest NA
Clindata$MENS1<-Clindata$pub112
Clindata$MENS2<-Clindata$pub212
matched<-which(!is.na(Clindata$MENS1))
Clindata$MENS2[matched]<-NA

Clindata$MENS3<-Clindata$pub312
combi<-which(!is.na(Clindata$MENS2))
matched<-c(matched,combi)
Clindata$MENS3[matched]<-NA

Clindata$MENS4<-Clindata$pub412
combi<-which(!is.na(Clindata$MENS3))
matched<-c(matched,combi)
Clindata$MENS4[matched]<-NA

Clindata$MENS5<-Clindata$pub512
combi<-which(!is.na(Clindata$MENS4))
matched<-c(matched,combi)
Clindata$MENS5[matched]<-NA

Clindata$MENS6<-Clindata$pub612
combi<-which(!is.na(Clindata$MENS5))
matched<-c(matched,combi)
Clindata$MENS6[matched]<-NA

Clindata$MENS7<-Clindata$pub712
combi<-which(!is.na(Clindata$MENS6))
matched<-c(matched,combi)
Clindata$MENS7[matched]<-NA

Clindata$MENS8<-Clindata$pub812
combi<-which(!is.na(Clindata$MENS7))
matched<-c(matched,combi)
Clindata$MENS8[matched]<-NA

Clindata$MENS9<-Clindata$pub912
combi<-which(!is.na(Clindata$MENS8))
matched<-c(matched,combi)
Clindata$MENS9[matched]<-NA

#Create a new column 'AGMENS'. Adds up columns MENS1-MENS9, excluding missing values. 
#As all should be NA except the earliest one, AGMENS=earliest age in months given for 1st period.
Clindata$AGMENS<-rowSums(Clindata[,2271:2279],na.rm=TRUE)
Clindata$AGMENS[Clindata$AGMENS==0]<-NA
summary(Clindata$AGMENS)

#Add columns to calculate age in years
#Age at male pubic hair stage 2 in years:MPY
#Age at male genitalia stage 2 in years:MGY
#Age at breast stage 2 in years: BY
#Age at female pubic hair stage 2 in years: FPY
Clindata$MPY <- Clindata$EAMPS2/12
Clindata$MGY <- Clindata$EAGS2/12
Clindata$BY <- Clindata$EABS2/12
Clindata$FPY <- Clindata$EAFPS2/12
Clindata$MENY <- Clindata$AGMENS/12

#Each new variable plotted separately
#Plot of male pubic hair stage 2 only
ggplot(Clindata,aes(x=Clindata$MPY))+
  geom_histogram(color="white",fill = "blue", alpha=0.7)+
  labs(x = "Age in years", y= "Number of Participants", title ="Earliest Questionnaire Response of Male Pubic Hair Stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"))

#Plot of genitalia stage 2 only
ggplot(Clindata,aes(x=Clindata$MGY))+
  geom_histogram(color="white",fill = "red", alpha=0.7)+
  labs(x = "Age in years", y= "Number of Participants", title ="Earliest Questionnaire Response of Male Genitalia Stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"))

#Plot of breast stage 2 only
ggplot(Clindata,aes(x=Clindata$BY))+
  geom_histogram(color="white",fill = "blue", alpha=0.7)+
  labs(x = "Age in years", y= "Number of Participants", title ="Earliest Questionnaire Response of Breast Stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"))

#Plot of female pubic hair stage 2 only
ggplot(Clindata,aes(x=Clindata$FPY))+
  geom_histogram(color="white",fill = "red", alpha=0.7)+
  labs(x = "Age in years", y= "Number of Participants", title ="Earliest Questionnaire Response of Female Pubic Hair Stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"))

#Plot of age at menarche
ggplot(Clindata,aes(x=Clindata$MENY))+
  geom_histogram(color="white",fill = "purple", alpha=0.7)+
  labs(x = "Age in years", y= "Number of Participants", title ="Age at menarche (earliest age given)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"))

#Male PS & GS plotted together
ggplot (Clindata) +
  labs(color="Variable name",x="Age in years",y= "Number of Participants", title="Earliest Questionnaire Response of Tanner Stage 2", fill="Puberty Marker")+
  geom_histogram(aes(x=Clindata$MPY, fill= "Male Pubic Hair"), position = "dodge",alpha = 0.6) + 
  geom_histogram(aes(x=Clindata$MGY, fill= "Male Genitalia"), position = "dodge", alpha = 0.6)+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"), legend.text=element_text(size=14))

#Breast stage & female pubic hair stage plotted together
ggplot (Clindata) +
  labs(color="Variable name",x="Age in years",y="Number of Participants", title="Earliest Questionnaire Response of Tanner Stage 2", fill="Puberty Marker")+
  geom_histogram(aes(x=Clindata$BY, fill= "Breast Development"), position = "dodge", alpha = 0.6) + 
  geom_histogram(aes(x=Clindata$FPY, fill= "Female Pubic Hair"), position = "dodge", alpha = 0.6) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"), legend.text=element_text(size=14))

#Save a csv file of the clinical Clindata with the added columns including ATS2
write.csv(Clindata,file="AgeatTS2_years.csv")

#Load transcriptomics file
load("~/Documents/Code/EXP.Rdata")

#New object with transcriptomic Clindata transposed
Expdata_t <- (t(data))
rm(data)

#remove row 642 and 728 due to -9; Patient IDs #-----A and ----A
which((rownames(Expdata_t)=="-----A"))
which((rownames(Expdata_t)=="----A"))
Expdata_t <- Expdata_t[-c(642,728),]

#Object with list of participant IDs with transcriptomics
TransPIDs <- rownames(Expdata_t)

#Positions of participant ID within Clindata for those participants with transcriptome data
Pos<- match(TransPIDs,rownames(Clindata))

#New object: just clinical data for those with transcriptomics 
ClindataT<- Clindata[Pos,]

#Subsets of males & females transcriptomics (females n=497, males n=450)
ClindataTF <- ClindataT[ClindataT$kz021==2,]
ClindataTM <- ClindataT[ClindataT$kz021==1,]

#Save csv files of the clinical data for those with transcriptomics 
write.csv(ClindataT,file="AgeatTS2_yearsT.csv")
write.csv(ClindataTF,file="AgeatTS2_years_TF.csv")
write.csv(ClindataTM,file="AgeatTS2_years_TM.csv")

