setwd("~/Documents/Code")
library(ggplot2)

#Calculate time from methylome sample to puberty----
#Load clinical data
TimeMtopuberty <- read.csv("~/Documents/Code/AgeatTS2_years.csv")

#Load methylome sample data
load("~/Documents/Code/samplesheet/samplesheet.Robj")

#Remove duplicate rows as per ALSPAC methylome user guide
Dups<-which(samplesheet$duplicate.rm=="Remove")
samplesheet<-samplesheet[-Dups,]

#Change row name to participant ID plus rank (combine columns 2 & 3, no spaces), as for Clin Data file
rownames(samplesheet)<-paste(samplesheet[,2],samplesheet[,3],sep="")

#Subset of clinical data to include only those with methylome samples
AgeatTS2_yearsM <- TimeMtopuberty[(match(rownames(samplesheet),TimeMtopuberty$X)),]

write.csv(AgeatTS2_yearsM,file="AgeatTS2_yearsM.csv")

TimeMtopuberty <- AgeatTS2_yearsM

#Create new variable for male & female age at pubic hair stage 2, if male use MPY, otherwise use FPY
TimeMtopuberty$MFPY <- ifelse(TimeMtopuberty$kz021==1,TimeMtopuberty$MPY, TimeMtopuberty$FPY)

#f7003c is age in months at Focus@7 visit (when sample collected for methylome)
summary(TimeMtopuberty$f7003c)

#Checking values for Focus@7 visit - start by duplicating column at end
TimeMtopuberty$f7003cdup <- TimeMtopuberty$f7003c

#Convert to age in months at Focus@7 visit into years
TimeMtopuberty$f7003cyrs <-TimeMtopuberty$f7003c/12
summary(TimeMtopuberty$f7003cyrs)

#MeP2MF=time from methylome sample to age at earliest questionnaire response when P2
TimeMtopuberty$MeP2MF <-TimeMtopuberty$MFPY-(TimeMtopuberty$f7003cyrs)
summary(TimeMtopuberty$MeP2MF)

ggplot(TimeMtopuberty,aes(x=MeP2MF))+
  geom_histogram(color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from methylome sample to pubic hair stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#MeP2M=time from methylome sample to age at earliest questionnaire response when P2 in males
TimeMtopuberty$MeP2M <-TimeMtopuberty$MPY-(TimeMtopuberty$f7003cyrs)
summary(TimeMtopuberty$MeP2M)

ggplot(TimeMtopuberty,aes(x=MeP2M))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from methylome sample to pubic hair stage 2 - males")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#MeP2F=time from methylome sample to age at earliest questionnaire response when P2 in females
TimeMtopuberty$MeP2F <-TimeMtopuberty$FPY-(TimeMtopuberty$f7003cyrs)
summary(TimeMtopuberty$MeP2F)

ggplot(TimeMtopuberty,aes(x=MeP2F))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from methylome sample to pubic hair stage 2 - females")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

TimeMtopuberty$kz021 <- as.factor(TimeMtopuberty$kz021)

ggplot(TimeMtopuberty,aes(x=MeP2MF,group=kz021,fill=kz021, after_stat(density)))+
  geom_histogram(binwidth=0.5, position = "identity", color="grey", alpha=0.5)+
  scale_fill_discrete(labels=c('Male', 'Female'))+
  labs(x = "Time in years", y= "Density", title ="Time from methylome sample to pubic hair stage 2", fill="Sex")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"), legend.text=element_text(size=14))+
  geom_density(alpha=0.25)

#MeB2=time from methylome sample to age at earliest questionnaire response when B2
TimeMtopuberty$MeB2 <-TimeMtopuberty$BY-(TimeMtopuberty$f7003cyrs)
summary(TimeMtopuberty$MeB2)

ggplot(TimeMtopuberty,aes(x=MeB2))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from methylome sample to breast stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#MeMen=time from methylome sample to age at earliest questionnaire response for age at menarche
TimeMtopuberty$MeMen <-TimeMtopuberty$MENY-(TimeMtopuberty$f7003cyrs)
summary(TimeMtopuberty$MeMen)

ggplot(TimeMtopuberty,aes(x=MeMen))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from methylome sample to menarche")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#Calculate quartiles----
#P2MF
#Remove rows with NA for MeP2MF
P2MFQ <- TimeMtopuberty[!is.na(TimeMtopuberty$MeP2MF),]

#Add quartile column & fill as NA
P2MFQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2 <- quantile(P2MFQ$MeP2MF)

#Populate DistancetoP column with 1st-4th quartile
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$MeP2MF<Quartiles.DistancetoP2[2],"1st Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$MeP2MF>=Quartiles.DistancetoP2[2] & P2MFQ$MeP2MF<Quartiles.DistancetoP2[3], "2nd Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$MeP2MF>=Quartiles.DistancetoP2[3] & P2MFQ$MeP2MF<Quartiles.DistancetoP2[4], "3rd Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$MeP2MF>=Quartiles.DistancetoP2[4], "4th Quartile", P2MFQ$DtoPQuartile)

summary(P2MFQ$f7003cyrs)

#Save a csv of clinical data with puberty quartiles for those with methylome samples
write.csv(P2MFQ,file="P2MFquartilesM.csv")

#P2M
#Remove rows with NA for MeP2M
P2MQ <- TimeMtopuberty[!is.na(TimeMtopuberty$MeP2M),]

#Add quartile column & fill as NA
P2MQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2M <- quantile(P2MQ$MeP2M)

#Populate DistancetoP column with 1st-4th quartile
P2MQ$DtoPQuartile <- ifelse(P2MQ$MeP2M<Quartiles.DistancetoP2M[2],"1st Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$MeP2M>=Quartiles.DistancetoP2M[2] & P2MQ$MeP2M<Quartiles.DistancetoP2[3], "2nd Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$MeP2M>=Quartiles.DistancetoP2M[3] & P2MQ$MeP2M<Quartiles.DistancetoP2[4], "3rd Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$MeP2M>=Quartiles.DistancetoP2M[4], "4th Quartile", P2MQ$DtoPQuartile)

summary(P2MQ$f7003cyrs)

#Save a csv of clinical data with puberty quartiles for those with methylome samples
write.csv(P2MQ,file="P2MquartilesM.csv")

#P2F
#Remove rows with NA for MeP2F
P2FQ <- TimeMtopuberty[!is.na(TimeMtopuberty$MeP2F),]

#Add quartile column & fill as NA
P2FQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2F <- quantile(P2FQ$MeP2F)

#Populate DistancetoP column with 1st-4th quartile
P2FQ$DtoPQuartile <- ifelse(P2FQ$MeP2F<Quartiles.DistancetoP2F[2],"1st Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$MeP2F>=Quartiles.DistancetoP2F[2] & P2FQ$MeP2F<Quartiles.DistancetoP2F[3], "2nd Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$MeP2F>=Quartiles.DistancetoP2F[3] & P2FQ$MeP2F<Quartiles.DistancetoP2F[4], "3rd Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$MeP2F>=Quartiles.DistancetoP2F[4], "4th Quartile", P2FQ$DtoPQuartile)

summary(P2FQ$f7003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(P2FQ,file="P2FquartilesM.csv")

#B2
#Remove rows with NA for MeB2
B2Q <- TimeMtopuberty[!is.na(TimeMtopuberty$MeB2),]

#Add quartile column & fill as NA
B2Q$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoB2 <- quantile(B2Q$MeB2)

#Populate DistancetoP column with 1st-4th quartile
B2Q$DtoPQuartile <- ifelse(B2Q$MeB2<Quartiles.DistancetoB2[2],"1st Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$MeB2>=Quartiles.DistancetoB2[2] & B2Q$MeB2<Quartiles.DistancetoB2[3], "2nd Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$MeB2>=Quartiles.DistancetoB2[3] & B2Q$MeB2<Quartiles.DistancetoB2[4], "3rd Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$MeB2>=Quartiles.DistancetoB2[4], "4th Quartile", B2Q$DtoPQuartile)

summary(B2Q$f7003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(B2Q,file="B2quartilesM.csv")

#Menarche
#Remove rows with NA for MeMen
MenQ <- TimeMtopuberty[!is.na(TimeMtopuberty$MeMen),]

#Add quartile column & fill as NA
MenQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoMen <- quantile(MenQ$MeMen)

#Populate DistancetoP column with 1st-4th quartile
MenQ$DtoPQuartile <- ifelse(MenQ$MeMen<Quartiles.DistancetoMen[2],"1st Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$MeMen>=Quartiles.DistancetoMen[2] & MenQ$MeMen<Quartiles.DistancetoMen[3], "2nd Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$MeMen>=Quartiles.DistancetoMen[3] & MenQ$MeMen<Quartiles.DistancetoMen[4], "3rd Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$MeMen>=Quartiles.DistancetoMen[4], "4th Quartile", MenQ$DtoPQuartile)

summary(MenQ$f7003cyrs)

#Save a csv of clinical data with puberty quartiles for those with methylome samples
write.csv(MenQ,file="MenquartilesM.csv")

