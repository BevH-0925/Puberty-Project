setwd("~/Documents/Code")
library(ggplot2)
library(mixOmics)
library(rgl)

#Calculate time from transcriptomic sample to puberty----
#Load clinical data for transcriptomic cohort
TimeTtopuberty <- read.csv("~/Documents/Code/AgeatTS2_yearsT.csv")

#Create new variable for male & female age at pubic hair stage 2, if male use MPY, otherwise use FPY
TimeTtopuberty$MFPY <- ifelse(TimeTtopuberty$kz021==1,TimeTtopuberty$MPY, TimeTtopuberty$FPY)

#f9003c is age in months at Focus@9 visit (when sample collected for transcriptomics)
summary(TimeTtopuberty$f9003c)

#Checking values for Focus@9 visit - start by duplicating column at end
TimeTtopuberty$f9003cdup <- TimeTtopuberty$f9003c

#Convert to age in months at Focus@9 visit into years
TimeTtopuberty$f9003cyrs <-TimeTtopuberty$f9003c/12
summary(TimeTtopuberty$f9003cyrs)

#TrP2MF=time from transcriptomics sample to age at earliest questionnaire response when P2
TimeTtopuberty$TrP2MF <-TimeTtopuberty$MFPY-(TimeTtopuberty$f9003cyrs)
summary(TimeTtopuberty$TrP2MF)

ggplot(TimeTtopuberty,aes(x=TrP2MF))+
  geom_histogram(color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from transcriptomics sample to pubic hair stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#TrP2M=time from transcriptomics sample to age at earliest questionnaire response when P2 in males
TimeTtopuberty$TrP2M <-TimeTtopuberty$MPY-(TimeTtopuberty$f9003cyrs)
summary(TimeTtopuberty$TrP2M)

ggplot(TimeTtopuberty,aes(x=TrP2M))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from transcriptomics sample to pubic hair stage 2 - males")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#TrP2F=time from transcriptomics sample to age at earliest questionnaire response when P2 in females
TimeTtopuberty$TrP2F <-TimeTtopuberty$FPY-(TimeTtopuberty$f9003cyrs)
summary(TimeTtopuberty$TrP2F)

ggplot(TimeTtopuberty,aes(x=TrP2F))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from transcriptomics sample to pubic hair stage 2 - females")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

TimeTtopuberty$kz021 <- as.factor(TimeTtopuberty$kz021)

ggplot(TimeTtopuberty,aes(x=TrP2MF,group=kz021,fill=kz021, after_stat(density)))+
  geom_histogram(binwidth=0.5, position = "identity", color="grey", alpha=0.5)+
  scale_fill_discrete(labels=c('Male', 'Female'))+
  labs(x = "Time in years", y= "Density", title ="Time from transcriptomics sample to pubic hair stage 2", fill="Sex")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=16,face="bold"), legend.text=element_text(size=14))+
  geom_density(alpha=0.25)

#TrB2=time from transcriptomics sample to age at earliest questionnaire response when B2
TimeTtopuberty$TrB2 <-TimeTtopuberty$BY-(TimeTtopuberty$f9003cyrs)
summary(TimeTtopuberty$TrB2)

ggplot(TimeTtopuberty,aes(x=TrB2))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from transcriptomics sample to breast stage 2")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#TrMen=time from transcriptomics sample to age at earliest questionnaire response for age at menarche
TimeTtopuberty$TrMen <-TimeTtopuberty$MENY-(TimeTtopuberty$f9003cyrs)
summary(TimeTtopuberty$TrMen)

ggplot(TimeTtopuberty,aes(x=TrMen))+
  geom_histogram(binwidth=0.5, color="white",fill = "blue")+
  labs(x = "Time in years", y= "Number of Participants", title ="Time from transcriptomics sample to menarche")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"), title = element_text(size=14,face="bold"))

#Calculate quartiles----
#P2MF
#Remove rows with NA for TrP2MF
P2MFQ <- TimeTtopuberty[!is.na(TimeTtopuberty$TrP2MF),]

#Add quartile column & fill as NA
P2MFQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2 <- quantile(P2MFQ$TrP2MF)

#Populate DistancetoP column with 1st-4th quartile
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$TrP2MF<Quartiles.DistancetoP2[2],"1st Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$TrP2MF>=Quartiles.DistancetoP2[2] & P2MFQ$TrP2MF<Quartiles.DistancetoP2[3], "2nd Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$TrP2MF>=Quartiles.DistancetoP2[3] & P2MFQ$TrP2MF<Quartiles.DistancetoP2[4], "3rd Quartile", P2MFQ$DtoPQuartile)
P2MFQ$DtoPQuartile <- ifelse(P2MFQ$TrP2MF>=Quartiles.DistancetoP2[4], "4th Quartile", P2MFQ$DtoPQuartile)

summary(P2MFQ$f9003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(P2MFQ,file="P2MFquartilesT.csv")

#P2M
#Remove rows with NA for TrP2M
P2MQ <- TimeTtopuberty[!is.na(TimeTtopuberty$TrP2M),]

#Add quartile column & fill as NA
P2MQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2M <- quantile(P2MQ$TrP2M)

#Populate DistancetoP column with 1st-4th quartile
P2MQ$DtoPQuartile <- ifelse(P2MQ$TrP2M<Quartiles.DistancetoP2M[2],"1st Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$TrP2M>=Quartiles.DistancetoP2M[2] & P2MQ$TrP2M<Quartiles.DistancetoP2[3], "2nd Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$TrP2M>=Quartiles.DistancetoP2M[3] & P2MQ$TrP2M<Quartiles.DistancetoP2[4], "3rd Quartile", P2MQ$DtoPQuartile)
P2MQ$DtoPQuartile <- ifelse(P2MQ$TrP2M>=Quartiles.DistancetoP2M[4], "4th Quartile", P2MQ$DtoPQuartile)

summary(P2MQ$f9003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(P2MQ,file="P2MquartilesT.csv")

#P2F
#Remove rows with NA for TrP2F
P2FQ <- TimeTtopuberty[!is.na(TimeTtopuberty$TrP2F),]

#Add quartile column & fill as NA
P2FQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoP2F <- quantile(P2FQ$TrP2F)

#Populate DistancetoP column with 1st-4th quartile
P2FQ$DtoPQuartile <- ifelse(P2FQ$TrP2F<Quartiles.DistancetoP2F[2],"1st Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$TrP2F>=Quartiles.DistancetoP2F[2] & P2FQ$TrP2F<Quartiles.DistancetoP2F[3], "2nd Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$TrP2F>=Quartiles.DistancetoP2F[3] & P2FQ$TrP2F<Quartiles.DistancetoP2F[4], "3rd Quartile", P2FQ$DtoPQuartile)
P2FQ$DtoPQuartile <- ifelse(P2FQ$TrP2F>=Quartiles.DistancetoP2F[4], "4th Quartile", P2FQ$DtoPQuartile)

summary(P2FQ$f9003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(P2FQ,file="P2FquartilesT.csv")

#B2
#Remove rows with NA for TrB2
B2Q <- TimeTtopuberty[!is.na(TimeTtopuberty$TrB2),]

#Add quartile column & fill as NA
B2Q$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoB2 <- quantile(B2Q$TrB2)

#Populate DistancetoP column with 1st-4th quartile
B2Q$DtoPQuartile <- ifelse(B2Q$TrB2<Quartiles.DistancetoB2[2],"1st Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$TrB2>=Quartiles.DistancetoB2[2] & B2Q$TrB2<Quartiles.DistancetoB2[3], "2nd Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$TrB2>=Quartiles.DistancetoB2[3] & B2Q$TrB2<Quartiles.DistancetoB2[4], "3rd Quartile", B2Q$DtoPQuartile)
B2Q$DtoPQuartile <- ifelse(B2Q$TrB2>=Quartiles.DistancetoB2[4], "4th Quartile", B2Q$DtoPQuartile)

summary(B2Q$f9003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(B2Q,file="B2quartilesT.csv")

#Menarche
#Remove rows with NA for TrMen
MenQ <- TimeTtopuberty[!is.na(TimeTtopuberty$TrMen),]

#Add quartile column & fill as NA
MenQ$DtoPQuartile <-NA 

#Calculate quartiles of distance to puberty & put in new object
Quartiles.DistancetoMen <- quantile(MenQ$TrMen)

#Populate DistancetoP column with 1st-4th quartile
MenQ$DtoPQuartile <- ifelse(MenQ$TrMen<Quartiles.DistancetoMen[2],"1st Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$TrMen>=Quartiles.DistancetoMen[2] & MenQ$TrMen<Quartiles.DistancetoMen[3], "2nd Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$TrMen>=Quartiles.DistancetoMen[3] & MenQ$TrMen<Quartiles.DistancetoMen[4], "3rd Quartile", MenQ$DtoPQuartile)
MenQ$DtoPQuartile <- ifelse(MenQ$TrMen>=Quartiles.DistancetoMen[4], "4th Quartile", MenQ$DtoPQuartile)

summary(MenQ$f9003cyrs)

#Save a csv of clinical data with puberty quartiles for those with transcriptomics
write.csv(MenQ,file="MenquartilesT.csv")

#Subset of transcriptomics data----
#Load transcriptomics file
load("~/Documents/Code/EXP.Rdata")

#New data frame with data transposed
Expdata_t <- t(data)
#remove row 642 and 728 due to -9; Patient IDs #-----A and ----A
Expdata_t <- Expdata_t[-c(642,728),]

rm(data)

#Subset transcriptomic data to only include participants with distance to puberty quartile.
Expdata_P2MF <- Expdata_t[match(P2MFQ$X,rownames(Expdata_t)),]

#Subset transcriptomic data to only include female participants with distance to puberty quartile.
Expdata_P2M <- Expdata_t[match(P2MQ$X,rownames(Expdata_t)),]

#Subset transcriptomic data to only include male participants with distance to puberty quartile.
Expdata_P2F <- Expdata_t[match(P2FQ$X,rownames(Expdata_t)),]

#Subset transcriptomic data to only include female participants with distance to puberty quartile.
Expdata_B2 <- Expdata_t[match(B2Q$X,rownames(Expdata_t)),]

#Subset transcriptomic data to only include male participants with distance to puberty quartile.
Expdata_Men <- Expdata_t[match(MenQ$X,rownames(Expdata_t)),]

#New object with just the participant ID & distance to P2 quartile
PubertyGpsP2MF <- P2MFQ[,c("X","DtoPQuartile")]

#New object with just the participant ID & distance to P2 quartile for males
PubertyGpsP2M <- P2MQ[,c("X","DtoPQuartile")]

#New object with just the participant ID & distance to P2 quartile for females
PubertyGpsP2F <- P2FQ[,c("X","DtoPQuartile")]

#New object with just the participant ID & distance to B2 quartile 
PubertyGpsB2 <- B2Q[,c("X","DtoPQuartile")]

#New object with just the participant ID & distance to Menarche quartile 
PubertyGpsMen <- MenQ[,c("X","DtoPQuartile")]

