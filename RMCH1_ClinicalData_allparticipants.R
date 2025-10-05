setwd("~/Documents/Code/RMCH Clinical Data")
library(tidyverse)
library(psych) 

#Load Data----
#See Excels on work laptop for notes on data cleanup
clindata <- read.csv("~/Documents/Code/RMCH Clinical Data/Data/Clindata_all.csv")

clindata$Sex<- as.factor(clindata$Sex)
clindata$Ethnicity<- as.factor(clindata$Ethnicity)
clindata$BreastStage <- as.factor(clindata$BreastStage)
clindata$GenitalStage<- as.factor(clindata$GenitalStage)
clindata$PubicHairStage<- as.factor(clindata$PubicHairStage)
clindata$AxillaryHairStage<- as.factor(clindata$AxillaryHairStage)
clindata$PituitaryMR<- as.factor(clindata$PituitaryMR)
clindata$BoneAgeXray<- as.factor(clindata$BoneAgeXray)
clindata$PelvicUltrasound<- as.factor(clindata$PelvicUltrasound)
clindata$GnRHInterpretation<- as.factor(clindata$GnRHInterpretation)
clindata$Group<- as.factor(clindata$Group)

clindata$HeightSDS <- as.numeric(clindata$HeightSDS)
clindata$WeightSDS <- as.numeric(clindata$WeightSDS)
clindata$LH30min <- as.numeric(clindata$LH30min)
clindata$LH60min <- as.numeric(clindata$LH60min)
clindata$FSH30min <- as.numeric(clindata$FSH30min)
clindata$FSH60min <- as.numeric(clindata$FSH60min)
clindata$PeakLH <- as.numeric(clindata$PeakLH)
clindata$PeakFSH <- as.numeric(clindata$PeakFSH)
clindata$PeakLHFSHratio <- as.numeric(clindata$PeakLHFSHratio)
clindata$TestosteroneConc <- as.numeric(clindata$TestosteroneConc)
clindata$OestradiolConc <- as.numeric(clindata$OestradiolConc)

clindata_males <- clindata[clindata$Sex=="Male",]
clindata_females <- clindata[clindata$Sex=="Female",]

clindata_prepub <- clindata[clindata$Group=="Prepubertal",]
clindata_pub <- clindata[clindata$Group=="Pubertal",]

clindata_prepub_males <- clindata_prepub[clindata_prepub$Sex=="Male",]
clindata_prepub_females <- clindata_prepub[clindata_prepub$Sex=="Female",]
clindata_pub_males <- clindata_pub[clindata_pub$Sex=="Male",]
clindata_pub_females <- clindata_pub[clindata_pub$Sex=="Female",]

#Summary tables----
Descriptive_all <- describe(clindata[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_all <- sapply((clindata[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
#CountNA only use for numeric variables
CountNA_all <- sapply(clindata, function(x) sum(is.na(x)))

write.csv(Descriptive_all,file="Descriptive_all.csv")
write.csv(Quartiles_all,file="Quartiles_all.csv")
write.csv(CountNA_all,file="CountNA_all.csv")

Descriptive_males <- describe(clindata_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_males <- sapply((clindata_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_males <- sapply(clindata_males, function(x) sum(is.na(x)))

write.csv(Descriptive_males,file="Descriptive_males.csv")
write.csv(Quartiles_males,file="Quartiles_males.csv")
write.csv(CountNA_males,file="CountNA_males.csv")

Descriptive_females<- describe(clindata_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_females <- sapply((clindata_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_females <- sapply(clindata_females, function(x) sum(is.na(x)))

write.csv(Descriptive_females,file="Descriptive_females.csv")
write.csv(Quartiles_females,file="Quartiles_females.csv")
write.csv(CountNA_females,file="CountNA_females.csv")

Descriptive_prepub <- describe(clindata_prepub[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_prepub <- sapply((clindata_prepub[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_prepub <- sapply(clindata_prepub, function(x) sum(is.na(x)))

write.csv(Descriptive_prepub,file="Descriptive_prepub.csv")
write.csv(Quartiles_prepub,file="Quartiles_prepub.csv")
write.csv(CountNA_prepub,file="CountNA_prepub.csv")

Descriptive_pub <- describe(clindata_pub[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_pub <- sapply((clindata_pub[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_pub <- sapply(clindata_pub, function(x) sum(is.na(x)))

write.csv(Descriptive_pub,file="Descriptive_pub.csv")
write.csv(Quartiles_pub,file="Quartiles_pub.csv")
write.csv(CountNA_pub,file="CountNA_pub.csv")

Descriptive_prepub_males <- describe(clindata_prepub_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_prepub_males <- sapply((clindata_prepub_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_prepub_males <- sapply(clindata_prepub_males, function(x) sum(is.na(x)))

write.csv(Descriptive_prepub_males,file="Descriptive_prepub_males.csv")
write.csv(Quartiles_prepub_males,file="Quartiles_prepub_males.csv")
write.csv(CountNA_prepub_males,file="CountNA_prepub_males.csv")

Descriptive_pub_males <- describe(clindata_pub_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_pub_males <- sapply((clindata_pub_males[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_pub_males <- sapply(clindata_pub_males, function(x) sum(is.na(x)))

write.csv(Descriptive_pub_males,file="Descriptive_pub_males.csv")
write.csv(Quartiles_pub_males,file="Quartiles_pub_males.csv")
write.csv(CountNA_pub_males,file="CountNA_pub_males.csv")

Descriptive_prepub_females <- describe(clindata_prepub_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_prepub_females <- sapply((clindata_prepub_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_prepub_females <- sapply(clindata_prepub_females, function(x) sum(is.na(x)))

write.csv(Descriptive_prepub_females,file="Descriptive_prepub_females.csv")
write.csv(Quartiles_prepub_females,file="Quartiles_prepub_females.csv")
write.csv(CountNA_prepub_females,file="CountNA_prepub_females.csv")

Descriptive_pub_females <- describe(clindata_pub_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')])
Quartiles_pub_females <- sapply((clindata_pub_females[,c('HeightSDS','WeightSDS','TesticularVolume', 'Age_yrs', 'LH0min','FSH0min','LH30min','FSH30min','LH60min','FSH60min','PeakLH','PeakFSH','BasalLHFSHratio','PeakLHFSHratio', 'TestosteroneConc','OestradiolConc')]), function(x) quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
CountNA_pub_females <- sapply(clindata_pub_females, function(x) sum(is.na(x)))

write.csv(Descriptive_pub_females,file="Descriptive_pub_females.csv")
write.csv(Quartiles_pub_females,file="Quartiles_pub_females.csv")
write.csv(CountNA_pub_females,file="CountNA_pub_females.csv")

#Plots----
ggplot(clindata, aes(x=Age_yrs, y=PeakLH, color=Group))+
geom_point(aes(shape=Sex))+
labs(x = "Age at GnRH test (years)", y = "Peak LH Concentration (IU/L)")+
scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18))+
scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
geom_abline(slope=0,intercept=6,lwd=1,lty=2,col="grey")+
theme_classic()
ggsave("AgeLHscatter.png",dpi=320)

clindataLH<- clindata[,c(16,18,20,29)]
clindataLHpivotlong <- clindataLH|> pivot_longer(
  cols = starts_with("LH"), 
  names_to = "LH", 
  values_to = "Concentration"
)
 
clindataLHpivotlong[clindataLHpivotlong == 'LH0min'] <- 'Basal LH'
clindataLHpivotlong[clindataLHpivotlong == 'LH30min'] <- '30 min LH'
clindataLHpivotlong[clindataLHpivotlong == 'LH60min'] <- '60 min LH'

ggplot(clindataLHpivotlong, aes(x = Group, y = Concentration, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "LH Concentration (IU/L)", fill="Puberty Classification")+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=6,lwd=1,lty=2,col="grey")+
  theme_classic()+
  guides(fill = FALSE)+
  facet_wrap(~factor(LH, levels=c('Basal LH', '30 min LH', '60 min LH')))
ggsave("GroupLHallBoxplot.png", dpi=320)   

clindataFSH<- clindata[,c(17,19,21,29)]
clindataFSHpivotlong <- clindataFSH|> pivot_longer(
  cols = starts_with("FSH"), 
  names_to = "FSH", 
  values_to = "Concentration"
)

clindataFSHpivotlong[clindataFSHpivotlong == 'FSH0min'] <- 'Basal FSH'
clindataFSHpivotlong[clindataFSHpivotlong == 'FSH30min'] <- '30 min FSH'
clindataFSHpivotlong[clindataFSHpivotlong == 'FSH60min'] <- '60 min FSH'

ggplot(clindataFSHpivotlong, aes(x = Group, y = Concentration, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "FSH Concentration (IU/L)", fill="Puberty Classification")+
  scale_y_continuous(breaks = c(0,10,20,30,40))+
  theme_classic()+
  guides(fill = FALSE)+
  facet_wrap(~factor(FSH, levels=c('Basal FSH', '30 min FSH', '60 min FSH')))
ggsave("GroupFSHallBoxplot.png", dpi=320)   

ggplot(clindata, aes(x = Group, y = Age_yrs, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "Age (years)", fill="Puberty Classification")+
  theme_classic()+
  guides(fill = FALSE)+
  facet_wrap(~Sex)
ggsave("GroupAgeSexBoxplot.png", dpi=320)   

ggplot(clindata, aes(x = Group, y = PeakLHFSHratio, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "Peak LH/FSH ratio", fill="Puberty Classification")+
  theme_classic()+
  guides(fill = FALSE)
ggsave("GroupPeakLHFSHBoxplot.png", dpi=320)  

ggplot(clindata_females, aes(x = Group, y = OestradiolConc, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "Oestradiol Concentration in females (pmol/L)", fill="Puberty Classification")+
  theme_classic()+
  guides(fill = FALSE)
ggsave("GroupOestFBoxplot.png", dpi=320)  

ggplot(clindata_males, aes(x = Group, y = TestosteroneConc, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "Testosterone Concentration in males (nmol/L)", fill="Puberty Classification")+
  theme_classic()+
  guides(fill = FALSE)
ggsave("GroupTestoMBoxplot.png", dpi=320)  

ggplot(clindata, aes(x = Group, y = BasalLHFSHratio, fill=Group)) +
  geom_boxplot()+
  labs(x = "Puberty Classification", y = "Basal LH/FSH ratio", fill="Puberty Classification")+
  theme_classic()+
  guides(fill = FALSE)
ggsave("GroupBasalLHFSHBoxplot.png", dpi=320) 

#Stats - male and female together----
#Fishers 
#Fisher's Exact test - Sex & Puberty
data <- clindata[,c("Sex", "Group")]
data <- na.omit(data)
data <- table(data$Sex, data$Group)
data
fisher.test(data)

#Fisher's Exact test - Ethnicity & Puberty
data <- clindata[,c("Ethnicity", "Group")]
data <- na.omit(data)
data <- table(data$Ethnicity, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Breast Stage & Puberty
data <- clindata[,c("BreastStage", "Group")]
data <- na.omit(data)
data <- table(data$BreastStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Genital Stage & Puberty
data <- clindata[,c("GenitalStage", "Group")]
data <- na.omit(data)
data <- table(data$GenitalStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pubic Hair Stage & Puberty
data <- clindata[,c("PubicHairStage", "Group")]
data <- na.omit(data)
data <- table(data$PubicHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Axillary Hair Stage & Puberty
data <- clindata[,c("AxillaryHairStage", "Group")]
data <- na.omit(data)
data <- table(data$AxillaryHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pituitary MR & Puberty
data <- clindata[,c("PituitaryMR", "Group")]
data <- na.omit(data)
data <- table(data$PituitaryMR, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Bone Age X-ray & Puberty
data <- clindata[,c("BoneAgeXray", "Group")]
data <- na.omit(data)
data <- table(data$BoneAgeXray, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pelvic US & Puberty
data <- clindata[,c("PelvicUltrasound", "Group")]
data <- na.omit(data)
data <- table(data$PelvicUltrasound, data$Group)
data
fisher.test(data)

#Measured Variables
#NAs automatically omitted from test
#Height SDS
Prepub_HeightSDS <- clindata_prepub$HeightSDS
Pub_HeightSDS <- clindata_pub$HeightSDS
#Perform the Mann Whitney U test (uses Wilcox function in R)
wilcox.test(Prepub_HeightSDS, Pub_HeightSDS, conf.int = TRUE, conf.level = 0.95)

#Weight SDS
Prepub_WeightSDS <- clindata_prepub$WeightSDS
Pub_WeightSDS <- clindata_pub$WeightSDS
wilcox.test(Prepub_WeightSDS, Pub_WeightSDS, conf.int = TRUE, conf.level = 0.95)

#Testicular Volume
Prepub_TV <- clindata_prepub$TesticularVolume
Pub_TV<- clindata_pub$TesticularVolume
wilcox.test(Prepub_TV, Pub_TV, conf.int = TRUE, conf.level = 0.95)

#Age 
Prepub_Age <- clindata_prepub$Age_yrs
Pub_Age <- clindata_pub$Age_yrs
wilcox.test(Prepub_Age, Pub_Age, conf.int = TRUE, conf.level = 0.95)

#LH0min
Prepub_LH0min <- clindata_prepub$LH0min
Pub_LH0min<- clindata_pub$LH0min
wilcox.test(Prepub_LH0min, Pub_LH0min, conf.int = TRUE, conf.level = 0.95)

#LH30min
Prepub_LH30min <- clindata_prepub$LH30min
Pub_LH30min<- clindata_pub$LH30min
wilcox.test(Prepub_LH30min, Pub_LH30min,conf.int = TRUE, conf.level = 0.95)

#LH60min
Prepub_LH60min <- clindata_prepub$LH60min
Pub_LH60min<- clindata_pub$LH60min
wilcox.test(Prepub_LH60min, Pub_LH60min,conf.int = TRUE, conf.level = 0.95)

#FSH0min
Prepub_FSH0min <- clindata_prepub$FSH0min
Pub_FSH0min<- clindata_pub$FSH0min
wilcox.test(Prepub_FSH0min, Pub_FSH0min,conf.int = TRUE, conf.level = 0.95)

#FSH30min
Prepub_FSH30min <- clindata_prepub$FSH30min
Pub_FSH30min<- clindata_pub$FSH30min
wilcox.test(Prepub_FSH30min, Pub_FSH30min,conf.int = TRUE, conf.level = 0.95)

#FSH60min
Prepub_FSH60min <- clindata_prepub$FSH60min
Pub_FSH60min<- clindata_pub$FSH60min
wilcox.test(Prepub_FSH60min, Pub_FSH60min,conf.int = TRUE, conf.level = 0.95)

#Peak LH
Prepub_PeakLH <- clindata_prepub$PeakLH
Pub_PeakLH <- clindata_pub$PeakLH
wilcox.test(Prepub_PeakLH, Pub_PeakLH,conf.int = TRUE, conf.level = 0.95)

#Peak FSH
Prepub_PeakFSH <- clindata_prepub$PeakFSH
Pub_PeakFSH <- clindata_pub$PeakFSH
wilcox.test(Prepub_PeakFSH, Pub_PeakFSH,conf.int = TRUE, conf.level = 0.95)

#Basal LH/FSH Ratio
Prepub_BasalLHFSHratio <- clindata_prepub$BasalLHFSHratio
Pub_BasalLHFSHratio <- clindata_pub$BasalLHFSHratio
wilcox.test(Prepub_BasalLHFSHratio, Pub_BasalLHFSHratio,conf.int = TRUE, conf.level = 0.95)

#Peak LH/FSH Ratio
Prepub_PeakLHFSHratio <- clindata_prepub$PeakLHFSHratio
Pub_PeakLHFSHratio <- clindata_pub$PeakLHFSHratio
wilcox.test(Prepub_PeakLHFSHratio, Pub_PeakLHFSHratio,conf.int = TRUE, conf.level = 0.95)

#Testosterone Conc
Prepub_Testo <- clindata_prepub$TestosteroneConc
Pub_Testo <- clindata_pub$TestosteroneConc
wilcox.test(Prepub_Testo, Pub_Testo,conf.int = TRUE, conf.level = 0.95)

#Oestradiol Conc
Prepub_Oest <- clindata_prepub$OestradiolConc
Pub_Oest <- clindata_pub$OestradiolConc
wilcox.test(Prepub_Oest, Pub_Oest,conf.int = TRUE, conf.level = 0.95)

#Stats - males only----
#Fishers 
#Fisher's Exact test - Sex & Puberty
data <- clindata_males[,c("Sex", "Group")]
data <- na.omit(data)
data <- table(data$Sex, data$Group)
data
fisher.test(data)

#Fisher's Exact test - Ethnicity & Puberty
data <- clindata_males[,c("Ethnicity", "Group")]
data <- na.omit(data)
data <- table(data$Ethnicity, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Genital Stage & Puberty
data <- clindata_males[,c("GenitalStage", "Group")]
data <- na.omit(data)
data <- table(data$GenitalStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pubic Hair Stage & Puberty
data <- clindata_males[,c("PubicHairStage", "Group")]
data <- na.omit(data)
data <- table(data$PubicHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Axillary Hair Stage & Puberty
data <- clindata_males[,c("AxillaryHairStage", "Group")]
data <- na.omit(data)
data <- table(data$AxillaryHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pituitary MR & Puberty
data <- clindata_males[,c("PituitaryMR", "Group")]
data <- na.omit(data)
data <- table(data$PituitaryMR, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Bone Age X-ray & Puberty
data <- clindata_males[,c("BoneAgeXray", "Group")]
data <- na.omit(data)
data <- table(data$BoneAgeXray, data$Group)
data
fisher.test(data)

#Measured Variables
#NAs automatically omitted from test
#Height SDS
Prepub_HeightSDS_males <- clindata_prepub_males$HeightSDS
Pub_HeightSDS_males <- clindata_pub_males$HeightSDS
#Perform the Mann Whitney U test (uses Wilcox function in R)
wilcox.test(Prepub_HeightSDS_males, Pub_HeightSDS_males, conf.int = TRUE, conf.level = 0.95)

#Weight SDS
Prepub_WeightSDS_males <- clindata_prepub_males$WeightSDS
Pub_WeightSDS_males <- clindata_pub_males$WeightSDS
wilcox.test(Prepub_WeightSDS_males, Pub_WeightSDS_males, conf.int = TRUE, conf.level = 0.95)

#Testicular Volume
Prepub_TV_males <- clindata_prepub_males$TesticularVolume
Pub_TV_males<- clindata_pub_males$TesticularVolume
wilcox.test(Prepub_TV_males, Pub_TV_males, conf.int = TRUE, conf.level = 0.95)

#Age 
Prepub_Age_males <- clindata_prepub_males$Age_yrs
Pub_Age_males <- clindata_pub_males$Age_yrs
wilcox.test(Prepub_Age_males, Pub_Age_males, conf.int = TRUE, conf.level = 0.95)

#LH0min
Prepub_LH0min_males <- clindata_prepub_males$LH0min
Pub_LH0min_males<- clindata_pub_males$LH0min
wilcox.test(Prepub_LH0min_males, Pub_LH0min_males, conf.int = TRUE, conf.level = 0.95)

#LH30min
Prepub_LH30min_males <- clindata_prepub_males$LH30min
Pub_LH30min_males<- clindata_pub_males$LH30min
wilcox.test(Prepub_LH30min_males, Pub_LH30min_males,conf.int = TRUE, conf.level = 0.95)

#LH60min
Prepub_LH60min_males <- clindata_prepub_males$LH60min
Pub_LH60min_males<- clindata_pub_males$LH60min
wilcox.test(Prepub_LH60min_males, Pub_LH60min_males,conf.int = TRUE, conf.level = 0.95)

#FSH0min
Prepub_FSH0min_males <- clindata_prepub_males$FSH0min
Pub_FSH0min_males<- clindata_pub_males$FSH0min
wilcox.test(Prepub_FSH0min_males, Pub_FSH0min_males,conf.int = TRUE, conf.level = 0.95)

#FSH30min
Prepub_FSH30min_males <- clindata_prepub_males$FSH30min
Pub_FSH30min_males<- clindata_pub_males$FSH30min
wilcox.test(Prepub_FSH30min_males, Pub_FSH30min_males,conf.int = TRUE, conf.level = 0.95)

#FSH60min
Prepub_FSH60min_males <- clindata_prepub_males$FSH60min
Pub_FSH60min_males<- clindata_pub_males$FSH60min
wilcox.test(Prepub_FSH60min_males, Pub_FSH60min_males,conf.int = TRUE, conf.level = 0.95)

#Peak LH
Prepub_PeakLH_males <- clindata_prepub_males$PeakLH
Pub_PeakLH_males <- clindata_pub_males$PeakLH
wilcox.test(Prepub_PeakLH_males, Pub_PeakLH_males,conf.int = TRUE, conf.level = 0.95)

#Peak FSH
Prepub_PeakFSH_males <- clindata_prepub_males$PeakFSH
Pub_PeakFSH_males <- clindata_pub_males$PeakFSH
wilcox.test(Prepub_PeakFSH_males, Pub_PeakFSH_males,conf.int = TRUE, conf.level = 0.95)

#Basal LH/FSH Ratio
Prepub_BasalLHFSHratio_males <- clindata_prepub_males$BasalLHFSHratio
Pub_BasalLHFSHratio_males <- clindata_pub_males$BasalLHFSHratio
wilcox.test(Prepub_BasalLHFSHratio_males, Pub_BasalLHFSHratio_males,conf.int = TRUE, conf.level = 0.95)

#Peak LH/FSH Ratio
Prepub_PeakLHFSHratio_males <- clindata_prepub_males$PeakLHFSHratio
Pub_PeakLHFSHratio_males <- clindata_pub_males$PeakLHFSHratio
wilcox.test(Prepub_PeakLHFSHratio_males, Pub_PeakLHFSHratio_males,conf.int = TRUE, conf.level = 0.95)

#Testosterone Conc
Prepub_Testo_males <- clindata_prepub_males$TestosteroneConc
Pub_Testo_males <- clindata_pub_males$TestosteroneConc
wilcox.test(Prepub_Testo_males, Pub_Testo_males,conf.int = TRUE, conf.level = 0.95)

#Stats - females only----
#Fishers 
#Fishers 
#Fisher's Exact test - Sex & Puberty
data <- clindata_females[,c("Sex", "Group")]
data <- na.omit(data)
data <- table(data$Sex, data$Group)
data
fisher.test(data)

#Fisher's Exact test - Ethnicity & Puberty
data <- clindata_females[,c("Ethnicity", "Group")]
data <- na.omit(data)
data <- table(data$Ethnicity, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Breast Stage & Puberty
data <- clindata_females[,c("BreastStage", "Group")]
data <- na.omit(data)
data <- table(data$BreastStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pubic Hair Stage & Puberty
data <- clindata_females[,c("PubicHairStage", "Group")]
data <- na.omit(data)
data <- table(data$PubicHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Axillary Hair Stage & Puberty
data <- clindata_females[,c("AxillaryHairStage", "Group")]
data <- na.omit(data)
data <- table(data$AxillaryHairStage, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pituitary MR & Puberty
data <- clindata_females[,c("PituitaryMR", "Group")]
data <- na.omit(data)
data <- table(data$PituitaryMR, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Bone Age X-ray & Puberty
data <- clindata_females[,c("BoneAgeXray", "Group")]
data <- na.omit(data)
data <- table(data$BoneAgeXray, data$Group)
data
fisher.test(data)

#Fisher's Exact test  - Pelvic US & Puberty
data <- clindata_females[,c("PelvicUltrasound", "Group")]
data <- na.omit(data)
data <- table(data$PelvicUltrasound, data$Group)
data
fisher.test(data)

#Measured Variables
#NAs automatically omitted from test
#Height SDS
Prepub_HeightSDS_females <- clindata_prepub_females$HeightSDS
Pub_HeightSDS_females <- clindata_pub_females$HeightSDS
#Perform the Mann Whitney U test (uses Wilcox function in R)
wilcox.test(Prepub_HeightSDS_females, Pub_HeightSDS_females, conf.int = TRUE, conf.level = 0.95)

#Weight SDS
Prepub_WeightSDS_females <- clindata_prepub_females$WeightSDS
Pub_WeightSDS_females <- clindata_pub_females$WeightSDS
wilcox.test(Prepub_WeightSDS_females, Pub_WeightSDS_females, conf.int = TRUE, conf.level = 0.95)

#Age 
Prepub_Age_females <- clindata_prepub_females$Age_yrs
Pub_Age_females <- clindata_pub_females$Age_yrs
wilcox.test(Prepub_Age_females, Pub_Age_females, conf.int = TRUE, conf.level = 0.95)

#LH0min
Prepub_LH0min_females <- clindata_prepub_females$LH0min
Pub_LH0min_females<- clindata_pub_females$LH0min
wilcox.test(Prepub_LH0min_females, Pub_LH0min_females, conf.int = TRUE, conf.level = 0.95)

#LH30min
Prepub_LH30min_females <- clindata_prepub_females$LH30min
Pub_LH30min_females<- clindata_pub_females$LH30min
wilcox.test(Prepub_LH30min_females, Pub_LH30min_females,conf.int = TRUE, conf.level = 0.95)

#LH60min
Prepub_LH60min_females <- clindata_prepub_females$LH60min
Pub_LH60min_females<- clindata_pub_females$LH60min
wilcox.test(Prepub_LH60min_females, Pub_LH60min_females,conf.int = TRUE, conf.level = 0.95)

#FSH0min
Prepub_FSH0min_females <- clindata_prepub_females$FSH0min
Pub_FSH0min_females<- clindata_pub_females$FSH0min
wilcox.test(Prepub_FSH0min_females, Pub_FSH0min_females,conf.int = TRUE, conf.level = 0.95)

#FSH30min
Prepub_FSH30min_females <- clindata_prepub_females$FSH30min
Pub_FSH30min_females<- clindata_pub_females$FSH30min
wilcox.test(Prepub_FSH30min_females, Pub_FSH30min_females,conf.int = TRUE, conf.level = 0.95)

#FSH60min
Prepub_FSH60min_females <- clindata_prepub_females$FSH60min
Pub_FSH60min_females<- clindata_pub_females$FSH60min
wilcox.test(Prepub_FSH60min_females, Pub_FSH60min_females,conf.int = TRUE, conf.level = 0.95)

#Peak LH
Prepub_PeakLH_females <- clindata_prepub_females$PeakLH
Pub_PeakLH_females <- clindata_pub_females$PeakLH
wilcox.test(Prepub_PeakLH_females, Pub_PeakLH_females,conf.int = TRUE, conf.level = 0.95)

#Peak FSH
Prepub_PeakFSH_females <- clindata_prepub_females$PeakFSH
Pub_PeakFSH_females <- clindata_pub_females$PeakFSH
wilcox.test(Prepub_PeakFSH_females, Pub_PeakFSH_females,conf.int = TRUE, conf.level = 0.95)

#Basal LH/FSH Ratio
Prepub_BasalLHFSHratio_females <- clindata_prepub_females$BasalLHFSHratio
Pub_BasalLHFSHratio_females <- clindata_pub_females$BasalLHFSHratio
wilcox.test(Prepub_BasalLHFSHratio_females, Pub_BasalLHFSHratio_females,conf.int = TRUE, conf.level = 0.95)

#Peak LH/FSH Ratio
Prepub_PeakLHFSHratio_females <- clindata_prepub_females$PeakLHFSHratio
Pub_PeakLHFSHratio_females <- clindata_pub_females$PeakLHFSHratio
wilcox.test(Prepub_PeakLHFSHratio_females, Pub_PeakLHFSHratio_females,conf.int = TRUE, conf.level = 0.95)

#Oestradiol Conc
Prepub_Oest_females <- clindata_prepub_females$OestradiolConc
Pub_Oest_females <- clindata_pub_females$OestradiolConc
wilcox.test(Prepub_Oest_females, Pub_Oest_females,conf.int = TRUE, conf.level = 0.95)

#Tests for normality----
#Shapiro-Wilks test of normality
set.seed(1)
shapiro.test(clindata$HeightSDS)
shapiro.test(clindata$WeightSDS)
shapiro.test(clindata$Age_yrs)
shapiro.test(clindata$LH0min)
shapiro.test(clindata$FSH0min)
shapiro.test(clindata$LH30min)
shapiro.test(clindata$FSH30min)
shapiro.test(clindata$LH60min)
shapiro.test(clindata$FSH60min)
shapiro.test(clindata$BasalLHFSHratio)
shapiro.test(clindata$PeakLHFSHratio)
shapiro.test(clindata$TestosteroneConc)
shapiro.test(clindata$OestradiolConc)
         
shapiro.test(clindata_females$HeightSDS)
shapiro.test(clindata_females$WeightSDS)
shapiro.test(clindata_females$Age_yrs)
shapiro.test(clindata_females$LH0min)
shapiro.test(clindata_females$FSH0min)
shapiro.test(clindata_females$LH30min)
shapiro.test(clindata_females$FSH30min)
shapiro.test(clindata_females$LH60min)
shapiro.test(clindata_females$FSH60min)
shapiro.test(clindata_females$BasalLHFSHratio)
shapiro.test(clindata_females$PeakLHFSHratio)
shapiro.test(clindata_females$OestradiolConc)

shapiro.test(clindata_males$HeightSDS)
shapiro.test(clindata_males$WeightSDS)
shapiro.test(clindata_males$Age_yrs)
shapiro.test(clindata_males$LH0min)
shapiro.test(clindata_males$FSH0min)
shapiro.test(clindata_males$LH30min)
shapiro.test(clindata_males$FSH30min)
shapiro.test(clindata_males$LH60min)
shapiro.test(clindata_males$FSH60min)
shapiro.test(clindata_males$BasalLHFSHratio)
shapiro.test(clindata_males$PeakLHFSHratio)
shapiro.test(clindata_males$TestosteroneConc)



#Tidy data for thesis tables----
Descriptive_prepub_sub <- Descriptive_prepub[,c(2,5)]
Quartiles_prepub_sub <- Quartiles_prepub[c(2:4),]
Quartiles_prepub_sub <- t(Quartiles_prepub_sub)
Prepub_sub <- cbind(Descriptive_prepub_sub,Quartiles_prepub_sub)
colnames(Prepub_sub) <- paste0(colnames(Prepub_sub), "_Prepub")

Descriptive_pub_sub <- Descriptive_pub[,c(2,5)]
Quartiles_pub_sub <- Quartiles_pub[c(2:4),]
Quartiles_pub_sub <- t(Quartiles_pub_sub)
Pub_sub <- cbind(Descriptive_pub_sub,Quartiles_pub_sub)
colnames(Pub_sub) <- paste0(colnames(Pub_sub), "_Pub")

Descriptive_prepub_sub_females <- Descriptive_prepub_females[,c(2,5)]
Quartiles_prepub_sub_females <- Quartiles_prepub_females[c(2:4),]
Quartiles_prepub_sub_females <- t(Quartiles_prepub_sub_females)
Prepub_sub_females <- cbind(Descriptive_prepub_sub_females,Quartiles_prepub_sub_females)
colnames(Prepub_sub_females) <- paste0(colnames(Prepub_sub_females), "_Prepub_females")

Descriptive_pub_sub_females <- Descriptive_pub_females[,c(2,5)]
Quartiles_pub_sub_females <- Quartiles_pub_females[c(2:4),]
Quartiles_pub_sub_females <- t(Quartiles_pub_sub_females)
Pub_sub_females <- cbind(Descriptive_pub_sub_females,Quartiles_pub_sub_females)
colnames(Pub_sub_females) <- paste0(colnames(Pub_sub_females), "_Pub_females")

Descriptive_prepub_sub_males <- Descriptive_prepub_males[,c(2,5)]
Quartiles_prepub_sub_males <- Quartiles_prepub_males [c(2:4),]
Quartiles_prepub_sub_males  <- t(Quartiles_prepub_sub_males)
Prepub_sub_males  <- cbind(Descriptive_prepub_sub_males,Quartiles_prepub_sub_males)
colnames(Prepub_sub_males) <- paste0(colnames(Prepub_sub_males), "_Prepub_males")

Descriptive_pub_sub_males <- Descriptive_pub_males[,c(2,5)]
Quartiles_pub_sub_males <- Quartiles_pub_males[c(2:4),]
Quartiles_pub_sub_males <- t(Quartiles_pub_sub_males)
Pub_sub_males <- cbind(Descriptive_pub_sub_males,Quartiles_pub_sub_males)
colnames(Pub_sub_males) <- paste0(colnames(Pub_sub_males), "_Pub_males")

TableData <- cbind(Prepub_sub,Pub_sub,Prepub_sub_males,Pub_sub_males,Prepub_sub_females,Pub_sub_females)
TableData <- TableData[,c(-4, -9, -14, -19, -24,-29)]

new_order <- c(4, 5:16,3,1,2 )

TableData <- TableData[new_order, ]

write.csv(TableData, "TableData_all.csv")

