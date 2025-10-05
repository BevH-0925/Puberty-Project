library(qqman)
library(ggrepel)
library(tidyverse)

Kent <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Kent_genes.csv")
Hollis <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Hollis_genes.csv")
Thompson <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Thompson_genes.csv")
Bessa <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Bessa_genes.csv")
Palumbo <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Palumbo_genes.csv")
Chen <- read.csv("~/Documents/Code/ALSPAC_furtheranalysis/01_tidy_data/Chen_genes.csv")

RMCHMethylome <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/DMP.csv")
RMCHMethylomePred <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/DMP_Pred.csv")
RMCHTranscriptome<- read.csv("~/Documents/Code/RNASeqBH1-BH32/Random Forest/Allcombined.csv")

RMCHMethylomePred_genes <- RMCHMethylomePred[,c(1,15)]
RMCHMethylomePred_genes <- RMCHMethylomePred_genes[RMCHMethylomePred_genes$gene != "", ]

GWAS <- rbind(Kent,Hollis)
GWAS <- unique(GWAS)
Allstudies <- rbind(Kent,Hollis,Chen,Palumbo,Thompson, Bessa)
Allstudies <- unique(Allstudies)

RMCHMethylome <- RMCHMethylome[,c(1, 11, 12, 5)]
colnames(RMCHMethylome)[colnames(RMCHMethylome) == "X"] <- "SNP"
colnames(RMCHMethylome)[colnames(RMCHMethylome) == "MAPINFO"] <- "BP"
colnames(RMCHMethylome)[colnames(RMCHMethylome) == "P.Value"] <- "P"

RMCHMeth_annotateGWAS <- RMCHMethylomePred_genes
RMCHMeth_annotateGWAS$is_annotate <- NA

RMCHMeth_annotateGWAS <- RMCHMeth_annotateGWAS %>%
  mutate(is_annotate = if_else(gene %in% GWAS$Gene, "Yes", "No"))
colnames(RMCHMeth_annotateGWAS)[colnames(RMCHMeth_annotateGWAS) == "Probe"] <- "SNP"
colnames(RMCHMeth_annotateGWAS)[colnames(RMCHMeth_annotateGWAS) == "gene"] <- "Gene"

RMCHMeth_matchestudies <- RMCHMethylomePred_genes
RMCHMeth_matchestudies$is_highlight2 <- NA
RMCHMeth_matchestudies <- RMCHMeth_matchestudies %>%
  mutate(is_highlight2 = if_else(gene %in% Allstudies$Gene, "Yes", "No"))
colnames(RMCHMeth_matchestudies)[colnames(RMCHMeth_matchestudies) == "Probe"] <- "SNP"

RMCHMeth_nomatchstudies

RMCHMeth_nomatchstudies <- RMCHMethylomePred_genes
RMCHMeth_nomatchstudies$is_highlight1 <- NA
RMCHMeth_nomatchstudies <- RMCHMeth_nomatchstudies %>%
  mutate(is_highlight1 = if_else(gene %in% Allstudies$Gene, "No", "Yes"))
colnames(RMCHMeth_nomatchstudies)[colnames(RMCHMeth_nomatchstudies) == "Probe"] <- "SNP"

#Change chr column to just number not chr1 etc.
RMCHMethylome$CHR <- substr(RMCHMethylome$CHR, 4, nchar(RMCHMethylome$CHR))
RMCHMethylome <- RMCHMethylome[RMCHMethylome$CHR !="0",]
RMCHMethylome <- RMCHMethylome[RMCHMethylome$CHR !="M",]

RMCHMethylome$CHR <- as.numeric(RMCHMethylome$CHR)

data_cum <- RMCHMethylome |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

RMCHMethylome <- RMCHMethylome |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- RMCHMethylome |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- RMCHMethylome |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(RMCHMethylome)


#Put all data frames into list
df_list <- list(RMCHMethylome, RMCHMeth_nomatchstudies, RMCHMeth_matchestudies,RMCHMeth_annotateGWAS)

#Merge all data frames in list
RMCHMethylome2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS or puberty methylome studies
GWASPMmatching <- subset(RMCHMeth_matchestudies, is_highlight2=="Yes")
GWASPMmatching <- GWASPMmatching[,1]
RMCHMethylome2$SNP[RMCHMethylome2$SNP %in% GWASPMmatching] <- RMCHMethylome2$Gene[RMCHMethylome2$SNP %in% GWASPMmatching]

# Make the plot
ggplot(RMCHMethylome2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis

  
  # Add highlighted points
  geom_point(data=subset(RMCHMethylome2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(RMCHMethylome2, is_highlight2=="Yes"), color="orange", size=2) +
  
  
  #All genes from studies labelled
  geom_label_repel(min.segment.length = 0, data=subset(RMCHMethylome2, is_highlight2=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_RMCHMethylome_allmatchgeneslabelled.png", width = 10, height = 5, dpi=320))

qq(RMCHMethylome$P)
