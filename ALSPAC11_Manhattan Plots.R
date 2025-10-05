setwd("~/Documents/ALSPAC Methylome/Data/Manhattan")

library(qqman)
library(ggrepel)
library(tidyverse)

#cB2----
DMPcB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcB2_Mhtn.csv")
cB2_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cB2_annotateGWAS.csv")
cB2_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cB2_Borg.csv")
cB2_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cB2_Borgo.csv")

data_cum <- DMPcB2_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPcB2_Mhtn <- DMPcB2_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPcB2_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPcB2_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPcB2_Mhtn)

#Put all data frames into list
df_list <- list(DMPcB2_Mhtn, cB2_Borg, cB2_Borgo,cB2_annotateGWAS)

#Merge all data frames in list
DMPcB2_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(cB2_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPcB2_Mhtn2$SNP[DMPcB2_Mhtn$SNP %in% GWASmatching] <- DMPcB2_Mhtn2$Gene[DMPcB2_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPcB2_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
 
  # Add highlighted points
  geom_point(data=subset(DMPcB2_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +

  geom_point(data=subset(DMPcB2_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPcB2_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="B2_Cl", x="Chromosome", y=expression(-log[10](p-value)))
  (ggsave("Manhattan_B2_ClMeth_nolines.png"))

qq(DMPcB2_Mhtn$P)

#fB2----
DMPfB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfB2_Mhtn.csv")
fB2_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fB2_annotateGWAS.csv")
fB2_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fB2_Borg.csv")
fB2_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fB2_Borgo.csv")

data_cum <- DMPfB2_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPfB2_Mhtn <- DMPfB2_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPfB2_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPfB2_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPfB2_Mhtn)

#Put all data frames into list
df_list <- list(DMPfB2_Mhtn, fB2_Borg, fB2_Borgo, fB2_annotateGWAS)

#Merge all data frames in list
DMPfB2_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(fB2_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPfB2_Mhtn2$SNP[DMPfB2_Mhtn$SNP %in% GWASmatching] <- DMPfB2_Mhtn2$Gene[DMPfB2_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPfB2_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPfB2_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPfB2_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPfB2_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="B2_Fu", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_B2_FuMeth_nolines.png"))

qq(DMPfB2_Mhtn$P)

#cP2M----
DMPcP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2M_Mhtn.csv")
cP2M_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2M_annotateGWAS.csv")
cP2M_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2M_Borg.csv")
cP2M_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2M_Borgo.csv")

data_cum <- DMPcP2M_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPcP2M_Mhtn <- DMPcP2M_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPcP2M_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPcP2M_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPcP2M_Mhtn)

#Put all data frames into list
df_list <- list(DMPcP2M_Mhtn, cP2M_Borg, cP2M_Borgo, cP2M_annotateGWAS)

#Merge all data frames in list
DMPcP2M_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(cP2M_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPcP2M_Mhtn2$SNP[DMPcP2M_Mhtn$SNP %in% GWASmatching] <- DMPcP2M_Mhtn2$Gene[DMPcP2M_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPcP2M_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPcP2M_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPcP2M_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPcP2M_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2M_Cl", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2M_ClMeth_nolines.png"))

qq(DMPcP2M_Mhtn$P)

#fP2M----
DMPfP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2M_Mhtn.csv")
fP2M_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2M_annotateGWAS.csv")
fP2M_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2M_Borg.csv")
fP2M_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2M_Borgo.csv")

data_cum <- DMPfP2M_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPfP2M_Mhtn <- DMPfP2M_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPfP2M_Mhtn |>
  group_by(CHR) |>
  summarise(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPfP2M_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPfP2M_Mhtn)

#Put all data frames into list
df_list <- list(DMPfP2M_Mhtn, fP2M_Borg, fP2M_Borgo, fP2M_annotateGWAS)

#Merge all data frames in list
DMPfP2M_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(fP2M_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPfP2M_Mhtn2$SNP[DMPfP2M_Mhtn$SNP %in% GWASmatching] <- DMPfP2M_Mhtn2$Gene[DMPfP2M_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPfP2M_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPfP2M_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPfP2M_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPfP2M_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2M_Fu", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2M_FuMeth_nolines.png"))

qq(DMPfP2M_Mhtn$P)
#cP2F----
DMPcP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2F_Mhtn.csv")
cP2F_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2F_annotateGWAS.csv")
cP2F_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2F_Borg.csv")
cP2F_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cP2F_Borgo.csv")

data_cum <- DMPcP2F_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPcP2F_Mhtn <- DMPcP2F_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPcP2F_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPcP2F_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPcP2F_Mhtn)

#Put all data frames into list
df_list <- list(DMPcP2F_Mhtn, cP2F_Borg, cP2F_Borgo, cP2F_annotateGWAS)

#Merge all data frames in list
DMPcP2F_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(cP2F_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPcP2F_Mhtn2$SNP[DMPcP2F_Mhtn$SNP %in% GWASmatching] <- DMPcP2F_Mhtn2$Gene[DMPcP2F_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPcP2F_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPcP2F_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPcP2F_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPcP2F_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2F_Cl", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2F_ClMeth_nolines.png"))

qq(DMPcP2F_Mhtn$P)

#fP2F----
DMPfP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2F_Mhtn.csv")
fP2F_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2F_annotateGWAS.csv")
fP2F_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2F_Borg.csv")
fP2F_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fP2F_Borgo.csv")

data_cum <- DMPfP2F_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPfP2F_Mhtn <- DMPfP2F_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPfP2F_Mhtn |>
  group_by(CHR) |>
  summarise(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPfP2F_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPfP2F_Mhtn)

#Put all data frames into list
df_list <- list(DMPfP2F_Mhtn, fP2F_Borg, fP2F_Borgo, fP2F_annotateGWAS)

#Merge all data frames in list
DMPfP2F_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(fP2F_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPfP2F_Mhtn2$SNP[DMPfP2F_Mhtn2$SNP %in% GWASmatching] <- DMPfP2F_Mhtn2$Gene[DMPfP2F_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPfP2F_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPfP2F_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPfP2F_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPfP2F_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2F_Fu", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2F_FuMeth_nolines.png"))

qq(DMPfP2F_Mhtn$P)

#cMen----
DMPcMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcMen_Mhtn.csv")
cMen_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cMen_annotateGWAS.csv")
cMen_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cMen_Borg.csv")
cMen_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/cMen_Borgo.csv")

data_cum <- DMPcMen_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPcMen_Mhtn <- DMPcMen_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPcMen_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPcMen_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPcMen_Mhtn)

#Put all data frames into list
df_list <- list(DMPcMen_Mhtn, cMen_Borg, cMen_Borgo, cMen_annotateGWAS)

#Merge all data frames in list
DMPcMen_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(cMen_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPcMen_Mhtn2$SNP[DMPcMen_Mhtn$SNP %in% GWASmatching] <- DMPcMen_Mhtn2$Gene[DMPcMen_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPcMen_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPcMen_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPcMen_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPcMen_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="Men_Cl", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_Men_ClMeth_nolines.png"))

qq(DMPcMen_Mhtn$P)

#fMen----
DMPfMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfMen_Mhtn.csv")
fMen_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fMen_annotateGWAS.csv")
fMen_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fMen_Borg.csv")
fMen_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/fMen_Borgo.csv")

data_cum <- DMPfMen_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPfMen_Mhtn <- DMPfMen_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPfMen_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPfMen_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

sig <- 0.05 / nrow(DMPfMen_Mhtn)

#Put all data frames into list
df_list <- list(DMPfMen_Mhtn, fMen_Borg, fMen_Borgo, fMen_annotateGWAS)

#Merge all data frames in list
DMPfMen_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(fMen_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPfMen_Mhtn2$SNP[DMPfMen_Mhtn$SNP %in% GWASmatching] <- DMPfMen_Mhtn2$Gene[DMPfMen_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPfMen_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPfMen_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPfMen_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPfMen_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="Men_Fu", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_Men_FuMeth_nolines.png"))

qq(DMPfMen_Mhtn$P)

#B2----
DMPcB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcB2_Mhtn.csv")
DMPfB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfB2_Mhtn.csv")
DMPB2_Mhtn <- rbind(DMPcB2_Mhtn,DMPfB2_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPB2_Mhtn <- DMPB2_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

B2_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/B2_annotateGWAS.csv")
B2_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/B2_Borg.csv")
B2_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/B2_Borgo.csv")

data_cum <- DMPB2_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPB2_Mhtn <- DMPB2_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPB2_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPB2_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPB2_Mhtn)

#Put all data frames into list
df_list <- list(DMPB2_Mhtn, B2_Borg, B2_Borgo, B2_annotateGWAS)

#Merge all data frames in list
DMPB2_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(B2_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPB2_Mhtn2$SNP[DMPB2_Mhtn$SNP %in% GWASmatching] <- DMPB2_Mhtn2$Gene[DMPB2_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPB2_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPB2_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPB2_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPB2_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="B2", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_B2Meth_nolines.png"))

qq(DMPB2_Mhtn$P)

#P2M----
DMPcP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2M_Mhtn.csv")
DMPfP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2M_Mhtn.csv")
DMPP2M_Mhtn <- rbind(DMPcP2M_Mhtn,DMPfP2M_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPP2M_Mhtn <- DMPP2M_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

P2M_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2M_annotateGWAS.csv")
P2M_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2M_Borg.csv")
P2M_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2M_Borgo.csv")

data_cum <- DMPP2M_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPP2M_Mhtn <- DMPP2M_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPP2M_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPP2M_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPP2M_Mhtn)

#Put all data frames into list
df_list <- list(DMPP2M_Mhtn, P2M_Borg, P2M_Borgo, P2M_annotateGWAS)

#Merge all data frames in list
DMPP2M_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(P2M_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPP2M_Mhtn2$SNP[DMPP2M_Mhtn$SNP %in% GWASmatching] <- DMPP2M_Mhtn2$Gene[DMPP2M_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPP2M_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPP2M_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPP2M_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPP2M_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2M", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2MMeth_nolines.png"))

qq(DMPP2M_Mhtn$P)
#P2F----
DMPcP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2F_Mhtn.csv")
DMPfP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2F_Mhtn.csv")
DMPP2F_Mhtn <- rbind(DMPcP2F_Mhtn,DMPfP2F_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPP2F_Mhtn <- DMPP2F_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

P2F_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2F_annotateGWAS.csv")
P2F_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2F_Borg.csv")
P2F_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/P2F_Borgo.csv")

data_cum <- DMPP2F_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPP2F_Mhtn <- DMPP2F_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPP2F_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPP2F_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPP2F_Mhtn)

#Put all data frames into list
df_list <- list(DMPP2F_Mhtn, P2F_Borg, P2F_Borgo, P2F_annotateGWAS)

#Merge all data frames in list
DMPP2F_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(P2F_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPP2F_Mhtn2$SNP[DMPP2F_Mhtn$SNP %in% GWASmatching] <- DMPP2F_Mhtn2$Gene[DMPP2F_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPP2F_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPP2F_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPP2F_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPP2F_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="P2F", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_P2FMeth_nolines.png"))

qq(DMPP2F_Mhtn$P)

#Men----
DMPcMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcMen_Mhtn.csv")
DMPfMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfMen_Mhtn.csv")
DMPMen_Mhtn <- rbind(DMPcMen_Mhtn,DMPfMen_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPMen_Mhtn <- DMPMen_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

Men_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/Men_annotateGWAS.csv")
Men_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/Men_Borg.csv")
Men_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/Men_Borgo.csv")

data_cum <- DMPMen_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPMen_Mhtn <- DMPMen_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPMen_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPMen_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPMen_Mhtn)

#Put all data frames into list
df_list <- list(DMPMen_Mhtn, Men_Borg, Men_Borgo, Men_annotateGWAS)

#Merge all data frames in list
DMPMen_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(Men_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPMen_Mhtn2$SNP[DMPMen_Mhtn$SNP %in% GWASmatching] <- DMPMen_Mhtn2$Gene[DMPMen_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPMen_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPMen_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPMen_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPMen_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2, max.overlaps =20) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="Men", x="Chromosome", y=expression(-log[10](p-value)))
(ggsave("Manhattan_MenMeth_nolines.png"))

qq(DMPMen_Mhtn$P)

#Closest----
DMPcB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcB2_Mhtn.csv")
DMPcP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2M_Mhtn.csv")
DMPcP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2F_Mhtn.csv")
DMPcMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcMen_Mhtn.csv")
DMPc_Mhtn <- rbind(DMPcB2_Mhtn, DMPcP2M_Mhtn, DMPcP2F_Mhtn, DMPcMen_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPc_Mhtn <- DMPc_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

c_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/c_annotateGWAS.csv")
c_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/c_Borg.csv")
c_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/c_Borgo.csv")

data_cum <- DMPc_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPc_Mhtn <- DMPc_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPc_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPc_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPc_Mhtn)

#Put all data frames into list
df_list <- list(DMPc_Mhtn, c_Borg, c_Borgo, c_annotateGWAS)

#Merge all data frames in list
DMPc_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(c_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPc_Mhtn2$SNP[DMPc_Mhtn$SNP %in% GWASmatching] <- DMPc_Mhtn2$Gene[DMPc_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPc_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPc_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPc_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPc_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2,max.overlaps = 25) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="Closest Combined", x="Chromosome", y=expression(-log[10](p-value)))
ggsave("Manhattan_cMeth_nolines.png")

qq(DMPc_Mhtn$P)

#Furthest----
DMPfB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfB2_Mhtn.csv")
DMPfP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2M_Mhtn.csv")
DMPfP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2F_Mhtn.csv")
DMPfMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfMen_Mhtn.csv")
DMPf_Mhtn <- rbind(DMPfB2_Mhtn, DMPfP2M_Mhtn, DMPfP2F_Mhtn, DMPfMen_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPf_Mhtn <- DMPf_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

f_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/f_annotateGWAS.csv")
f_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/f_Borg.csv")
f_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/f_Borgo.csv")

data_cum <- DMPf_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPf_Mhtn <- DMPf_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPf_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPf_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#sig <- 0.05 / nrow(DMPf_Mhtn)

#Put all data frames into list
df_list <- list(DMPf_Mhtn, f_Borg, f_Borgo, f_annotateGWAS)

#Merge all data frames in list
DMPf_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(f_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPf_Mhtn2$SNP[DMPf_Mhtn$SNP %in% GWASmatching] <- DMPf_Mhtn2$Gene[DMPf_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPf_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPf_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPf_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPf_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2,max.overlaps = 25) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="Furthest Combined", x="Chromosome", y=expression(-log[10](p-value)))
ggsave("Manhattan_fMeth_nolines.png")

qq(DMPf_Mhtn$P)

#All combined----
DMPcB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcB2_Mhtn.csv")
DMPcP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2M_Mhtn.csv")
DMPcP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcP2F_Mhtn.csv")
DMPcMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPcMen_Mhtn.csv")
DMPfB2_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfB2_Mhtn.csv")
DMPfP2M_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2M_Mhtn.csv")
DMPfP2F_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfP2F_Mhtn.csv")
DMPfMen_Mhtn <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/DMPfMen_Mhtn.csv")

DMPall_Mhtn <- rbind(DMPcB2_Mhtn, DMPcP2M_Mhtn, DMPcP2F_Mhtn, DMPcMen_Mhtn, DMPfB2_Mhtn, DMPfP2M_Mhtn, DMPfP2F_Mhtn, DMPfMen_Mhtn)

#Remove duplicate cgs - keep lowest P value
DMPall_Mhtn <- DMPall_Mhtn %>%
  group_by(SNP) %>%
  filter(P== min(P)) %>%
  ungroup()

All_annotateGWAS <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/All_annotateGWAS.csv")
All_Borg <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/All_Borg.csv")
All_Borgo <- read.csv("~/Documents/ALSPAC Methylome/Data/Manhattan/All_Borgo.csv")

data_cum <- DMPall_Mhtn |>
  group_by(CHR) |>
  summarise(max_bp = max(BP)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

DMPall_Mhtn <- DMPall_Mhtn |>
  inner_join(data_cum, by = "CHR")|>
  mutate(bp_cum = BP + bp_add)

axis_set <- DMPall_Mhtn |>
  group_by(CHR) |>
  summarize(center = mean(bp_cum))

#Smaller y axis
ylim <- DMPall_Mhtn |>
  filter(P == min(P)) |>
  mutate(ylim = abs(floor(log10(P)))) |>
  pull(ylim)

#Put all data frames into list
df_list <- list(DMPall_Mhtn, All_Borg, All_Borgo, All_annotateGWAS)

#Merge all data frames in list
DMPall_Mhtn2 <- df_list %>% reduce(full_join, by='SNP')

#Rename any cgs with matching gene name from GWAS
GWASmatching <- subset(All_annotateGWAS, is_annotate=="Yes")
GWASmatching <- GWASmatching[,1]
DMPall_Mhtn2$SNP[DMPall_Mhtn$SNP %in% GWASmatching] <- DMPall_Mhtn2$Gene[DMPall_Mhtn2$SNP %in% GWASmatching]

# Make the plot
ggplot(DMPall_Mhtn2, aes(x=bp_cum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axis_set$CHR, breaks= axis_set$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +    # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(DMPall_Mhtn2, is_highlight1=="Yes"), color="violet", size=2) +
  
  geom_point(data=subset(DMPall_Mhtn2, is_highlight2=="Yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(min.segment.length = 0, data=subset(DMPall_Mhtn2, is_annotate=="Yes"), aes(label=SNP), size=2, max.overlaps = 45) +
  
  # Custom the theme:
  theme_classic() +
  theme( 
    legend.position="none",
  )+
  labs(title="All groups Combined", x="Chromosome", y=expression(-log[10](p-value)))
ggsave("Manhattan_allMeth_nolines.png")

qq(DMPall_Mhtn$P)