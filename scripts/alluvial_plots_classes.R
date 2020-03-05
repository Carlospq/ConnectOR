library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(grid)
library(ggplotify)
library(colorspace)

#Packages for synteny plots
#install.packages("RIdeogram")
#install.packages("rsvg")
#library(RIdeogram)
#library(rsvg)

##################################################################################
new_levels <- c("IG_C_gene",
                "IG_V_gene",
                #"NOVEL",
		"lncRNA",
                "polymorphic_pseudogene",
                "protein_coding",
                "TR_V_gene",
                "class1",
                "class2",
                "class3",
                "class5",
                "class6",
                "lncRNA lncRNA",
                "lncRNA lncRNA_pc",
                "lncRNA lncRNAs",
                "lncRNA_other lncRNA_other",
                "lncRNA_pc lncRNA_pc",
                "lncRNAs lncRNA_pc",
                "lncRNAs lncRNAs",
                "pc lncRNA_pc",
                "pc pc",
                "pc pc_other",
                "pc pcs",
                "pc_other lncRNA_pc",
                "pc_other pc_other",
                "pcs lncRNA_pc",
                "pcs pcs",
                "none lncRNA",
                "none lncRNA_pc",
                "none lncRNAs",
                "none none",
                "none pc",				
                "none pcs",
                "other other",
                "others others")

new_scale <- c(#"NOVEL"="grey50",
  "lncRNA"="grey50",
  "protein_coding"="grey80",
  "class1"="darkgreen",
  "class2"="greenyellow",
  "class3"="gold",
  "class5"="darkorange",
  "class6"="indianred1",
  "lncRNA lncRNA"="#023FA5",
  "lncRNA lncRNA_pc"="#3B53A6",
  "lncRNA lncRNAs"="#5666AB",
  "lncRNA_other lncRNA_other"="#6C79B3",
  "lncRNA_pc lncRNA_pc"="#818ABB",
  "lncRNAs lncRNA_pc"="#939BC2",
  "lncRNAs lncRNAs"="#A5AAC9",
  "none lncRNA"="#B5B8D0",
  "none lncRNA_pc"="#C3C5D6",
  "none lncRNAs"="#CFD1DB",
  "none none"="#D9DADF",
  "none pc"="#E1E1E2",
  "none pcs"="#E2E0E1",
  "other other"="#E0D9DA",
  "others others"="#DCCED0",
  "pc lncRNA_pc"="#D8C1C4",
  "pc pc"="#D2B1B7",
  "pc pc_other"="#CCA0A8",
  "pc pcs"="#C58E98",
  "pc_other lncRNA_pc"="#BC7A87",
  "pc_other pc_other"="#B26575",
  "pcs lncRNA_pc"="#A84D63",
  "pcs pcs"="#9B324F")

ggplot_function <- function(counts_long){
  ggplot(data = counts_long,
         aes(x = class,
             stratum = stratum,
             alluvium = alluvium,
             fill = stratum,
             y = COUNT, label = stratum)) +
    geom_flow(aes(fill=as.factor(id)), width = 1/50, alpha = 0.6) +
    scale_fill_manual(values = new_scale) +
    geom_stratum(color = "black", size = 0.1, width = 1/50) +
    geom_label_repel(aes(fill=NULL),
                     stat = "stratum",
                     size = 5,
                     direction = "y",
                     nudge_x = 0,
                     hjust = -1.01,
                     point.padding = 0,
                     segment.size = 0.8,
                     show.legend = FALSE)
}

alluvial_plot <- function(sp1, sp2, path = "./classification/"){
  #setwd("X:/p283/projects/carlos/ConnectOR.v0.01/")
  file_name <- paste0(path, sp1, "to", sp2, ".classification.classes")
  df_complete <- read.csv(file_name, sep = "\t", header = TRUE)
  
  df <- df_complete[,c(2,9,10)]
  counts <- as.data.frame(df %>% group_by_all() %>% summarise(COUNT = n()))
  counts$id <- as.factor(counts[,2])
  #counts <- counts[counts[,1] %in% c("protein_coding", "NOVEL"), ]
  counts$sp1_btype <- as.character(counts$sp1_btype)
  counts[counts$sp1_btype %in% c("antisense", "lincRNA", "processed_transcript", "sense_intronic", "sense_overlapping"), "sp1_btype"] = "lncRNA"
  counts$sp1_btype <- as.factor(counts$sp1_btype)
  counts <- counts[counts[,1] %in% c("protein_coding", "lncRNA"), ]
  counts_long <- to_lodes_form(counts,
                               key = "class",
                               axes = 1:3)
  labels_n <- length(unique(as.data.frame(counts_long)[,5]))
  stratum_counts <- counts_long %>%
    group_by(stratum) %>%
    summarize(new_counts = sum(COUNT))
  counts_long$stratum <- factor(counts_long$stratum, levels = new_levels)
  
  p <- ggplot_function(counts_long)
  plot_name <- paste0("./plots/", sp1, "to", sp2, ".classification.classes", ".pdf")
  ggsave(plot_name, p, width=6, height=4, units="in", scale=3)

  p <- bar_plot(sp1, sp2)
  plot_name <- paste0("./plots/", sp1, "to", sp2, ".classification.hist", ".pdf")
  ggsave(plot_name, p, width=6, height=4, units="in", scale=3)
}

bar_plot <- function(sp1, sp2, path = "./classification/"){
  #setwd("X:/p283/projects/carlos/ConnectOR.v0.01/")
  file_name <- paste0(path, sp1, "to", sp2, ".classification.classes")
  df_complete <- read.csv(file_name, sep = "\t", header = TRUE)
  
  df <- df_complete[,c(2,3,4)]
  counts <- as.data.frame(df %>% group_by_all() %>% summarise(COUNT = n()))
  counts$sp1_btype <- as.character(counts$sp1_btype)
  counts[counts$sp1_btype %in% c("antisense", "lincRNA", "processed_transcript", "sense_intronic", "sense_overlapping"), "sp1_btype"] = "lncRNA"
  counts$sp1_btype <- as.factor(counts$sp1_btype)

  sp1lncRNAs <- counts[counts$sp1_btype=="lncRNA",]
	
  sp2none      <- sp1lncRNAs$eclass=="none" & sp1lncRNAs$gclass=="none"
  sp2hc_lncRNA <- sp1lncRNAs$eclass %in% c("lncRNA", "lncRNAs", "lncRNA_other", "lncRNA_sncRNA")
  sp2lc_lncRNA <- sp1lncRNAs$eclass=="none" & sp1lncRNAs$gclass %in% c("lncRNA", "lncRNAs", "lncRNA_other", "lncRNA_sncRNA")
  sp2hc_pc     <- sp1lncRNAs$eclass %in% c("pc", "pcs", "lncRNA_PC")
  sp2lc_pc     <- sp1lncRNAs$eclass=="none" & sp1lncRNAs$gclass %in% c("pc", "pcs", "lncRNA_PC")
  others       <- !(sp2none | sp2hc_lncRNA | sp2lc_lncRNA | sp2hc_pc | sp2lc_pc)

  count_sp2none      <- sum(sp1lncRNAs[sp2none,]$COUNT)
  count_sp2hc_lncRNA <- sum(sp1lncRNAs[sp2hc_lncRNA,]$COUNT)
  count_sp2lc_lncRNA <- sum(sp1lncRNAs[sp2lc_lncRNA,]$COUNT)
  count_sp2hc_pc     <- sum(sp1lncRNAs[sp2hc_pc,]$COUNT)
  count_sp2lc_pc     <- sum(sp1lncRNAs[sp2lc_pc,]$COUNT)
  count_others       <- sum(sp1lncRNAs[others,]$COUNT)
  
  df_hist <- data.frame("Predictions"=c("hc_lncRNA", "lc_lncRNA", "hc_pc", "lc_pc", "others", "none"))
  df_hist$Predictions <- factor(x = c("hc_lncRNA", "lc_lncRNA", "hc_pc", "lc_pc", "others", "none"),
                            levels = c("hc_lncRNA", "lc_lncRNA", "hc_pc", "lc_pc", "others", "none"))

  df_hist$names <- "Counts"
  df_hist$value <- c(count_sp2hc_lncRNA, count_sp2lc_lncRNA, count_sp2hc_pc, count_sp2lc_pc, count_others, count_sp2none)
  df_hist$value <- round(df_hist$value*100/sum(df_hist$value), digits = 1)

  ggplot(data = df_hist, aes(x = names, y=value, label=paste0(value,"%"), fill=Predictions)) +
		geom_bar(aes(fill = Predictions),
		           stat="identity",
        		   width=0.2) +
  		scale_fill_manual(values = c("springgreen4", "palegreen3", "steelblue4", "steelblue3", "grey50", "grey80")) +
		geom_text(position = position_stack(vjust = 0.5))
}
##################################################################################
args = commandArgs(trailingOnly=TRUE)
sp1=args[1]
sp2=args[2]
alluvial_plot(sp1,sp2)





                                


