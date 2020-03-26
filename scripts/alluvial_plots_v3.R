library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(grid)
library(ggplotify)
library(colorspace)
library(gridExtra)
#library(ggpubr)

#Packages for synteny plots
#install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#install.packages("RIdeogram")
#install.packages("rsvg")
#library(RIdeogram)
#library(rsvg)

##################################################################################
cols = c("lncRNA"="springgreen4",
         "protein_coding"="steelblue4",
         "pseudogene" = "hotpink3",
         "lncRNA(h.c.)" = "springgreen4",
         "lncRNA(l.c.)" = "palegreen3",
         "pc(h.c.)" = "steelblue4",
         "pc(l.c.)" = "steelblue3",
         "pseudogene(l.c.)" = "hotpink3",
         "none" = "grey")

q <-  theme_minimal() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(size = 25, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"))

counts_table <- function(sp1, sp2, path = "./classification/"){
  
  file_name <- paste0(path, sp1, "to", sp2, ".classification")
  df_complete <- read.csv(file_name, sep = "\t", header = FALSE)
  
  df <- df_complete[,c(2,3)]
  counts <- as.data.frame(df %>% group_by_all() %>% summarise(COUNT = n()))
  counts$id <- counts$V3
  counts_long <- to_lodes_form(counts, key = "class", axes = 1:2)
  return(counts_long)
  
}

alluvial_plot <- function(sp1, sp2, path = "./classification/"){ 
  
  counts_long <- counts_table(sp1, sp2, path)
  p <- ggplot(data = counts_long,
              aes(x = class,
                  y = COUNT,
                  stratum = stratum,
                  alluvium = alluvium,
                  label = stratum)) +
    geom_flow(aes(fill=stratum), width = 1/20, alpha = 0.6) +
    geom_stratum(aes(fill=stratum), width = 1/20, size = 0.5) +
    scale_fill_manual(values = cols) +
    geom_label_repel(aes(fill=NULL),
                     stat = "stratum",
                     size = 5,
                     direction = "y",
                     nudge_x = 0.05,
                     hjust = -0.05,
                     point.padding = 0,
                     segment.size = 0.8,
                     show.legend = FALSE) +
    scale_x_discrete(labels=c(sp1, sp2)) +
    guides(fill="none")

    return(p+q+ggtitle("Orthology predictions", paste0(sp1, " to ", sp2)))
}


bar_plot <- function(sp1, sp2, path = "./classification/"){
  
  counts_long <- counts_table(sp1, sp2, path)
  counts_bar <- counts_long[counts_long$stratum =="lncRNA" & counts_long$class=="V2" ,]
  counts_bar$COUNT <- round(counts_bar$COUNT*100/sum(counts_bar$COUNT), digits = 2)
  p <- ggplot(data = counts_bar, aes(x = class, y=COUNT, label=paste0(COUNT,"%"), fill=id)) +
    geom_bar(aes(fill = id),
             position="stack",
             stat="identity",
             width=1) +
    scale_fill_manual(values = cols)+#c("springgreen4", "palegreen3", "steelblue4", "steelblue3", "grey50", "grey80")) +
    geom_text(position = position_stack(vjust = 0.5)) +
    guides(fill="none") +
    scale_x_discrete(labels = "lncRNA")
  
  return(p+q+ggtitle("1", "2")+theme(plot.title = element_text(colour = "white"),
                                     plot.subtitle = element_text(colour = "white")))
}
##################################################################################
args = commandArgs(trailingOnly=TRUE)
sp1=args[1]
sp2=args[2]

all <- alluvial_plot(sp1, sp2, path = "./classification/")
bar <- bar_plot(sp1, sp2, path = "./classification/")

plot <- grid.arrange(bar, all, ncol=2, widths=c(1, 9), nrow = 1)
plot_name <- paste0("./plots/", sp1, "to", sp2, ".plot", ".pdf")
ggsave(plot_name, plot, width=6, height=4, units="in", scale=3)

