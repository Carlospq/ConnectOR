#install.packages("ggalluvial")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(grid)
library(ggplotify)

#####################################################################################################
ggplot_function <- function(counts_long, feature){
  #dev.new(width = 3, height = 4, unit = "in")
  labels_n <- length(unique(as.data.frame(counts_long)[,5]))
  stratum_counts <- counts_long %>%
    group_by(stratum) %>%
    summarize(new_counts = sum(COUNT))
  p <- ggplot(data = counts_long,
         aes(x = BioType,
             stratum = stratum,
             alluvium = alluvium,
             fill = stratum,
             y = COUNT, label = stratum)) +
    geom_alluvium(aes(fill = id), width = 1/50, alpha = 0.6) +
    geom_stratum(color = "black", size = 0.1, width = 1/50)  +
    scale_fill_manual(values = c("NOVEL"="grey50",
                                 "protein_coding"="grey80",
                                 "lncRNA"="darkgreen",
                                 "lncRNA_other"="chartreuse3",
                                 "lncRNA_pc"="darkorange2",
                                 "lncRNAs"="greenyellow",
                                 "none"="gold",
                                 "other"="darkorange",
                                 "others"="darkorange4",
                                 "pc"="indianred1",
                                 "pc_other"="firebrick1",
                                 "pcs"="brown4")) +
    geom_label_repel(aes(fill=NULL),
                     stat = "stratum",
                     size = 5,
                     direction = "y",
                     nudge_x = c(rep(0,2), rep(1,labels_n-2)),
                     hjust = c(rep(1.01,2), rep(0,labels_n-2)),
                     vjust = c(rep(0.01,2), rep(0,labels_n-2)),
                     point.padding = 0,
                     segment.size = 0.8,
                     show.legend = FALSE) +
    theme_minimal() +
    ggtitle("Orthology predictions", paste0(sp1, " to ", sp2, " at ", feature," level")) +
    theme(
      plot.title = element_text(size = 25, hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5),
      panel.grid = element_blank(),
      axis.ticks.y = element_line(size= 0.5),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15)
    )
  
  print(p)
  grid.force()
  
  stratum_counts <- counts_long %>%
    group_by(stratum) %>%
    summarize(new_counts = sum(COUNT))
  kids <- childNames(grid.get("labelrepeltree", grep = TRUE))
  d <- do.call(rbind, lapply(kids, function(n) {
    x <- grid.get(n)
    data.frame(
      x = convertX(x$x, "native", valueOnly = TRUE),
      y = convertY(x$y, "native", valueOnly = TRUE),
      x.orig = convertX(x$x.orig, "native", valueOnly = TRUE),
      y.orig = convertY(x$y.orig, "native", valueOnly = TRUE)
    )
  }))
  d <- rbind(d[order(d[d$x<0.5,]$y.orig, decreasing = TRUE),],
             d[order(d[d$x>0.5,]$y.orig, decreasing = TRUE)+2,])
  p + geom_label_repel(data=stratum_counts,
                       mapping = aes(x=c(rep(sp1,2), rep(sp2,labels_n-2)),
                                     y=(d$y.orig*(sum(counts_long$COUNT)/2))+((d$y.orig-0.5)*(sum(counts_long$COUNT))*0.04),
                                     alluvium=NULL,
                                     label=new_counts),
                            direction = "y",
                            nudge_x = c(rep(0.1,2), rep(0.1,labels_n-2)),
                            hjust = c(rep(0.2,2), rep(0.3,labels_n-2)),
                            point.padding = 0,
                            segment.color = NA,
                            show.legend = FALSE)
}

alluvial_plot <- function(sp1, sp2, path = "./classification/"){
  
  # Read input file with orthology predictions
  # setwd("X:/p283/projects/carlos/ConnectOR.v0.01/")
  file_name <- paste0(path, sp1, "to", sp2, ".classification")
  df_complete <- read.csv(file_name, sep = "\t", header = FALSE)
  
  # Prepare DF for plotting
  for (i in c(3,4)){
    df <- df_complete[,c(2,i)]
    colnames(df) <- c(sp1, sp2)
    
    counts <- as.data.frame(df %>% group_by_all() %>% summarise(COUNT = n()))
    counts$id <- counts[,2]
    counts <- counts[counts[,sp1] %in% c("protein_coding", "NOVEL"), ]
    counts_long <- to_lodes_form(counts,
                                  key = "BioType",
                                  axes = 1:2)

    feature <- if (i==3) { "exon" } else { "gene" }
    p <- ggplot_function(counts_long, feature)
    
    plot_name <- paste0("./plots/", sp1, "to", sp2, ".classification.", feature, ".pdf")
    ggsave(plot_name, p, width=3, height=3, units="in", scale=3)
    
  }
}
#####################################################################################################
#sp1="hg38"
#sp2="mm10"

args = commandArgs(trailingOnly=TRUE)
sp1=args[1]
sp2=args[2]
alluvial_plot(sp1,sp2)
