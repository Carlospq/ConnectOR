library(ggplot2)
library(gridExtra)
library(ggplotify)
library(dplyr)
library(scales)
library(ggrepel)

##### Arguments #####
args <- commandArgs(trailingOnly = TRUE)
gene_level <- as.logical(args[1])

#####   Theme   #####
plots_theme <-  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),            #Title
        legend.position="top",                                        #Legenesd
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 60, hjust = 1), #Lables
        axis.text.y = element_text(size = 14, hjust = 1),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 14),
        strip.background = element_rect(colour = NA, fill=c(NA,NA))) 


##### Functions #####
read_df <- function(df_type, sp1=FALSE, sp2=FALSE, gl=gene_level){
  if (sp1 != FALSE){
    f <- paste0(sp1,"_", sp2, "_")
  }else{
    f <- ""
  }
  file_name1 <- paste0("counts/", f, df_type, "_exon", ".csv")
  df1 <- read.csv(file = file_name1, sep = ",", header = TRUE)
  df1$Level <- "exon"
  if (gl){
    file_name2 <- paste0("counts/", f, df_type, "_gene", ".csv")
    df2 <- read.csv(file = file_name2, sep = ",", header = TRUE)
    df2$Level <- "gene"
    
    df <- rbind(df1, df2)
  } else {
    df <- df1
  }

  df <- data.frame(table(df[,c("Species", "Biotype", "Cluster.type", "Level")]))
  levels(df$Biotype)[levels(df$Biotype)=="ncRNA"] <- "lncRNA"
  levels(df$Biotype)[levels(df$Biotype)=="pc"] <- "Protein coding"
  df <- df[df$Species!="",]
  df <- df[df$Biotype!="",]
  df <- df[df$Cluster.type!="",]

  df$Cluster.type <- factor(df$Cluster.type, levels = c("Not lifted", "One to none", "Many to many", "One to many", "One to half", "One to one"))
  df <- df %>%
    group_by(Level, Species, Biotype, Cluster.type) %>%
    summarise(Count = sum(Freq)) %>%
    mutate( Ratio = Count / sum(Count),
            pos = (cumsum(Ratio) - 0.5 * Ratio))
           #, label = percent(Ratio %>% round(2)))
  df$Cluster.type <- factor(df$Cluster.type, levels = c("One to one", "One to half", "One to many", "Many to many", "One to none", "Not lifted"))
  #df[df$Count==0, "Count"] <- ""
  df <- df[df$Count!=0, ]
  
  return(df)
}
prepare_ticks_df <- function(df){
  valid_types <- c("Many to many", "One to many", "One to half", "One to one")
  ticks_df <- df[df$Cluster.type %in% valid_types, c("Level", "Species", "Biotype", "Ratio")] %>%
    mutate( ticks = 1, #round(1-(sum(Ratio)), 1),
            ticks_label = round(sum(Ratio)*100,1))
  ticks_df <- unique(ticks_df[,-4])
  #ticks_df$ticks <- round(ticks_df$ticks, 0)
  ticks_df$ticks_label <- format(as.double(ticks_df$ticks_label), nsmall = 1)
  ticks_df[ticks_df == " 0.0"] <- ""
  ticks_df[ticks_df == " NaN"] <- ""
  ticks_df$Cluster.type = "One to one"
  ticks_df$Cluster.type <- factor(ticks_df$Cluster.type, levels = c("One to one", "One to half", "One to many", "Many to many", "One to none", "Not lifted"))
  return(ticks_df)
}
insert_minor <- function(major_labs, n_minor) {
  labs <- c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
  labs[1:(length(labs)-n_minor)]
}

##### Get data #####
genes_stats_df <- read_df("genes_stats")
ticks_df <- prepare_ticks_df(genes_stats_df)


#####   Plot   #####
pdf(paste0("./plots/", "gene_stats_all.pdf"), width=9, height=12)
p <- ggplot(genes_stats_df, aes(x=Species, y=Ratio, label=Count, fill=Cluster.type)) +
            geom_bar(stat='identity') +
            labs(y="Percent genes", x="Reference annotations") +
            # Scale & Colors
            scale_y_continuous(labels = insert_minor(seq(0,100,25),4), breaks = seq(0,1,length.out = 21)) +
            #scale_fill_brewer("Orthology", palette = "BuGn") +
            scale_fill_manual(values = c("#238B45", "#74C476", "#BAE4B3", "#EDF8E9", "Grey", "Darkgrey")) +
            # Facets (splitting the plot)
            facet_grid(Biotype~Level, scales = "free", space='free_x') +
            # Labels (Total number of genes per category)
            geom_label_repel(aes(x=Species, y=Ratio, label=Count, group=Biotype), 
                             position = position_stack(vjust = .5),
                             show.legend = FALSE) +
            # Labels for total percent of genes with predicted orthology
            geom_text(data = ticks_df, aes(x = Species, y = ticks, label=ticks_label,
                                           group=Biotype, color=Species, fill=Cluster.type), 
                                       hjust = -1.65,
                                       show.legend = FALSE,
                                       fontface = "bold",
                                       check_overlap = T) +
            plots_theme
print(p)
dev.off()
