library(ggplot2)
library(gridExtra)
library(ggplotify)
library(dplyr)

##### Functions #####
read_df <- function(df_type, sp1=FALSE, sp2=FALSE){
  if (sp1 != FALSE){
    f <- paste0(sp1,"_", sp2, "_")
  }else{
    f <- ""
  }
  file_name1 <- paste0("counts/", f, "cluster_", df_type, "_exon", ".csv")
  file_name2 <- paste0("counts/", f, "cluster_", df_type, "_gene", ".csv")
  df1 <- read.csv(file = file_name1, sep = ",", header = TRUE)
  df2 <- read.csv(file = file_name2, sep = ",", header = TRUE)
  df1$Level <- "exon"
  df2$Level <- "gene"
  df <- rbind(df1, df2)
  return(df)
}
get_names <- function(df){
  sp_names <- levels(df$Species)
  return(c(sp_names))
}
lifted_nonlifted <- function(names){
  counter <- data.frame()
  for (sp1 in names){
    biotype_name <- paste0("maps/", sp1, ".geneID_geneName_geneType_map.txt")
    biotype <- read.csv(file = biotype_name, sep = "\t", header = FALSE)
    for (sp2 in names){
      if (sp1 == sp2){next}
      for (level in c("exon", "gene")){
        lifted_name <- paste0("liftovers/", sp1, "to", sp2, ".", level, "s.liftover")
        nonlifted_name <- paste0("liftovers/", sp1, "to", sp2, ".", level, "s.unmapped")
        
        lifted <- read.csv(file = lifted_name, sep = "\t", header = FALSE)
        lifted$level <- level
        lifted$lifted <- "Lifted"
        lifted <- merge(lifted[, c("V4", "level", "lifted")], biotype[,c("V2", "V3")], by.x = 1, by.y = 1)
        lifted <- lifted %>% distinct()
        
        nonlifted <- read.csv(file = nonlifted_name, sep = "\t", header = FALSE)
        nonlifted$level <- level
        nonlifted$lifted <- "NonLifted"
        nonlifted <- merge(nonlifted[, c("V4", "level", "lifted")], biotype[,c("V2", "V3")], by.x = 1, by.y = 1)
        ids_lifted <- unique(lifted$V4)
        nonlifted <- nonlifted[! nonlifted$V4 %in% ids_lifted, ]
        nonlifted <- nonlifted %>% distinct()
        
        counter_lifted <- data.frame(table(lifted[, c("level", "lifted", "V3")]))
        counter_nonlifted <- data.frame(table(nonlifted[, c("level", "lifted", "V3")]))
        counter_sp_level <- rbind(counter_lifted, counter_nonlifted)
        counter_sp_level$Sp <- sp1
        counter_sp_level$SpTo <- sp2
        
        counter <- rbind(counter, counter_sp_level)
      }
    }
  }
  return(counter)
}

cluster_stats_counts_all <- function(){
  df <- read_df("stats")
  df <- data.frame(table(df[,c("Cluster.type", "Biotypes", "Level")]))
  return(df)
}
cluster_genes_counts_all <- function(){
  df <- read_df("genes")
  df <- data.frame(table(df[,c("Species", "Biotype", "Cluster.type", "Level")]))
  return(df)
}

cluster_stats_counts_sps <- function(){
  files <- list.files("./counts/")[!startsWith(list.files("./counts/"), "cluster")]
  files_sps <- unique(lapply(files, function(x) strsplit(x, "_")[[1]][c(1,2)]))
  
  df_list <- list()
  for (f in files_sps){
    sp1 <- f[1]
    sp2 <- f[2]
    
    df <- read_df("stats", sp1, sp2)
    df <- data.frame(table(df[,c("Biotypes", "Cluster.type", "Level")]))
    
    pdf(paste0("./plots/clusters_stats_", sp1, "_", sp2, ".pdf"))
    p <- ggplot(data=df, aes(x=Cluster.type, y=Freq, fill = Level, label = Freq)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_manual(name = "Prediction Level", values=color_labels) +
      ggtitle(paste0("Clusters specificity: ", sp1, " - ", sp2)) +
      labs(y="Number of clusters", x = "") +
      geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
      facet_grid(Biotypes~., scales = "free", space='free_x') +
      plots_theme
    print(p)
    dev.off()
  }
}
cluster_genes_counts_sps <- function(){
  files <- list.files("./counts/")[!startsWith(list.files("./counts/"), "cluster")]
  files_sps <- unique(lapply(files, function(x) strsplit(x, "_")[[1]][c(1,2)]))
  
  df_list <- list()
  for (f in files_sps){
    sp1 <- f[1]
    sp2 <- f[2]

    df <- read_df("genes", sp1, sp2)
    df <- data.frame(table(df[,c("Species", "Biotype", "Cluster.type", "Level")]))
    
    pdf(paste0("./plots/clusters_genes_", sp1, "_", sp2, ".pdf"))
    p <- ggplot(data=df, aes(x=Cluster.type, y=Freq, fill = Species, label = Freq)) +
      geom_bar(stat="identity", position=position_dodge()) +
      scale_fill_manual(name = "Prediction Level", values=color_labels) +
      ggtitle(paste0("Gene's orthology type: ", sp1, " - ", sp2)) +
      labs(y="Number of genes", x = "") +
      geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
      facet_grid(Biotype~Level, scales = "free", space='free_x') +
      plots_theme
    p
    print(p)
    dev.off()
  }
}
#####################


##### Get data #####
cluster_genes_df <- read_df("genes")
sps_names <- get_names(cluster_genes_df)

counter_lifted <- lifted_nonlifted(sps_names)

cluster_stats_all <- cluster_stats_counts_all()
cluster_genes_all <- cluster_genes_counts_all()
####################


###### Plots ######
color_labels <- c("goldenrod1","darkseagreen4")
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

#All clusters stats and genes
pdf(paste0("./plots/", "clusters_stats_all.pdf"))
p <- ggplot(data=cluster_stats_all, aes(x=Cluster.type, y=Freq, fill = Level, label = Freq)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(name = "Prediction Level", values=color_labels) +
  ggtitle(paste0("Clusters specificity")) +
  labs(y="Number of clusters", x = "") +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), hjust=0.1, vjust=0.5, angle = 35) +
  facet_grid(Biotypes~., scales = "free", space='free_x') +
  plots_theme
print(p)
dev.off()

pdf(paste0("./plots/", "clusters_genes_all.pdf"))
p <- ggplot(data=cluster_genes_all, aes(x=Cluster.type, y=Freq, fill = Species, label = Freq)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(name = "Prediction Level", values=color_labels) +
  ggtitle(paste0("Gene's orthology type")) +
  labs(y="Number of genes", x = "") +
  facet_grid(Biotype~Level, scales = "free", space='free_x') +
  plots_theme +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), hjust=0.1, vjust=0.5, angle = 35)
print(p)
dev.off()

#Per pair os species clusters stats and genes
cluster_stats_counts_sps()
cluster_genes_counts_sps()






