orthology_analysis <- function(file, sp1="", sp2="", path=""){
  library(dplyr)
  library(ggplot2)
  library(reshape)
  setwd(path)
  
  change_NAs <- function(table, value){
    for (name in c("e_class", "g_class")){
      a <- is.na(as.character(table[,name]))
      b <- as.character(table[,name])
      b[a] <- value
      table[,name] <- as.factor(b)
      
      a <- table[,name]=="."
      b <- as.character(table[,name])
      b[a] <- value
      table[,name] <- as.factor(b)
    }
    return(table)
  }
  ## READ orthology results
  orthology <- read.table(file, sep = "\t")
  colnames(orthology) <- c("gene", "e_class", "g_class", "e_ortho", "e_counts", "g_ortho", "g_counts")
  
  ## READ biotype of transcripts from sp1 # before geneName_geneType_geneID_map
  btype <- as.data.frame(read.table(paste0("./maps/",sp1,".geneID_geneName_geneType_map.txt")))
  colnames(btype) <- c("geneID", "gene", "Biotype")
  btype[,3] <- as.character(btype$Biotype)
  PC <- c("processed_transcript,protein_coding", "protein_coding", "protein_coding,lincRNA", "protein_coding,transcribed_unprocessed_pseudogene")
  lncRNA <- c("processed_transcript,transcribed_unprocessed_pseudogene", "processed_transcript,antisense", "processed_transcript,lincRNA", "snoRNA,lincRNA", "ribozyme,lincRNA", "misc_RNA,antisense", "TEC,lincRNA", "antisense,lincRNA", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense", "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna", "lncRNA")
  sncRNA <- c("scaRNA,miRNA", "misc_RNA,TEC", "snoRNA,miRNA", "misc_RNA,miRNA", "scRNA", "snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", "scaRNA", "vaultRNA")
  others <- c("translated_processed_pseudogene", "rRNA_pseudogene", "IG_J_pseudogene", "unprocessed_pseudogene,miRNA", "unprocessed_pseudogene,processed_pseudogene", "IG_LV_gene,IG_V_pseudogene", "TR_C_gene", "IG_C_gene", "IG_J_gene", "TR_D_gene", "IG_D_gene", "processed_pseudogene,lincRNA", "TR_J_gene", "IG_V_gene", "TR_V_gene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_LV_gene", "IG_V_pseudogene", "IG_pseudogene", "IG_V_pseudogene", "polymorphic_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "unprocessed_pseudogene", "unprocessed_pseudogene")
  novel <- c("NOVEL")
  
  ## MERGE sp1 biotypes with orthology results
  t <- merge(x = btype, y = orthology, by.x = "gene", by.y = "gene", all = TRUE)
  t2 <- t[,c(1,3,4,5)]
  
  ## SIMPLIFYING sp1 Biotype
  t2[t2$Biotype %in% PC,"Biotype"] <- "PC"
  t2[t2$Biotype %in% lncRNA,"Biotype"] <- "lncRNA"
  t2[t2$Biotype %in% sncRNA,"Biotype"] <- "sncRNA"
  t2[t2$Biotype %in% others,"Biotype"] <- "others"
  t2[t2$Biotype %in% novel,"Biotype"] <- "novel"
  t2 <- change_NAs(t2,"no_lifted")
  t2[is.na(t2$Biotype),"Biotype"] <- "novel"
  t2 <- melt(t2, id=c("gene","Biotype"))
  
  ## SIMPLIFYING classess from sp2
  lncRNA <- c("lncRNA", "lncRNAs", "lncRNA_stringtie", "lncRNA_sncRNA", "stringtie")
  sncRNA <- c("sncRNA", "sncRNAs")
  pc <- c("pc", "pcs")
  other <- c("other", "others") 
  t2[t2$value %in% lncRNA,"value"] <- "lncRNA"
  t2[t2$value %in% sncRNA,"value"] <- "sncRNA"
  t2[t2$value %in% pc,"value"] <- "pc"
  t2[t2$value %in% other,"value"] <- "other"
  
  ## INCLUDE SCORE for P-Index // lncRNA, sncRNA, lncRNA_other, pc, lncRNA_PC, other, none, no_lifted 
  a <- c("lncRNA") #1
  b <- c("lncRNA_other", "lncRNA_PC") #0.5
  c <- c("none") #0.25
  d <- c("sncRNA", "pc", "other", "no_lifted") #0.1
  t2$score <- NA #blank column 
  t2[t2$value %in% a,"score"] <- 1
  t2[t2$value %in% b,"score"] <- 0.5
  t2[t2$value %in% c,"score"] <- 0.25
  t2[t2$value %in% d,"score"] <- 0.1
  t2$value <- t2$value[,drop=TRUE]
  #levels(t2$value) <- unique(t2$value) #c("lncRNA","none","pc","other","lncRNA_PC","no_lifted","sncRNA","lncRNA_other")

  te <- t2[t2$variable=="e_class",]
  tg <- t2[t2$variable=="g_class",]
  write.table(te, file = paste0("classification/",sp1,"to",sp2,".classification.exons.csv"), quote = FALSE, sep = ",",
              dec = ".", row.names = FALSE, col.names = TRUE)
  write.table(tg, file = paste0("classification/",sp1,"to",sp2,".classification.genes.csv"), quote = FALSE, sep = ",",
              dec = ".", row.names = FALSE, col.names = TRUE)
              
  t2 <- t2[t2$Biotype!="others",]
  
  xflags = c("lncRNA", "sncRNA", "lncRNA_other", "pc", "lncRNA_PC", "other", "none", "no_lifted")
  dfl <- t2 %>% 
    group_by(Biotype,variable,value) %>% 
    summarise(n=n()) %>% 
    group_by(Biotype,variable) %>% 
    mutate(perc=100*n/sum(n))
  
  pdf(paste0("classification/",sp1,"to",sp2,".pdf")) 
  print(ggplot(dfl, aes(x=value, y=perc, color=variable, fill=value, width=.85)) +
    scale_fill_brewer(palette="YlGnBu") +
    scale_color_manual(values = c("black","darkgrey")) +
    geom_bar(stat="identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label=round(perc)), position=position_dodge(width = 0.9), vjust=-0.25) +
    ylab("percent") +
    scale_x_discrete(limits=xflags) +
    facet_wrap(~ Biotype) +
    theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()))
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)
file=args[1]
sp1=args[2]
sp2=args[3]
path=args[4]
orthology_analysis(file=file, sp1=sp1, sp2=sp2, path=path)
#orthology_analysis(file="Orthology/danrer11tohg38.classification.txt", sp1="danrer11", sp2="hg38", path="X:/p283/projects/carlos/barbara_orthologs/")
#orthology_analysis(file="Orthology/danrer11tomm10.classification.txt", sp1="danrer11", sp2="mm10", path="X:/p283/projects/carlos/barbara_orthologs/")

