#!/usr/bin/env R
library(data.table)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage Rscript plot_classes.R DIR TITLE SUBTITLE", call. = FALSE)
} else {
  dir <- args[1]
  title <- args[2] # "bigTranscriptome ConnectOR predictions"
  subtitle <- args[3] # "Danio rerio to Homo Sapiens"
}

getClassificationStats <- function (dir = "./classification") {
  files <- list.files(path = dir, pattern = ".classification", full.names = TRUE)
  stats <- lapply(
    X = files, 
    FUN = function (x) {
      specie <- gsub(".*/([a-z]+[0-9]+)to.*.classification$", "\\1", x)
      dt <- fread(file = x, header = FALSE, sep = "\t")
      colnames(dt) <-
        c(
          "gene_name",
          "unknown1",
          "prediction",
          "unknown2",
          "unknown3",
          "unknown4",
          "unknown5",
          "unknown6",
          "unknown7"
        )
      dt[, .(count = .N, frequency = .N / nrow(dt), specie = specie),  by = c("prediction")]
    }
  )
  rbindlist(stats)
}

barplot_data <- getClassificationStats(dir)
ggplot(data = barplot_data, aes(x=prediction, y=frequency, fill=prediction)) +
  geom_bar(stat="identity", width = 0.5) + 
  geom_text(
    mapping = aes(label = count),
    position = position_dodge(width = 0.5),
    vjust = -0.25
  ) +
  labs(title = title, subtitle = subtitle) +
  ylab("Percentage of genes") + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~specie)

filename <- paste0(tolower(gsub(" ", "_", paste(title, subtitle))), "_bar_plot.png")
ggsave(filename = paste0(dir, filename))