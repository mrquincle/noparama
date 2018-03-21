#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
source("read.octave.R")

if (length(args) < 2) {
	stop("Usage: Rscript --vanilla visualize_input.R [input file] [output folder].\n", call.=FALSE)
} 
	
file_in <- args[1]
dir_out <- args[2]

bname <- basename(file_in)

datapoints <- read.table(file_in, sep="\t")

dataframe <- data.frame(
	datapoints
)

pdf(NULL)

p <- ggplot(data=dataframe, aes(x=V1, y=V2, color=factor(V3))) + 
	geom_point() + xlab("") + ylab("") +
	theme(legend.position="none") 

bold.text <- element_text(face = "bold")

dir.create(file.path(dir_out), showWarnings = FALSE)
ggsave(file.path(dir_out, paste("original_", bname, ".png", sep="")))
