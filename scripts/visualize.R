#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)

if (length(args)<3) {
	stop("Usage: Rscript --vanilla visualize.R [method] [input folder] [output folder].\n", call.=FALSE)
} 
	
method <- args[1]
dir_in <- args[2]
dir_out <- args[3]

#method <- 'triadic-backup'
#method <- 'algorithm8'
#method <- 'jain-neal'

if (method == 'triadic') {
	method_title = "Triadic MCMC sampler"
} else if (method == 'algorithm8') {
	method_title = "Auxiliary variable MCMC sampler"
} else if (method == 'jain-neal') {
	method_title = "Jain/Neal Split-Merge MCMC sampler"
} else {
	method_title = "Unknown"
	print("Unknown method");
}

wpath <- file.path(dir_in)

purity <- read.table(file.path(wpath, 'purity.txt'), header=FALSE)
rand <- read.table(file.path(wpath, 'rand.txt'), header=FALSE)
adjusted_rand <- read.table(file.path(wpath, 'adjusted_rand.txt'), header=FALSE)

dataframe <- data.frame(
	purity,
	rand,
	adjusted_rand
)
names(dataframe) <- c("Purity", "Rand Index", "Adjusted Rand Index")

dataframe.m <- reshape2::melt(dataframe, id.vars = NULL)

bold.text <- element_text(face = "bold")

pdf(NULL)

ggplot(dataframe.m, aes(x = variable, y = value)) + 
	geom_violin() + 
	ggtitle(paste("Line estimation with the", method_title, sep=" ")) +
	scale_x_discrete(name = "Clustering metric") + 
	scale_y_continuous(name = "Clustering performance [0..1]", limits = c(0,1) ) +
	theme(title = bold.text, axis.title = bold.text) + 
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_boxplot(width=.1, outlier.size=0)

dir.create(file.path(dir_out), showWarnings = FALSE)
ggsave(file.path(dir_out, paste(method, "png", sep=".")))
