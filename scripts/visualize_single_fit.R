#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
source("read.octave.R")

if (length(args) < 2) {
	stop("Usage: Rscript --vanilla visualize_single_fit.R [input folder] [output folder].\n", call.=FALSE)
} 
	
dir_in <- args[1]
dir_out <- args[2]

bname <- basename(dir_in)
wpath <- file.path(dir_in)

results <- read.octave(file.path(wpath, 'results.txt'), options.skip_line=0)

dataframe <- data.frame(
	results
)

mu <- data.matrix( dataframe[c(1,2)] )

K <- nrow(mu)

index <- 1:K

# Get all data points from 1 to K from the individual files
datalist = list()
for (i in index) {
	j <- i - 1
	fname = paste('results', j, '.txt', sep='')
	dat <- read.table(file.path(wpath, fname))
	dat$index <- i
	datalist[[i]] <- dat
}

ldata = do.call(rbind, datalist)

mu_df <- data.frame ( 
	mu,
	index
)

names(mu_df) <- c("intercept", "slope")

pdf(NULL)

p <- ggplot(data=ldata, aes(x=V2, y=V3, color=factor(index))) + 
	geom_point() + xlab("") + ylab("") + 
	theme(legend.position="none") + 
	geom_abline(data=mu_df, aes(intercept=intercept, slope=slope, color=factor(index)))

bold.text <- element_text(face = "bold")

dir.create(file.path(dir_out), showWarnings = FALSE)
ggsave(file.path(dir_out, paste("fit_", bname, ".png", sep="")))
