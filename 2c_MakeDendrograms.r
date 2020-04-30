#! /usr/bin/Rscript

# Convert a FASTA file of aligned SNPs into dendrograms

# Libraries
library(argparse)
library(ape)
library(RColorBrewer)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="TASSEL Distance matrix file")
parser$add_argument("-p", "--popkey", help="Population key")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/2_Diversity/')
# args=parser$parse_args(c("-i","2b_distances.txt", "-o", "99_tmp", "-p", "2b_mds.txt"))


# Load data
cat("Loading distance matrix from", args$infile, "\n")
dist = read.delim(args$infile, skip=5, row.names=1, header=F)
key = read.delim(args$popkey)
popkey = key$pop
names(popkey) = rownames(key)


# Prettify sample names
samples = rownames(dist)
samples = sub(samples, pattern=".+/", repl="")
samples = sub(samples, pattern=".bam$", repl="")
samples = sub(samples, pattern="\\..+", repl="")
key = samples
names(key) = rownames(dist)

# Calculate tree
tree = nj(as.dist(dist))

# Update tip labels and colors
pops = popkey[tree$tip.label]
tree$tip.label = key[tree$tip.label]

# Write text tree output
write.tree(tree, file=paste(args$outprefix, ".tre", sep=""))

# Plot tree output
colors = brewer.pal(length(unique(pops)), name="Set1")
colorkey = colors[as.numeric(factor(pops))]
colorkey[pops=='none'] = 'gray'

png(paste(args$outprefix, ".png", sep=""), width=5, height=5, units="in", res=300)
    plot(tree, "unrooted", cex=0.5, edge.width=0.5, lab4ut='axial', tip.color=colorkey)
dev.off()

svg(paste(args$outprefix, ".svg", sep=""), width=5, height=5)
    plot(tree, "unrooted", cex=0.5, edge.width=0.5, lab4ut='axial', tip.color=colorkey)
dev.off()