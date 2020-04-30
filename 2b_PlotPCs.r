#! /usr/bin/Rscript

# Plot principal coordinates of population

# Libraries
library(argparse)
library(ggplot2)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="TASSEL distance matrix")
parser$add_argument("-q", "--popfile", help="File with taxon names and fastStructure population assignments")
parser$add_argument("-m", "--min-pop-identity", type='double', default=0.6, help="Minimum fraction of population assignment to get counted as that population")
parser$add_argument("-o", "--outprefix",  help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/2_Diversity/')
# args=parser$parse_args(c("-i","2b_distances.txt", "-q", "2a_structure_calls.txt", "-o", "99_tmp"))


# Load data
cat("Plotting PCs for",args$infile,"\n")
dist=read.delim(args$infile, header=F, skip=5, row.names=1)
q=read.table(args$popfile, header=F, row.names=1)

# Calculate PCs
mds=cmdscale(dist)
colnames(mds)=c('PC1','PC2')

# Get best population assignment
names(q) = paste("Pop", 1:ncol(q), sep="")
best=apply(q, MARGIN=1, FUN=which.max)
best_score = apply(q, MARGIN=1, FUN=max)
q$best = names(q)[best]

# Make sure pass min fraction
q$best[best_score < args$min_pop_identity] = "none"

# Plot results
mds=as.data.frame(mds)
mds$pop = q[rownames(mds), 'best']
myplot = ggplot(mds, mapping=aes(x=PC1, y=PC2, color=pop)) +
    geom_point(size=5, alpha=0.5) +
    labs(x="PC1", y="PC2", color="Population")
    
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=5, height=4)

# Write out PCs with pop assignments
write.table(mds, file=paste(args$outprefix, ".txt", sep=""), sep='\t', quote=F, row.names=T, col.names=T)
