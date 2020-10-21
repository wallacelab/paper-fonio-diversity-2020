#! /usr/bin/Rscript

# Plot principal coordinates of population

# Libraries
library(argparse)
library(ggplot2)
library(RColorBrewer)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="TASSEL distance matrix")
parser$add_argument("-q", "--popfile", help="File with taxon names and fastStructure population assignments")
parser$add_argument("-m", "--min-pop-identity", type='double', default=0.6, help="Minimum fraction of population assignment to get counted as that population")
parser$add_argument("-o", "--outprefix",  help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/2_Diversity.het01/')
# args=parser$parse_args(c("-i","2b_distances.txt", "-q", "2a_structure_calls.txt", "-o", "99_tmp"))


# Load data
cat("Plotting PCs for",args$infile,"\n")
dist=read.delim(args$infile, header=F, skip=5, row.names=1)
q=read.table(args$popfile, header=F, row.names=1)

# Clean distance matrix of any missing values (for very low het cutoffs)
if(any(is.na(dist))){
    cat("Removing NA values from the distance matrix\n")
    orig_count = nrow(dist)
    dist=as.matrix(dist)
    while(any(is.na(dist))){
        to.remove = which.max(rowSums(is.na(dist)))
        dist=dist[-to.remove, -to.remove]
    }
    cat("\tRemoved", orig_count - nrow(dist), "samples with missing data;", nrow(dist), "remain\n")
}

# Calculate PCs
rescaled = cmdscale(dist, eig=TRUE)
eig=rescaled$eig
mds = rescaled$points
colnames(mds)=c('PC1','PC2')

# Get best population assignment
names(q) = paste("Pop", 1:ncol(q), sep="")
best=apply(q, MARGIN=1, FUN=which.max)
best_score = apply(q, MARGIN=1, FUN=max)
q$best = names(q)[best]

# Make sure pass min fraction
q$best[best_score < args$min_pop_identity] = "none"

# Add pop to principal coordinates
mds=as.data.frame(mds)
mds$pop = q[rownames(mds), 'best']

# Make a color scheme
mypops = unique(mds$pop[mds$pop != "none"])
colors = brewer.pal(length(mypops), name="Set1")
names(colors) = mypops
colors = c(colors, "none" = "gray")

# Add color scheme to PCA data so can be used by later scripts
mds$color = colors[mds$pop]
mds$color[mds$pop=='none'] = 'gray'

# Plot results
percent_1 = round(eig[1]/sum(eig) * 100, digits=1)
percent_2 = round(eig[2]/sum(eig) * 100, digits=1)
myplot = ggplot(mds, mapping=aes(x=PC1, y=PC2, color=pop)) +
    geom_point(size=5, alpha=0.5) +
    labs(x=paste("PC1 ", percent_1, "%", sep=""), y=paste("PC2 ", percent_2, "%", sep=""), color="Population") +
    scale_color_manual(values=colors)
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=5, height=4)
ggsave(myplot, file=paste(args$outprefix, ".svg", sep=""), width=5, height=4)

# Write out PCs with pop assignments
write.table(mds, file=paste(args$outprefix, ".txt", sep=""), sep='\t', quote=F, row.names=T, col.names=T)
