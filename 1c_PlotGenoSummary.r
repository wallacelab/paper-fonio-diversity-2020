#! /usr/bin/Rscript

# Test for distortion in the kmers made from each sequence

# Libraries
library(argparse)
library(ggplot2)
library(gridExtra)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-s", "--sitefile", help="TASSEL genotype summary file of sites")
parser$add_argument("-t", "--taxafile", help="TASSEL genotype summary file of taxa")
parser$add_argument("-d", "--depthfile", help="VCFTools depth output file")
parser$add_argument("-o", "--outfile", help="Output graphic")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/PearlMillet/BoubacarPearlMilletPops/GenomeLinkageMaps/2_Analysis/')
# args=parser$parse_args(c("-s","2a_sitesummary.txt", "-t", "2a_taxasummary.txt", "-d", "2a_depth.txt", "-o", "99_tmp.png"))


plots = list()


###
# Depth summary 
###


if(!is.null(args$depthfile)){
    
    cat("Processing depth data from",args$depthfile, "\n")
    depthvals = read.delim(args$depthfile, check.names=F)$INFO
    depthvals = sort(depthvals)
    
    # Depth - TODO: As histogram or as dot plot? (Prob dots to match % het)
    depth_raw = qplot(1:length(depthvals), depthvals, color=I('darkgoldenrod')) +
        labs(title = "Sites - Read Depth", x="Site # (sorted)", y="Read Depth") 

    depth_log = qplot(1:length(depthvals), log10(depthvals), color=I('darkgoldenrod4')) +
        labs(title = "Sites - Read Depth (Log scale)", x="Site # (sorted)", y="Log10 Read Depth") 
    
    
    plots = c(plots, list(depth_raw, depth_log))
    
}


###
# Site summary
###

if(!is.null(args$sitefile)){
    
    cat("Processing site data from",args$sitefile, "\n")
    
    sites = read.delim(args$sitefile, check.names=F)
    
    # Percent Heterozygous; for some reason qplot() was not handling this properly so did full ggplot()
    hetvals = sort(sites[,'Proportion Heterozygous'])
    hetdata=data.frame(xval=1:length(hetvals), yval=hetvals)
    hets = ggplot(hetdata, mapping=aes(x=xval, y=yval)) +
        geom_point(color="darkred") +
        labs(title="Sites - Percent Heterozygous", x="Sites (sorted)", y="Fraction Het")
    
    # Percent Missing
    missing = qplot(sites[,'Proportion Missing'], bins=50, fill=I('darkgreen')) +
        labs(title="Sites - Percent Missing", x="Fraction Missing", y="# Sites")
    
    # Minor allele frequency
    maf = qplot(x=sites[,'Minor Allele Frequency'], bins=50, fill=I('darkblue')) +
        labs(title="Sites - Minor Allele Frequency", x="MAF", y="# Sites")
    

    plots = c(plots, list(hets, missing, maf))
        
}


###
# Taxa summary 
###

if(!is.null(args$taxafile)){
    
    cat("Processing taxa data from",args$sitefile, "\n")
    
    taxa = read.delim(args$taxafile, check.names=F)
    
    # Percent Heterozygous
    hetvals = sort(taxa[,'Proportion Heterozygous'])
    hets = qplot(hetvals, bins=50, fill=I('coral')) +
        labs(title = "Taxa - Percent Heterozygous", x="Fraction Het", y="# Taxa") 
        
    # Percent Missing
    missing = qplot(taxa[,'Proportion Missing'], bins=50, fill=I('aquamarine3')) +
        labs(title="Taxa - Percent Missing", x="Fraction Missing", y="# Taxa")
    
    
    
    
    plots = c(plots, list(hets, missing))
    
}


###
# Output
###

nplots = length(plots)
ncol=2
nrow = ceiling(nplots/ncol)

cat("Writing graphic to",args$outfile, "\n")
png(args$outfile, height=3*nrow, width=4*ncol, units="in", res=150)
    grid.arrange(grobs=plots, ncol=ncol, nrow=nrow, as.table=TRUE)
dev.off()
