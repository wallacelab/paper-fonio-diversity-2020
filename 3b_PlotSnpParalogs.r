#! /usr/bin/Rscript

# Convert a FASTA file of aligned SNPs into dendrograms

# Libraries
library(argparse)
library(ggplot2)
library(reshape2)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-a", "--a", help="Hapmap file of first set of SNPs")
parser$add_argument("-b", "--b", help="Hapmap file of second set of SNPs")
parser$add_argument("-p", "--pop-structure", help="Population structure results from step 2_")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/3_Paralogs/')
# args=parser$parse_args(c("-a","3a_snps15.hmp.txt", "-b", "3a_snps25.hmp.txt", "-o", "99_tmp", "-p", "../2_Diversity/2b_mds.txt"))

# Load data
cat("Comparing SNP sets from",args$a, "and",args$b, "\n")
a = read.delim(args$a, check.names=F)
b = read.delim(args$b, check.names=F)
pops = read.delim(args$pop_structure)
het_calls = c("K","M","R","S","W","Y") # IUPAC heterozygous codes

# Convert hapmap format to just genotype calls
hapmap_to_calls=function(x){
    hmp = x[, 12:ncol(x)]
    rownames(hmp) = x[,1]  # SNP name
    hmp = as.matrix(hmp)
    
    # Simplify missing and het calls
    hmp[hmp=="N"] = NA
    hmp[hmp %in% het_calls] = "Het"
    
    # Simplify other calls to just major/minor
    for(row in 1:nrow(hmp)){
        alleles = table(hmp[row,])
        alleles = alleles[names(alleles) %in% c("A","C","G","T")]
        alleles = sort(alleles, decreasing = TRUE)
        alleles = names(alleles)
        
        mysite = hmp[row,]
        # Major allele
        if(length(alleles) >= 1){
            hmp[row, mysite==alleles[1]] = "Major"
        }
        # Minor allele
        if(length(alleles) >= 2){
            hmp[row, mysite==alleles[2]] = "Minor"
        }
        # Set tertiary alleles to missing
        if(length(alleles) >= 3){
            hmp[row, mysite %in% alleles[3:length(alleles)]] = NA
        }
    }
    return(hmp)
}
a.calls = hapmap_to_calls(a)
b.calls = hapmap_to_calls(b)

# Make sure everything using the same samples
samples = intersect(colnames(a.calls), colnames(b.calls))
samples = intersect(samples, rownames(pops))
a.calls = a.calls[,samples]
b.calls = b.calls[,samples]
pops = pops[samples,]

# Find the SNPs that are unique to one of the hapmaps
all_sites = union(rownames(a.calls), rownames(b.calls))
common_sites = intersect(rownames(a.calls), rownames(b.calls))
unique_sites = setdiff(all_sites, common_sites)

# Make matrix of just the unique calls
a.unique = a.calls[rownames(a.calls) %in% unique_sites,]
b.unique = b.calls[rownames(b.calls) %in% unique_sites,]
unique.calls = rbind(a.unique, b.unique)

# Count % het calls in each pop (for later sorting)
percent_in_pop = function(mycalls, mypops, target_pop=NULL, count_minor=FALSE){
    # Make a matrix of counts
    mycount = matrix(0, nrow=nrow(mycalls), ncol=ncol(mycalls))
    mycount[is.na(mycalls)] = NA # Propagate missing values
    if(count_minor){  # Count occurances of minor alleles
        mycount[mycalls=="Minor"] = 1
        mycount[mycalls=="Het"] = 0.5
    } else{  # Count only hets
        mycount[mycalls=="Het"] = 1
    }
    
    # If given a target pop, subset down to just that pop
    if(!is.null(target_pop)){
        mycount = mycount[, mypops == target_pop]
        #cat("Subsetting to", target_pop,"\n")
    }
    
    # Return count
    return(rowSums(mycount, na.rm=T) / rowSums(!is.na(mycount)))
}
all_pops = sort(unique(pops$pop))
all_pops = all_pops[all_pops != "none"]
het.counts = lapply(all_pops, function(target_pop){
    percent_in_pop(unique.calls, mypops = pops$pop, target_pop = target_pop)
})
het.counts = as.data.frame(het.counts, col.names = all_pops)

# Sort sites by their het count in the different populations
site_order = order(het.counts$Pop1, het.counts$Pop2, het.counts$Pop3)
unique.calls = unique.calls[site_order,]

# Sort samples by their assigned population
if(!identical(rownames(pops), colnames(unique.calls))){ # Check that things match up
    stop("Sample names for calls and population assignments do not match!")
}
unique.calls = unique.calls[,order(pops$pop)]

################
# Actual graphic
################

# Add pop to sample names so ggplot sorts them correctly
tmp = unique.calls
colnames(tmp) = paste(pops[colnames(tmp), "pop"], colnames(tmp), sep="%")

# Reshape data into long form
plotdata = melt(tmp)
#plotdata = melt(tmp[1:20, ]) # DEBUGGING - REMOVE
names(plotdata) = c("site", "sample", "call")

# Extract population back out
plotdata$pop = factor(sub(plotdata$sample, pattern="%.+", repl=""))

# Set "Major" to missing (so doesn't add noise to plot)
to_keep=c('Minor', "Het")
plotdata$call[is.na(plotdata$call)] = "Missing"
plotdata$call[!plotdata$call %in% to_keep] = NA 
plotdata = subset(plotdata, !is.na(plotdata$call)) # Make smaller; less computation
plotdata$call = factor(as.character(plotdata$call), levels=c("Major", "Het", "Minor", "Missing"))
plotdata=droplevels(plotdata) # Cleaning factors up

# Make pop color scale
popcolors = unique(pops[,c("pop", "color")])
colorkey = popcolors$color
names(colorkey) = popcolors$pop

# Actual plot
myplot = ggplot(plotdata, mapping=aes(x=site, y=sample, fill=call)) +
    geom_raster() +
    geom_point(x=0.5, mapping=aes(y=sample, color=pop), 
               inherit.aes=FALSE, shape="-", size=6) +
    scale_y_discrete(limits=rev(levels(plotdata$sample))) +  # Get order of pops right
    coord_cartesian(clip='off') +  # So pop points don't get clipped
    labs(x="Sites (sorted)", y="Samples", fill="Call", color="Population") +
    scale_fill_brewer(palette="Dark2") +
    scale_color_manual(values=colorkey[levels(plotdata$pop)]) +
    theme_minimal() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(),
          axis.title = element_text(face='bold', size=14))
   
    
# Save data
write.csv(unique.calls, file=paste(args$outprefix, ".calls.csv", sep=""))
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=7, height=3)
