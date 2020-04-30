#! /usr/bin/Rscript

# Convert a FASTA file of aligned SNPs into dendrograms

# Libraries
library(argparse)
library(ggplot2)
library(ggmap)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-k", "--keyfile", help="CSV keyfile of passport data for the different accessions")
parser$add_argument("-p", "--popkey", help="Population key (generated as part of making PCA plots)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/2_Diversity/')
# args=parser$parse_args(c("-k","../fonio.master_key.csv", "-o", "99_tmp", "-p", "2b_mds.txt"))


# Load data
key = read.csv(args$keyfile, row.names="Corrected.Sample.IDs..avoid.duplications.")
popdata = read.delim(args$popkey)

# Match samples together
popdata$sample = sub(rownames(popdata), pattern=".+/([^\\.]+)\\..+", repl="\\1")
popdata = cbind(popdata, key[popdata$sample,c("LATITUDE", "LONGITUDE")])

# Helper function to get x and y limits for plotting
get_limits = function(points, tolerance=0.1){
    lims = range(points, na.rm=TRUE)
    sep = lims[2] - lims[1]
    lims[1] = lims[1] - sep * tolerance
    lims[2] = lims[2] + sep * tolerance
    return(lims)
}

# Plot on map
toplot = unique(popdata[,c("pop", "LATITUDE", "LONGITUDE")])
xlim = get_limits(toplot$LONGITUDE, tolerance=0.25)
ylim = get_limits(toplot$LATITUDE, tolerance=1)
box = make_bbox(lon=xlim, lat=ylim)
terrain = get_map(location = box,  maptype = "terrain", source = "google", zoom = 8)
mymap = ggmap(terrain) + 
    geom_point(data=toplot, mapping=aes(x=LONGITUDE, y=LATITUDE, color=pop), size=2) +
    labs(x="Longitude", y="Latitude") +
    guides(color=FALSE)

# Save 
ggsave(mymap, file=paste(args$outprefix, ".png", sep=""))
ggsave(mymap, file=paste(args$outprefix, ".svg", sep=""))
