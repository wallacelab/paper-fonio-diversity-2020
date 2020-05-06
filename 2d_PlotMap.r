#! /usr/bin/Rscript

# Convert a FASTA file of aligned SNPs into dendrograms

# Libraries
library(argparse)
library(ggplot2)
library(dplyr)
library(ggmap)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-k", "--keyfile", help="CSV keyfile of passport data for the different accessions")
parser$add_argument("-p", "--popkey", help="Population key (generated as part of making PCA plots)")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-z", "--zoom", type='integer', default=7, help="Zoom level for Google maps")
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

# Make color key
colors = unique(popdata[,c("pop", "color")])
colorkey = colors$color
names(colorkey) = colors$pop

# Plot data
toplot = unique(popdata[,c("pop", "LATITUDE", "LONGITUDE")])
xlim = get_limits(toplot$LONGITUDE, tolerance=0.5)
ylim = get_limits(toplot$LATITUDE, tolerance=1)
box = make_bbox(lon=xlim, lat=ylim)
world=map_data("world")
world$lat = world$lat - 0.075 # Slight nudge b/c Stamen maps and World map data off by just a smidge
world_labels = world %>% group_by(region) %>%  summarise(long = mean(long), lat = mean(lat))
terrain = get_stamenmap(bbox = box,  maptype = "watercolor", zoom = args$zoom)

# Make plot
mymap = ggmap(terrain, extent='normal') + 
    geom_polygon(data=world, mapping=aes(x=long, y=lat, group=group), fill=NA, color="black", size=0.1) + 
    geom_text(data=world_labels, mapping=aes(x=long, y=lat, label = region), size = 3, hjust = 0.5, fontface=2) +
    geom_point(data=toplot, mapping=aes(x=LONGITUDE, y=LATITUDE, color=pop), size=2.5, alpha=0.5) +
    labs(x="Longitude", y="Latitude") +
    guides(fill=FALSE) +
    coord_fixed(1.3, xlim=xlim, ylim=ylim) +
    scale_color_manual(values=colorkey) +
    theme_bw() +
    theme(panel.background = element_rect(fill='cornflowerblue'), panel.grid.minor=element_blank())

# Save 
ggsave(mymap, file=paste(args$outprefix, ".png", sep=""), dpi=600)
ggsave(mymap, file=paste(args$outprefix, ".svg", sep=""), dpi=600)
ggsave(mymap, file=paste(args$outprefix, ".pdf", sep=""), dpi=600)

