#! /usr/bin/Rscript

# Apply filters to determine which sites/taxa to keep

# Libraries
library(argparse)
library(seqinr)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="TASSEL genotype summary file of sites")
parser$add_argument("-f", "--fastafile", help="FASTA-formatted genome to look up SNP context")
parser$add_argument("-o", "--outfile", help="Output file")
parser$add_argument("-p", "--pad", type='integer', default=50, help="How far to go from the SNP to get genomic context")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/MinorMillets/ICRISAT_Fonio/2020_04_DiversityAnalysis/1_CallSnps/')
# args=parser$parse_args(c("-i","1h_sitesummary.pretty.txt", "-f", "~/Projects/0_RawData/FonioGenome/FN_CANU_PILON.v1.0.fasta.gz", "-o", "99_tmp.txt"))

# Load data
cat("Making SNP summary table\n")
sites = read.delim(args$infile, check.names=F)
genome = read.fasta(args$fastafile, as.string=T, forceDNAtolower = F, seqonly=F)

# Get basic SNP summary data
snpdata = sites[,c("Site Name", "Chromosome", "Physical Position", "Major Allele", "Minor Allele", "Minor Allele Frequency")]
names(snpdata) = c("Name", "Chromosome", "Position", "MajorAllele", "MinorAllele", "MinorAlleleFrequency")

# Correct capitalization TASSEL forces onto chromsomes
snpdata$Chromosome = sub(tolower(snpdata$Chromosome), pattern='version', repl='Version')


# Get genome context
context = lapply(1:nrow(snpdata), function(i){
    
    # Get start and stop coordinates
    mypos = snpdata$Position[i]
    start = mypos - args$pad
    if(start < 1){
        start=1
    }
    stop = mypos + args$pad
    
    # Extract flanking sequence
    mychrom = snpdata$Chromosome[i]
    left = substr(genome[[mychrom]], start=start, stop=mypos-1)
    right = substr(genome[[mychrom]], start=mypos+1, stop=stop)
    
    # Make context
    return(data.frame(left=as.character(left), right=as.character(right)))
})
context = do.call(rbind, context)

# Add context to summary
snpdata$Context = paste(context$left, "[", snpdata$MajorAllele, "/", snpdata$MinorAllele, "]", context$right, sep="")


# Write out
write.table(snpdata, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)
