#! /bin/bash 

# During review, some concerns came up about heterozygosity and how it was affecting our diversity analysis. Moving to a kmer-based approach to confirm.
# Jellyfish 2.2.4 was installed using Ubuntu's apt repository, `sudo apt-get install jellyfish`, because I couldn't find the right conda repository for it.

# Directories
datadir=0_data
workdir=4_Kmers
jellydir=$workdir/4a_jellyfish
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $jellydir ]; then mkdir $jellydir; fi

# Variables
kmer_size=31  # Size of kmers to use
num_kmer_markers=50000
ncores=7    # Number of parallel cores to use

##############
# CONDA ENVIRONMENT
##############

# Conda environment (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_env=proj-fonio-diversity-2020 
TASSEL5="run_pipeline.pl -Xms10g -Xmx50g"   # Shortcut TASSEL command

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate $conda_env


###############
# Kmer estimation
###############


# Use Jellyfish to count kmers in each sample
for infile in $datadir/*.fq.gz; do
    sample=`basename $infile`
    sample=${sample/\.digested*/}
    echo "Counting kmers in $sample"

    jellyfish count --mer-len $kmer_size --size 100M --threads $ncores --canonical -o $jellydir/$sample.jf --out-counter-len 1   <(zcat $infile) # 1 byte output b/c don't need to know if have >255 instances of a kmer
    jellyfish dump --column --tab -o $jellydir/$sample.tsv $jellydir/$sample.jf
    
#     break
done
# Zip up text files from Jellyfish
ls $jellydir/*.tsv | parallel gzip


# Tally up kmers and find the ones that are potentials to use (at least 2 reads in at least 2 samples), and select a random subset of ~50k to use as markers
python3 4b_SelectKmers.py -i $jellydir/*.tsv.gz --all $workdir/4b_kmers.potential.txt --subset $workdir/4b_kmers.selected.txt  -n $num_kmer_markers --min-count 2 --min-samples 2 #--debug
gzip $workdir/4b_kmers.potential.txt

# Make presence/absence genotype file
python3 4c_MakeKmerGenotypes.py -k $workdir/4b_kmers.selected.txt -i $jellydir/*.tsv.gz -o $workdir/4c_kmers.presence.txt --hapmap $workdir/4c_kmers.hmp.txt #--debug
$TASSEL5 -h $workdir/4c_kmers.hmp.txt -export $workdir/4c_kmers.vcf.gz -exportType VCF

##############
# BASIC FILTERING
##############

# Get reports
$TASSEL5 -h $workdir/4c_kmers.hmp.txt -genotypeSummary site -export $workdir/4d_sitesummary.txt
$TASSEL5 -h $workdir/4c_kmers.hmp.txt -genotypeSummary taxa -export $workdir/4d_taxasummary.txt
Rscript -e "kmer=read.csv('$workdir/4c_kmers.presence.txt', row.names=1); depth=data.frame(ID=rownames(kmer), INFO=rowSums(kmer))" \
    -e "write.table(depth, file='$workdir/4d_depth.txt', quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\\\t')"

# Plot 
Rscript 1c_PlotGenoSummary.r --sitefile $workdir/4d_sitesummary.txt --taxafile $workdir/4d_taxasummary.txt --depthfile $workdir/4d_depth.txt -o $workdir/4d_summary_plots.png


# Based on the above, use the below filters
site_max_depth=300 # About 10^2.5, which is where the depth inflection point is
site_min_maf=0.025
Rscript 1d_GetFilterLists.r --sitefile $workdir/4d_sitesummary.txt --taxafile $workdir/4d_taxasummary.txt --depthfile $workdir/4d_depth.txt \
    --site-max-depth $site_max_depth --site-min-maf $site_min_maf \
    --outtaxa $workdir/4e_taxa_to_keep.txt --outsites $workdir/4e_sites_to_keep.coords.txt --outsitenames $workdir/4e_sites_to_keep.names.txt

# Perform actual filtering
$TASSEL5 -h $workdir/4c_kmers.hmp.txt -includeSiteNamesInFile $workdir/4e_sites_to_keep.names.txt -includeTaxaInFile $workdir/4e_taxa_to_keep.txt -export $workdir/4f_kmers.filtered.vcf.gz -exportType VCF



# Try filtering for only high-confidence kmers instead of filtering them out.
# Based on the above, use the below filters
site_min_depth=250
site_min_maf=0.025
Rscript 1d_GetFilterLists.r --sitefile $workdir/4d_sitesummary.txt --taxafile $workdir/4d_taxasummary.txt --depthfile $workdir/4d_depth.txt \
    --site-min-depth $site_min_depth --site-min-maf $site_min_maf \
    --outtaxa $workdir/4g_taxa_to_keep.highcount.txt --outsites $workdir/4g_sites_to_keep.coords.highcount.txt --outsitenames $workdir/4g_sites_to_keep.names.highcount.txt
$TASSEL5 -h $workdir/4c_kmers.hmp.txt -includeSiteNamesInFile  $workdir/4g_sites_to_keep.names.highcount.txt -includeTaxaInFile $workdir/4g_taxa_to_keep.highcount.txt -export $workdir/4h_kmers.filtered.highcount.vcf.gz -exportType VCF
