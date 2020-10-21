#! /bin/bash

# Call SNPs on the fonio genome assembly using the raw data from Data2Bio

# Variables
genome=/home/jgwall/Projects/0_RawData/FonioGenome/FN_CANU_PILON.v1.0.fasta.gz
max_depth=1000  # Max depth recorded when calling SNPs
min_base_quality=20    # Min base quality when SNP calling
num_threads=7

# Directories
datadir=0_data
workdir=1_CallSnps.max_1_path
aligndir=$workdir/1a_alignment
snpdir=$workdir/1b_snps
if [ ! -e $workdir ]; then mkdir $workdir; fi
if [ ! -e $aligndir ]; then mkdir $aligndir; fi
if [ ! -e $snpdir ]; then mkdir $snpdir; fi

##############
# CONDA ENVIRONMENT
##############

# Conda environment (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_env=proj-fonio-diversity-2020 
TASSEL5="run_pipeline.pl -Xms10g -Xmx50g"   # Shortcut TASSEL command

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate $conda_env

##############
# CALL SNPS
##############


# # Build genome database for GSNAP
# gmap_build --dir $workdir --db fonio_v1 --gunzip $genome

# # Align reads with GSNAP
# for infile in $datadir/*.fq.gz; do
#     sample=`basename $infile`
#     sample=${sample/\.digested*/}
#     echo "Aligning $sample"

#     # Align with max 1 path allowed
#     gsnap --dir $workdir --db fonio_v1 --gunzip --nthreads $num_threads --format sam  $infile --npaths 1 --quiet-if-excessive | samtools view -S -h -b - > $aligndir/$sample.tmp.max_1_path.bam
    
    
#     # Sort and index with max 1 path
#     samtools sort $aligndir/$sample.tmp.max_1_path.bam $aligndir/$sample.max_1_path
#     samtools index $aligndir/$sample.max_1_path.bam
#     rm $aligndir/$sample.tmp.max_1_path.bam

    #     break
# done

# # Index genome for samtools
# samtools faidx $genome  
# 
# # Call SNPs on all samples; do in parallel so easier to see progress. (Also --threads option doesn't seem as efficient as it should be)
# # Note: Chained commands to keep things from being crazy large, since mpileup by default lists any base with any reads at it
# commands=$workdir/1b_call_commands.txt
# echo "echo 'Calling SNPs in parallel'" > $commands
# for chrom in `cut -f1 $genome.fai`; do
#     sample=${chrom/|*/}
# 
#     echo "bcftools mpileup -d $max_depth -f $genome -Q $min_base_quality --output-type u $aligndir/*.bam -r '$chrom' | \
#     bcftools call --multiallelic-caller --output-type u | \
#     bcftools view --min-alleles 2 --output-type b --output-file $snpdir/1c_snps.$sample.bcf" >> $commands
# 
# done
# cat $commands | parallel --progress

  
# # Combine all called SNPs into one file, excluding other polymorphism types
# bcftools concat --threads $num_threads --output-type u  $snpdir/*.bcf | bcftools view --types snps --output-type z --output-file $workdir/1b_snps_combined.vcf.gz 



##############
# BASIC FILTERING
##############

# Get reports
$TASSEL5 -vcf $workdir/1b_snps_combined.vcf.gz -genotypeSummary site -export $workdir/1c_sitesummary.txt
$TASSEL5 -vcf $workdir/1b_snps_combined.vcf.gz -genotypeSummary taxa -export $workdir/1c_taxasummary.txt
zcat $workdir/1b_snps_combined.vcf.gz | cut -f3,8 | grep -v "^#" | sed -r -e "s|DP=([0-9]+).+|\1|" > $workdir/1c_depth.txt

# Plot 
Rscript 1c_PlotGenoSummary.r --sitefile $workdir/1c_sitesummary.txt --taxafile $workdir/1c_taxasummary.txt --depthfile $workdir/1c_depth.txt -o $workdir/1c_summary_plots.png


# Stopped here because SNP distributions were not very different from the default pipeline
