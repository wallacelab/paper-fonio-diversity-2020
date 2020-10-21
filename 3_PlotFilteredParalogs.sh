#! /bin/bash

# Plot the SNPs filtered out by different levels of heterozygosity

# Directories
snpdir=1_CallSnps
divdir=2_Diversity
workdir=3_Paralogs
if [ ! -e $workdir ]; then mkdir $workdir; fi

# Variables
snp10=$snpdir/1k_genos_filtered.het10.vcf.gz    # 10% het cutoff
snp15=$snpdir/1l_genos_filtered.het15.vcf.gz    # 15% het cutoff
snp25=$snpdir/1e_genos_filtered.vcf.gz          # 25% het cutoff (default)
pop_structure=$divdir/2b_mds.txt    # MDS results that include population structure assignment
num_threads=7



##############
# CONDA ENVIRONMENT
##############

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_env=proj-fonio-diversity-2020 
TASSEL5="run_pipeline.pl -Xms10g -Xmx50g"   # Shortcut TASSEL command

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED

conda activate $conda_env

##############
# Plot SNPs
##############

# # Convert SNPs to hapmaps
# $TASSEL5 -vcf $snp10 -export $workdir/3a_snps10.hmp.txt
# $TASSEL5 -vcf $snp15 -export $workdir/3a_snps15.hmp.txt
# $TASSEL5 -vcf $snp25 -export $workdir/3a_snps25.hmp.txt

# # Plot paralogs
# Rscript 3b_PlotSnpParalogs.r -a $workdir/3a_snps25.hmp.txt -b $workdir/3a_snps10.hmp.txt --pop-structure $pop_structure --outprefix $workdir/3b_compare_25_10
# Rscript 3b_PlotSnpParalogs.r -a $workdir/3a_snps25.hmp.txt -b $workdir/3a_snps15.hmp.txt --pop-structure $pop_structure --outprefix $workdir/3b_compare_25_15
# Rscript 3b_PlotSnpParalogs.r -a $workdir/3a_snps15.hmp.txt -b $workdir/3a_snps10.hmp.txt --pop-structure $pop_structure --outprefix $workdir/3b_compare_15_10


# I thought about trying a random subset to check, but that one will also include the ones that separate, so I'll still get clusters. Not sure the best way to check that.
