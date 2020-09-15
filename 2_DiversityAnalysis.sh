#! /bin/bash

# Run Fonio diversity analysis

# Software
STRUCTURE="python2 /home/jgwall/Software/Phylogeny/fastStructure/"
PLINK="/home/jgwall/Software/GWAS/plink-1.9-beta-x86_64/plink"

# Directories
datadir=0_data
snpdir=1_CallSnps
workdir=2_Diversity
if [ ! -e $workdir ]; then mkdir $workdir; fi

# Variables
snps=$snpdir/1e_genos_filtered.vcf.gz
num_threads=7


# Alternative values for doing checks requested by Reviewer #2
# # 10% heterozygosity cutoff
snps=$snpdir/1k_genos_filtered.het10.vcf.gz
workdir=2_Diversity.het10
if [ ! -e $workdir ]; then mkdir $workdir; fi

# # # 15% heterozygosity cutoff
# snps=$snpdir/1l_genos_filtered.het15.vcf.gz
# workdir=2_Diversity.het15
# if [ ! -e $workdir ]; then mkdir $workdir; fi





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
# Population structure
##############

structdir=$workdir/2a_structure
if [ ! -e $structdir ]; then mkdir $structdir; fi

# Convert SNPs to PLINK format
$TASSEL5 -vcf $snps -export $structdir/2a_snps.plink -exportType Plink
$PLINK --file $structdir/2a_snps.plink.plk --make-bed --out $structdir/2a_snps.plink.plk --allow-extra-chr

# Run fastStructure analysis for k from 1 to 10 (=1-10 subpops)
for k in `seq 1 10`; do
    echo "Processing fastStructure for k=$k"
    $STRUCTURE/structure.py -K $k --input $structdir/2a_snps.plink.plk --output $structdir/2b_structure --cv 0 --seed $k
    $STRUCTURE/distruct.py -K $k --input $structdir/2b_structure --output $structdir/2b_structure.$k.png
#     break
done
$STRUCTURE/chooseK.py --input $structdir/2b_structure | tee $structdir/2c_best_structure.txt

# Based on the above, the best structure is k=3
best_k=3
cut -f2 $structdir/2a_snps.plink.plk.ped | paste - $structdir/2b_structure.$best_k.meanQ > $workdir/2a_structure_calls.txt


##############
# Principal coordinates
##############

$TASSEL5 -vcf $snps -distanceMatrix -export $workdir/2b_distances.txt
Rscript 2b_PlotPCs.r -i $workdir/2b_distances.txt -q $workdir/2a_structure_calls.txt -o $workdir/2b_mds


##############
# Dendrogram
##############

Rscript 2c_MakeDendrograms.r -i $workdir/2b_distances.txt --popkey $workdir/2b_mds.txt -o $workdir/2c_snp_tree

##############
# Map
##############

Rscript 2d_PlotMap.r -k fonio.master_key.csv --popkey $workdir/2b_mds.txt -o $workdir/2d_map
