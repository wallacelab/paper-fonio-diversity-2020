#! /usr/bin/env bash 

# Create the Conda environment for this analysis; only has to be done at the beginning

# Environment name
conda_env=proj-fonio-diversity-2020 

# # Create environment TODO: confirm this works; argparse acting up previously
conda create -n $conda_env bcftools=1.9 vcftools=0.1.16 tassel=5.2.40 r-base=3.6.2 r-essentials=3.6

# 
# # To be able to load the conda environment within scripts, you need to include the following line of code:
# . $(conda info --root)/etc/profile.d/conda.sh 
# 
# # Activate said environment
# conda activate $conda_env
