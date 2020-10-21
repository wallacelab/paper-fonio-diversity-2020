#! /usr/bin/env bash 

# Create the Conda environment for this analysis; only has to be done at the beginning

# Environment name
conda_env=proj-fonio-diversity-2020 

#################
# Original Environment creation commands - DO NOT USE
#################

# # Create environment (ORIGINAL COMMANDS - DO NOT USE)
# conda create -n $conda_env bcftools=1.9 vcftools=0.1.16 tassel=5.2.40 r-base=3.6.2 r-essentials=3.6 pandas=1.0.3 matplotlib=3.2.1 sourmash=3.5.0

# # Export environment for others to load (needs to be done with the conda environment active)
# conda list --explicit > conda_environment.txt


#################
# Create a new environment using these commands
#################

# # Create environment with explicit specification
# conda create --name $conda_env --file conda_environment.txt


#################
# Loading environment in scripts
#################

# To be able to load the conda environment within scripts, you need to include the following line of code:
. $(conda info --root)/etc/profile.d/conda.sh 

# Activate said environment
conda activate $conda_env






