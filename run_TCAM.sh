#!/bin/bash

# source conda
source ~/opt/miniconda3/etc/profile.d/conda.sh

# Activate spyder environment, which has the necessary packages
conda activate spyder-env



# Navigate to TCAM outputs directory
cd Outputs/beta_diversity/longitudinal/tcam/

# run python script.

# Process all data types.

python3 ../../../../code/TCAM.py  -o './taxa' -t 'taxa' -i 'taxa/table_dfb_meta.tsv'
python3 ../../../../code/TCAM.py  -o './koe' -t 'koe' -i 'koe/table_dfb_meta.tsv'
python3 ../../../../code/TCAM.py  -o './vfs' -t 'vfs' -i 'vfs/table_dfb_meta.tsv'
python3 ../../../../code/TCAM.py  -o './amr' -t 'amr' -i 'amr/table_dfb_meta.tsv'
python3 ../../../../code/TCAM.py  -o './caz' -t 'caz' -i 'caz/table_dfb_meta.tsv'
