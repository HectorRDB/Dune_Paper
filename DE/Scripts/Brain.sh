#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
loc="/scratch/users/singlecell/MiniAtlas/data/rds/"
out="/accounts/projects/epurdom/singlecell/Dune_Paper/DE/Data/Brain.csv"
m1="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/singleTree/"
m2="/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/Dune/"
Rscript --verbose Main.R -l $loc --f "SMARTer_nuclei_MOp" -s "SMARTer_cells_MOp" \
     -m $m1 -n $m2 -o $out > brain.out 2>&1
