#!/bin/bash
#SBATCH --mail-user=hector.rouxdebezieux@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
loc="/scratch/users/singlecell/Pancreas/ProcessedData/"
out="/accounts/projects/epurdom/singlecell/Dune_Paper/Mapping/Scripts/Pancreas.csv"
Rscript --verbose Main.R -l $loc --name1="baron" --name2="segerstolpe"\
     -m "~/Pancreas/Data/Dune/" -o $out > pancreas.out 2>&1
