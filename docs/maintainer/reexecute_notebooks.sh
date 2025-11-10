#!/bin/bash
#SBATCH -A IscrC_aGAL1
#SBATCH -p boost_usr_prod
#SBATCH --time=05:00:00
#SBATCH --job-name=nb_run
#SBATCH --output=nb_run_%j.out
#SBATCH --error=nb_run_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Run a notebook (with outputs saved into the same .ipynb)
for f in estimators projections rmsd NTL9_plots villin_plots villin_data ; do 
  uv run jupyter  nbconvert --execute --inplace ../$f.ipynb &
done

wait
