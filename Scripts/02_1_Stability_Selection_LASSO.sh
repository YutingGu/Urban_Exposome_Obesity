#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -N dict

cd /rds/general/project/hda_23-24/live/CompEpi/group1/Scripts

module load anaconda3/personal
source activate r413

Rscript 02_stability_selection_LASSO.R
