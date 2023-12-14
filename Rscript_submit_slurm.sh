#!/bin/bash
#
#SBATCH --job-name=Rscript_submit_himem_1node_12cpus_mem500_231125
#SBATCH --output=/home/fmbuga/tools/slurm_scripts/slurm_out_err/Rscript_submit_himem_1node_12cpus_mem500_231125_%j.out
#SBATCH --error=/home/fmbuga/tools/slurm_scripts/slurm_out_err/Rscript_submit_himem_1node_12cpus_mem500_231125_%j.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=500G
#SBATCH --partition=himem


datenow=$(date)
echo $datenow
srun hostname

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME


eval "$(conda shell.bash hook)"
conda activate seurat5
echo ""
conda info
echo ""
conda list
echo ""

echo ""
echo "Path to R used: "
which R
R --version


###

killall R

Rscript \
    --save \
    --verbose \
    $1 \
    $2 \
    $3 \
    $4

###

echo ""
srun hostname


echo ""
end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"


