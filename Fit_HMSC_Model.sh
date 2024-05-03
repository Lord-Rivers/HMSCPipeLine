#!/bin/bash
#SBATCH --job-name=Tensorflow
#SBATCH  -M ukko
#SBATCH --partition=gpu 
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=2
#SBATCH --time=00:30:00
#SBATCH --array=0-3
#SBATCH --output=R-%x.%j.out

module purge
module load Python cuDNN
source $HOME/myenv/bin/activate

model_name=$NAME

Area=$AREA
samples=$SAMP
thin=$THIN
chains=4
verbose=100
divide=$samples*$thin; by=2; (( transient=(divide+by-1)/by ))
data_path=$(printf "%s/%s" $Area $model_name)
screen_text=$(printf "#######\n Model Run starting\n Samples: %.4d\n Thining: %.2d Transient steps: %.d\n######\n" $samples $thin $transient)

input_path=$data_path/$(printf "INIT/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds" $samples $thin $chains)
output_path=$data_path/$(printf "Sampled/HPC_samples_%.4d_thin_%.2d_chain_%.1d.rds" $samples $thin $SLURM_ARRAY_TASK_ID)

echo $screen_text
echo $data_path
echo $input_path
echo $output_path

srun python Gibs_sampling_tensor_flow.py --samples $samples --transient $transient --thin $thin --verbose $verbose --input $input_path --output $output_path --chain $SLURM_ARRAY_TASK_ID