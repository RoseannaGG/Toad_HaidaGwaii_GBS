#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16000
#SBATCH --job-name=jobsub_index_BWA_Anaxyrusboreas
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_index_BWA_Anaxyrusboreas_%j.out

module load StdEnv/2020
module load bwa/0.7.17


echo “Starting run at: `date`”

genome_fa=/scratch/roseanna/Anaxyrus_boreas_genome/final_assembly.fasta.gz

bwa index -p /scratch/roseanna/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas $genome_fa &> /scratch/roseanna/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas_genome_index.oe

echo “Job finished with exit code $? at: `date`”