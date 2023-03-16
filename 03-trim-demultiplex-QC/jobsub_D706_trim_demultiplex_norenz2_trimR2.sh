#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=1000
#SBATCH --job-name=jobsub_D706_trim_demultiplex_norenz2_trimR2
#SBATCH --mail-user=roseanna.gamlen.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_D706_trim_demultiplex_norenz2_trimR2_%j.out

module load StdEnv/2020
module load trimmomatic/0.39
module load stacks/2.60

echo “Starting trim at: `date`”

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 4 -trimlog /scratch/roseanna/Trimmed_reverseplates/D706_P6_R2_trimmed_log.txt /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D706.Hamelin_202010_plate__2_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D706.Hamelin_202010_plate__2_R2.trimmed.fastq.gz HEADCROP:3


echo “Starting demultiplex,clean at: `date`”


process_radtags -1 /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D706.Hamelin_202010_plate__2_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D706.Hamelin_202010_plate__2_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_D706.txt -o /scratch/roseanna/Demultiplexing_stacks/D706_norenz2_trimR2 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 sbfI --quality --rescue --barcode_dist_1 1 -D &> process_radtags_standoutputerror_D706_norenz2_trimR2.oe 

echo “Job finished with exit code $? at: `date`”