#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=1000
#SBATCH --job-name=jobsub_lane3_plate1_trim_demultiplex_norenz2_trimR2
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate1_trim_demultiplex_norenz2_trimR2_%j.out

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/Trim_demultiplex_lane3/jobsub_lane3_plate1_trim_demultiplex_norenz2_trimR2.sh /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/Trim_demultiplex_lane3/stacks_barcode_lane3_plate1.txt roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/ --progress


module load StdEnv/2020
module load trimmomatic/0.39
module load stacks/2.62

cpu=1



echo “Starting run at: `date`”

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 1 -trimlog /scratch/roseanna/Trimmed_reverseplates/lane3_plate1_R2_trimmed_log.txt /home/roseanna/scratch/rawreads_lane3_2023/NS.2053.003.B715---D502.Hamelin__20230109-Plate-1_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/NS.2053.003.B715---D502.Hamelin__20230109-Plate-1_R2.trimmed.fastq.gz HEADCROP:3


echo “Starting demultiplex,clean at: `date`”

process_radtags -1 /home/roseanna/scratch/rawreads_lane3_2023/NS.2053.003.B715---D502.Hamelin__20230109-Plate-1_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.2053.003.B715---D502.Hamelin__20230109-Plate-1_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_lane3_plate1.txt -o /scratch/roseanna/Demultiplexing_stacks/lane3_plate1 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 pstI --quality --rescue --barcode-dist-1 1 --threads $cpu -D &> process_radtags_standoutputerror_lane3_plate1.oe 

echo “Job finished with exit code $? at: `date`”