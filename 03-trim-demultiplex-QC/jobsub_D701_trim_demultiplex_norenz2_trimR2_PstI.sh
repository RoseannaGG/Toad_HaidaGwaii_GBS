#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=1000
#SBATCH --job-name=jobsub_D701_trim_demultiplex_norenz2_trimR2_PstI
#SBATCH --mail-user=roseanna.gamlen.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_D701_trim_demultiplex_norenz2_trimR2_PstI_%j.out


#rsync -zv /drives/f/GBS_data_03_02_21/Demultiplexing_stacks/stacks_barcode_D701.txt roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/ --progress

#rsync -zv /drives/f/GBS_data_03_02_21/Trim_demultiplex_Oct2021/jobsub_D701_trim_demultiplex_norenz2_trimR2_PstI.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/ --progress

module load StdEnv/2020
module load trimmomatic/0.39
module load stacks/2.60


echo “Starting run at: `date`”

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 4 -trimlog /scratch/roseanna/Trimmed_reverseplates/D701_P1_R2_trimmed_log.txt /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D701.Hamelin_202009_P1_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D701.Hamelin_202009_P1_R2.trimmed.fastq.gz HEADCROP:3


echo “Starting demultiplex,clean at: `date`”

process_radtags -1 /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D701.Hamelin_202009_P1_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D701.Hamelin_202009_P1_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_D701.txt -o /scratch/roseanna/Demultiplexing_stacks/D701_norenz2_trimR2_PstI -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 pstI --quality --rescue --barcode_dist_1 1 -D &> process_radtags_standoutputerror_D701_norenz2_trimR2_PstI.oe 

echo “Job finished with exit code $? at: `date`”

