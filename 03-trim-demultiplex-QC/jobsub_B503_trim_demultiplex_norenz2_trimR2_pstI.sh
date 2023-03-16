#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=1000
#SBATCH --job-name=jobsub_B503_trim_demultiplex_norenz2_trimR2_pstI
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.cpom
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_B503_trim_demultiplex_norenz2_trimR2_pstI_%j.out


# rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/jobsub_B503_trim_demultiplex_norenz2_trimR2_pstI.sh /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/stacks_barcode_B503.txt roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/ --progress



module load StdEnv/2020
module load stacks/2.60


echo “Starting demultiplex,clean at: `date`”

process_radtags -1 /project/def-saitken/roseanna/rawreads_april152022/NS.1760.001.B709---B503.Hamelin_202111_pool1_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.1760.001.B709---B503.Hamelin_202111_pool1_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_B503.txt -o /scratch/roseanna/Demultiplexing_stacks/B503_norenz2_trimR2 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 pstI --quality --rescue --barcode_dist_1 1 -D &> process_radtags_standoutputerror_B503_norenz2_trimR2_pstI.oe 

echo “Job finished with exit code $? at: `date`”


