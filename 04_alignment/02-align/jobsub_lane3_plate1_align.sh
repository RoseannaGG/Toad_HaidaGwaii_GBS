#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate1_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate1_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R03-TO-ML-01
R03-TO-ML-02
R03-TO-ML-03
R03-TO-ML-04
R03-TO-ML-05
R03-TO-ML-06
R03-TO-ML-07
R03-TO-ML-08
R03-TO-ML-09
R03-TO-ML-10
R03-TO-ML-11
R03-TO-ML-12
R03-TO-ML-13
R03-TO-ML-14
R03-TO-ML-15
R03-TO-ML-16
R03-TO-ML-17
R03-TO-ML-18
R03-TO-ML-19
R03-TO-ML-20
R03-TO-ML-21
R03-TO-ML-22
R03-TO-ML-23
R03-TO-ML-24
R03-TO-ML-25
R03-TO-ML-26
R03-TO-ML-27
R03-TO-ML-28
R03-TO-ML-29
R03-TO-ML-30
R03-TO-KA-01
R03-TO-KA-02
R03-TO-KA-03
R03-TO-KA-04
R03-TO-KA-05
R03-TO-KA-06
R03-TO-KA-07
R03-TO-KA-08
R03-TO-KA-09
R03-TO-KA-10
R03-TO-KA-11
R03-TO-KA-12
R03-TO-KA-13
R03-TO-KA-14
R03-TO-KA-15
R03-TO-KA-16
R03-TO-KA-17
R03-TO-KA-18
R03-TO-KA-19
R03-TO-KA-20
R03-TO-KA-21
R03-TO-KA-22
R03-TO-KA-23
R03-TO-KA-24
R03-TO-KA-25
R03-TO-KA-26
R03-TO-KA-27
R03-TO-KA-28
R03-TO-KA-29
R03-TO-KA-30
R04-KB-SU-01
R04-KB-SU-02
R04-KB-SU-03
R04-KB-SU-04
R04-KB-SU-05
R04-KB-SU-06
R04-KB-SU-07
R04-KB-SU-08
R04-KB-SU-09
R04-KB-SU-10
R04-KB-SU-11
R04-KB-SU-12
R04-KB-SU-13
R04-KB-SU-14
R04-KB-SU-15
R04-KB-SU-16
R04-KB-SU-17
R04-KB-SU-18
R04-KB-SU-19
R04-KB-SU-20
R04-KB-SU-21
R04-KB-SU-22
R04-KB-SU-23
R04-KB-SU-24
R04-KB-SU-25
R04-KB-SU-26
R04-KB-SU-27
R04-KB-SU-28
R04-KB-SU-29
R04-KB-SU-30
R04-KB-FI-01
R04-KB-FI-02
R04-KB-FI-03
R04-KB-FI-04
R04-KB-FI-05
R04-KB-FI-06"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate1_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate1/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate1/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate1.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




