#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate6_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate6_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R04-KB-MC-25
R04-KB-MC-26
R04-KB-MC-27
R04-KB-MC-28
R04-KB-MC-29
R04-KB-MC-30
R04-KB-CO-01
R04-KB-CO-02
R04-KB-CO-03
R04-KB-CO-04
R04-KB-CO-05
R04-KB-CO-06
R04-KB-CO-07
R04-KB-CO-08
R04-KB-CO-09
R04-KB-CO-10
R04-KB-CO-11
R04-KB-CO-12
R04-KB-CO-13
R04-KB-CO-14
R04-KB-CO-15
R04-KB-CO-16
R04-KB-CO-17
R04-KB-CO-18
R04-KB-CO-19
R04-KB-CO-20
R04-KB-CO-21
R04-KB-CO-22
R04-KB-CO-23
R04-KB-CO-24
R04-KB-CO-25
R04-KB-CO-26
R04-KB-CO-27
R04-KB-CO-28
R04-KB-CO-29
R04-KB-CO-30
R04-KB-GO-01
R04-KB-GO-02
R04-KB-GO-03
R04-KB-GO-04
R04-KB-GO-05
R04-KB-GO-06
R04-KB-GO-07
R04-KB-GO-08
R04-KB-GO-09
R04-KB-GO-10
R04-KB-GO-11
R04-KB-GO-12
R04-KB-GO-13
R04-KB-GO-14
R04-KB-GO-15
R04-KB-GO-16
R04-KB-GO-17
R04-KB-GO-18
R04-KB-GO-19
R04-KB-GO-20
R04-KB-GO-21
R04-KB-GO-22
R04-KB-GO-23
R04-KB-GO-24
R04-KB-GO-25
R04-KB-GO-26
R04-KB-GO-27
R04-KB-GO-28
R04-KB-GO-29
R04-KB-GO-30"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate* roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate6/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate6/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate6.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




