#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate5_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate5_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R07-OM-EE-11
R07-OM-EE-12
R07-OM-EE-13
R07-OM-EE-14
R07-OM-EE-15
R07-OM-EE-16
R07-OM-EE-17
R07-OM-EE-18
R07-OM-EE-19
R07-OM-EE-20
R07-OM-EE-21
R07-OM-EE-22
R05-CB-ES-01
R05-CB-ES-02
R05-CB-ES-03
R05-CB-ES-04
R05-CB-ES-05
R05-CB-ES-06
R05-CB-ES-07
R05-CB-ES-08
R05-CB-ES-09
R05-CB-ES-10
R05-CB-ES-11
R05-CB-ES-12
R05-CB-ES-13
R05-CB-ES-14
R05-CB-ES-15
R05-CB-ES-16
R05-CB-ES-17
R05-CB-ES-18
R05-CB-ES-19
R05-CB-ES-20
R05-CB-ES-21
R05-CB-ES-22
R05-CB-ES-23
R05-CB-ES-24
R05-CB-ES-25
R05-CB-ES-26
R05-CB-ES-27
R05-CB-ES-28
R05-CB-ES-29
R05-CB-ES-30
R04-KB-MC-01
R04-KB-MC-02
R04-KB-MC-03
R04-KB-MC-04
R04-KB-MC-05
R04-KB-MC-06
R04-KB-MC-07
R04-KB-MC-08
R04-KB-MC-09
R04-KB-MC-10
R04-KB-MC-11
R04-KB-MC-12
R04-KB-MC-13
R04-KB-MC-14
R04-KB-MC-15
R04-KB-MC-16
R04-KB-MC-17
R04-KB-MC-18
R04-KB-MC-19
R04-KB-MC-20
R04-KB-MC-21
R04-KB-MC-22
R04-KB-MC-23
R04-KB-MC-24"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate5_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate5/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate5/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate5.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




