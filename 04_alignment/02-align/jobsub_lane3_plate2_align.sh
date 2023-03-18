#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate2_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate2_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R04-KB-FI-07
R04-KB-FI-08
R04-KB-FI-09
R04-KB-FI-10
R04-KB-FI-11
R04-KB-FI-12
R04-KB-FI-13
R04-KB-FI-14
R04-KB-FI-15
R04-KB-FI-16
R04-KB-FI-17
R04-KB-FI-18
R04-KB-FI-19
R04-KB-FI-20
R04-KB-FI-21
R04-KB-FI-22
R04-KB-FI-23
R04-KB-FI-24
R04-KB-FI-25
R04-KB-FI-26
R04-KB-FI-27
R04-KB-FI-28
R04-KB-FI-29
R04-KB-FI-30
R05-CB-MR-01
R05-CB-MR-02
R05-CB-MR-03
R05-CB-MR-04
R05-CB-MR-05
R05-CB-MR-06
R05-CB-MR-07
R05-CB-MR-08
R05-CB-MR-09
R05-CB-MR-10
R05-CB-MR-11
R05-CB-MR-12
R05-CB-MR-13
R05-CB-MR-14
R05-CB-MR-15
R05-CB-MR-16
R05-CB-MR-17
R05-CB-MR-18
R05-CB-MR-19
R05-CB-MR-20
R05-CB-MR-21
R05-CB-MR-22
R05-CB-MR-23
R05-CB-MR-24
R05-CB-MR-25
R05-CB-MR-26
R05-CB-MR-27
R05-CB-MR-28
R05-CB-MR-29
R05-CB-DU-01
R05-CB-DU-02
R05-CB-DU-03
R05-CB-DU-04
R05-CB-DU-05
R05-CB-DU-06
R05-CB-DU-07
R05-CB-DU-08
R05-CB-DU-09
R05-CB-DU-10
R05-CB-DU-11
R05-CB-DU-12
R05-CB-DU-13
R05-CB-DU-14
R05-CB-DU-15
R05-CB-DU-16
R05-CB-DU-17
R05-CB-DU-18
R05-CB-DU-19
R05-CB-DU-20
R05-CB-DU-21
R05-CB-DU-22
R05-CB-DU-23
R05-CB-DU-24
R05-CB-DU-25
R05-CB-DU-26
R05-CB-DU-27
R05-CB-DU-28
R05-CB-DU-29
R05-CB-DU-30
R05-CB-HU-01
R05-CB-HU-02
R05-CB-HU-03
R05-CB-HU-04
R05-CB-HU-05
R05-CB-HU-06
R05-CB-HU-07
R05-CB-HU-08
R05-CB-HU-09
R05-CB-HU-10
R05-CB-HU-11
R05-CB-HU-12
R05-CB-HU-13"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate2_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate2/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate2.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




