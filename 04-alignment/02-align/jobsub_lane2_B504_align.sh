#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_B504_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_B504_align_%j.out

# rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/jobsub_lane2_B504_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/  --progress

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R07-NE-AS-20
R07-NE-AS-21
R07-NE-AS-22
R07-NE-AS-23
R07-NE-AS-24
R07-NE-AS-25
R07-NE-AS-26
R07-NE-AS-27
R07-NE-AS-28
R07-NE-AS-29
R07-NE-AS-30
R07-NE-GS-01
R07-NE-GS-02
R07-NE-GS-03
R07-NE-GS-04
R07-NE-GS-05
R07-NE-GS-06
R07-NE-GS-07
R07-NE-GS-08
R07-NE-GS-09
R07-NE-GS-10
R07-NE-GS-11
R07-NE-GS-12
R07-NE-GS-13
R07-NE-GS-14
R07-NE-GS-15
R07-NE-GS-16
R07-NE-GS-17
R07-NE-GS-18
R07-NE-GS-19
R07-NE-GS-20
R07-NE-GS-21
R07-NE-GS-22
R07-NE-GS-23
R07-NE-GS-24
R07-NE-GS-25
R07-NE-GS-26
R07-NE-GS-27
R07-NE-GS-28
R07-NE-GS-29
R07-NE-GS-30
R05-CB-DR-01
R05-CB-DR-02
R05-CB-DR-03
R05-CB-DR-04
R05-CB-DR-05
R05-CB-DR-06
R05-CB-DR-07
R05-CB-DR-08
R05-CB-DR-09
R05-CB-DR-10
R05-CB-DR-11
R05-CB-DR-12
R05-CB-DR-13
R05-CB-DR-14
R05-CB-DR-15"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/B504_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/B504_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_B504.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




