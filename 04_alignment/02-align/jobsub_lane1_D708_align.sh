#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D708_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D708_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-SS-FA-18
R02-SS-FA-19
R02-SS-FA-20
R02-SS-FA-21
R02-SS-FA-22
R02-SS-FA-23
R02-SS-FA-24
R02-SS-FA-25
R02-SS-FA-26
R02-SS-FA-27
R02-SS-FA-28
R02-SS-FA-29
R02-SS-FA-30
R02-SS-CR-01
R02-SS-CR-02
R02-SS-CR-03
R02-SS-CR-04
R02-SS-CR-05
R02-SS-CR-06
R02-SS-CR-07
R02-SS-CR-08
R02-SS-CR-09
R02-SS-CR-10
R02-SS-CR-11
R02-SS-CR-12
R02-SS-CR-13
R02-SS-CR-14
R02-SS-CR-15
R02-SS-CR-16
R02-SS-CR-17
R02-SS-CR-18
R02-SS-CR-19
R02-SS-CR-20
R02-SS-CR-21
R02-SS-CR-22
R02-SS-CR-23
R02-SS-CR-24
R02-SS-CR-25
R02-SS-CR-26
R02-SS-CR-27
R02-SS-CR-28
R02-SS-CR-29
R02-SS-CR-30
R01-SI-MO-01
R01-SI-MO-02
R01-SI-MO-03
R01-SI-MO-04
R01-SI-MO-05
R01-SI-MO-06
R01-SI-MO-07
R01-SI-MO-08
R01-SI-MO-09
R01-SI-MO-10
R01-SI-MO-11
R01-SI-MO-12
R01-SI-MO-13
R01-SI-MO-14
R01-SI-MO-15
R01-SI-MO-16
R01-SI-MO-17
R01-SI-MO-18
R01-SI-MO-19
R01-SI-MO-20
R01-SI-MO-21
R01-SI-MO-22
R01-SI-MO-23
R01-SI-MO-24
R01-SI-MO-25
R01-SI-MO-26
R01-SI-MO-27
R01-SI-MO-28
R01-SI-MO-29
R01-SI-MO-30
R01-SI-FR-01
R01-SI-FR-02
R01-SI-FR-03
R01-SI-FR-04
R01-SI-FR-05
R01-SI-FR-06
R01-SI-FR-07
R01-SI-FR-08
R01-SI-FR-09
R01-SI-FR-10
R01-SI-FR-11
R01-SI-FR-12
R01-SI-FR-13
R01-SI-FR-14
R01-SI-FR-15
R01-SI-FR-16
R01-SI-FR-17
R01-SI-FR-18
R01-SI-FR-19
R01-SI-FR-20
R01-SI-FR-21
R01-SI-FR-22
R01-SI-FR-23"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D708_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D708_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D708.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




