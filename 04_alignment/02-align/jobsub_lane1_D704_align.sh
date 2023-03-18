#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D704_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D704_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-SS-LT-19
R02-SS-LT-20
R02-SS-LT-21
R02-SS-LT-22
R02-SS-LT-23
R02-SS-LT-24
R02-SS-LT-25
R02-SS-LT-26
R02-SS-LT-27
R02-SS-LT-28
R02-SS-LT-29
R02-SS-LT-30
R06-GI-MY-11
R06-GI-MY-12
R06-GI-MY-13
R06-GI-MY-14
R06-GI-MY-15
R06-GI-MY-16
R06-GI-MY-17
R06-GI-MY-18
R06-GI-MY-19
R06-GI-MY-20
R06-GI-MY-21
R06-GI-MY-22
R06-GI-MY-23
R06-GI-MY-24
R06-GI-MY-25
R06-GI-MY-26
R06-GI-MY-27
R06-GI-MY-28
R06-GI-MY-29
R06-GI-MY-30
R01-SI-LL-11
R01-SI-LL-12
R01-SI-LL-13
R01-SI-LL-14
R01-SI-LL-15
R01-SI-LL-16
R01-SI-LL-17
R01-SI-LL-18
R01-SI-LL-19
R01-SI-LL-20
R01-SI-LL-21
R01-SI-LL-22
R01-SI-LL-23
R01-SI-LL-24
R01-SI-LL-25
R01-SI-LL-26
R01-SI-LL-27
R01-SI-LL-28
R01-SI-LL-29
R01-SI-LL-30
R01-CR-CE-11
R01-CR-CE-12
R01-CR-CE-13
R01-CR-CE-14
R01-CR-CE-15
R01-CR-CE-16
R01-CR-CE-17
R01-CR-CE-18
R01-CR-CE-19
R01-CR-CE-20
R01-CR-CE-21
R01-CR-CE-22
R01-CR-CE-23
R01-CR-CE-24
R01-CR-CE-25
R01-CR-CE-26
R01-CR-CE-27
R01-CR-CE-28
R01-CR-CE-29
R01-CR-CE-30
R01-SI-CA-11
R01-SI-CA-12
R01-SI-CA-13
R01-SI-CA-14
R01-SI-CA-15
R01-SI-CA-16
R01-SI-CA-17
R01-SI-CA-18
R01-SI-CA-19
R01-SI-CA-20
R01-SI-CA-21
R01-SI-CA-22
R01-SI-CA-23
R01-SI-CA-24
R01-SI-CA-25
R01-SI-CA-26
R01-SI-CA-27
R01-SI-CA-28
R01-SI-CA-29
R01-SI-CA-30
R02-CW-MI-11
R02-CW-MI-12
R02-CW-MI-13
R02-CW-MI-14"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D704_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D704_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D704.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




