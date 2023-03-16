#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_D504_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_D504_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R01-SI-HR-13
R01-SI-HR-14
R01-SI-HR-15
R01-SI-HR-16
R01-SI-HR-17
R01-SI-HR-18
R01-SI-HR-19
R01-SI-HR-20
R01-SI-HR-21
R01-SI-HR-22
R01-SI-HR-23
R01-SI-HR-24
R01-SI-HR-25
R01-SI-HR-26
R01-SI-HR-27
R01-SI-HR-28
R01-SI-HR-29
R01-SI-HR-30
R02-CW-AG-01
R02-CW-AG-02
R02-CW-AG-03
R02-CW-AG-04
R02-CW-AG-05
R02-CW-AG-06
R02-CW-AG-07
R02-CW-AG-08
R02-CW-AG-09
R02-CW-AG-10
R02-CW-AG-11
R02-CW-AG-12
R02-CW-AG-13
R02-CW-AG-14
R02-CW-AG-15
R02-CW-AG-16
R02-CW-AG-17
R02-CW-AG-18
R02-CW-AG-19
R02-CW-AG-20
R02-CW-AG-21
R02-CW-AG-22
R02-CW-AG-23
R02-CW-AG-24
R02-CW-AG-25
R02-CW-AG-26
R02-CW-AG-27
R02-CW-AG-28
R02-CW-AG-29
R02-CW-AG-30
R02-CW-RY-01
R02-CW-RY-02
R02-CW-RY-03
R02-CW-RY-04
R02-CW-RY-05
R02-CW-RY-06
R02-CW-RY-07
R02-CW-RY-08
R02-CW-RY-09
R02-CW-RY-10
R02-CW-RY-11
R02-CW-RY-12
R02-CW-RY-13
R02-CW-RY-14
R02-CW-RY-15
R02-CW-RY-16
R02-CW-RY-17
R02-CW-RY-18
R02-CW-RY-19
R02-CW-RY-20
R02-CW-RY-21
R02-CW-RY-22
R02-CW-RY-23
R02-CW-RY-24
R02-CW-RY-25
R02-CW-RY-26
R02-CW-RY-27
R02-CW-RY-28
R02-CW-RY-29
R02-CW-DE-01
R02-CW-DE-02
R02-CW-DE-03
R02-CW-DE-04
R02-CW-DE-05
R02-CW-DE-06
R02-CW-DE-07
R02-CW-DE-08
R02-CW-DE-09
R02-CW-DE-10
R02-CW-DE-11
R02-CW-DE-12
R02-CW-DE-13
R02-CW-DE-14
R02-CW-DE-15
R02-CW-DE-16
R02-CW-DE-17
R02-CW-DE-18
R02-CW-DE-19"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D504_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D504_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_D504.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




