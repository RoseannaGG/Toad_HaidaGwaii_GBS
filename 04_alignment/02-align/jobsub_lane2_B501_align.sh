#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_B501_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_B501_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-CW-DE-20
R02-CW-DE-21
R02-CW-DE-22
R02-CW-DE-23
R02-CW-DE-24
R02-CW-DE-25
R02-CW-DE-26
R02-CW-DE-27
R02-CW-DE-28
R02-CW-DE-29
R02-CW-DE-30
R02-CW-SI-01
R02-CW-SI-02
R02-CW-SI-03
R02-CW-SI-04
R02-CW-SI-05
R02-CW-SI-06
R02-CW-SI-07
R02-CW-SI-08
R02-CW-SI-09
R02-CW-SI-10
R02-CW-SI-11
R02-CW-SI-12
R02-CW-SI-13
R02-CW-SI-14
R02-CW-SI-15
R02-CW-SI-16
R02-CW-SI-17
R02-CW-SI-18
R02-CW-SI-19
R02-CW-SI-20
R02-CW-SI-21
R02-CW-SI-22
R02-CW-SI-23
R02-CW-SI-24
R02-CW-SI-25
R02-CW-SI-26
R02-CW-SI-27
R02-CW-SI-28
R02-CW-SI-29
R02-CW-SI-30
R02-SC-HA-01
R02-SC-HA-02
R02-SC-HA-03
R02-SC-HA-04
R02-SC-HA-05
R02-SC-HA-06
R02-SC-HA-07
R02-SC-HA-08
R02-SC-HA-09
R02-SC-HA-10
R02-SC-HA-11
R02-SC-HA-12
R02-SC-HA-13
R02-SC-HA-14
R02-SC-HA-15
R02-SC-HA-16
R02-SC-HA-17
R02-SC-HA-18
R02-SC-HA-19
R02-SC-HA-20
R02-SC-HA-21
R02-SC-HA-22
R02-SC-HA-23
R02-SC-HA-24
R02-SC-HA-25
R02-SC-HA-26
R02-SC-HA-27
R02-SC-HA-28
R02-SC-HA-29
R02-SC-HA-30
R02-CW-JO-01
R02-CW-JO-02
R02-CW-JO-03
R02-CW-JO-04
R02-CW-JO-05
R02-CW-JO-06
R02-CW-JO-07
R02-CW-JO-08
R02-CW-JO-09
R02-CW-JO-10
R02-CW-JO-11
R02-CW-JO-12
R02-CW-JO-13
R02-CW-JO-14
R02-CW-JO-15
R02-CW-JO-16
R02-CW-JO-17
R02-CW-JO-18
R02-CW-JO-19
R02-CW-JO-20
R02-CW-JO-21
R02-CW-JO-22
R02-CW-JO-23
R02-CW-JO-24
R02-CW-JO-25"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1

#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/B501_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/B501_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_B501.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




