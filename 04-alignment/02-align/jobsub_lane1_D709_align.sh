#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D709_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D709_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R01-SI-FR-24
R01-SI-FR-25
R01-SI-FR-26
R01-SI-FR-27
R01-SI-FR-28
R01-SI-FR-29
R01-SI-FR-30
R01-SI-SJ-01
R01-SI-SJ-02
R01-SI-SJ-03
R01-SI-SJ-04
R01-SI-SJ-05
R01-SI-SJ-06
R01-SI-SJ-07
R01-SI-SJ-08
R01-SI-SJ-09
R01-SI-SJ-10
R01-SI-SJ-11
R01-SI-SJ-12
R01-SI-SJ-13
R01-SI-SJ-14
R01-SI-SJ-15
R01-SI-SJ-16
R01-SI-SJ-17
R01-SI-SJ-18
R01-SI-SJ-19
R01-SI-SJ-20
R01-SI-SJ-21
R01-SI-SJ-22
R01-SI-SJ-23
R01-SI-SJ-24
R01-SI-SJ-25
R01-SI-SJ-26
R01-SI-SJ-27
R01-SI-SJ-28
R01-SI-SJ-29
R01-SI-SJ-30
R01-NI-CX-01
R01-NI-CX-02
R01-NI-CX-03
R01-NI-CX-04
R01-NI-CX-05
R01-NI-CX-06
R01-NI-CX-07
R01-NI-CX-08
R01-NI-CX-09
R01-NI-CX-10
R01-NI-CX-11
R01-NI-CX-12
R01-NI-CX-13
R01-NI-CX-14
R01-NI-CX-15
R01-NI-CX-16
R01-NI-CX-17
R01-NI-CX-18
R01-NI-CX-19
R01-NI-CX-20
R01-NI-CX-21
R01-NI-CX-22
R01-NI-CX-23
R01-NI-CX-24
R01-NI-CX-25
R01-NI-CX-26
R01-NI-CX-27
R01-NI-CX-28
R01-NI-CX-29
R01-NI-CX-30"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D709_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D709_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D709.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




