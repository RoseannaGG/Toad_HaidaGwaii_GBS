#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_D503_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_D503_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R01-SI-GL-07
R01-SI-GL-08
R01-SI-GL-09
R01-SI-GL-10
R01-SI-GL-11
R01-SI-GL-12
R01-SI-GL-13
R01-SI-GL-14
R01-SI-GL-15
R01-SI-GL-16
R01-SI-GL-17
R01-SI-GL-18
R01-SI-GL-19
R01-SI-GL-20
R01-SI-GL-21
R01-SI-GL-22
R01-SI-GL-23
R01-SI-GL-24
R01-SI-GL-25
R01-SI-GL-26
R01-SI-GL-27
R01-SI-GL-28
R01-SI-GL-29
R01-SI-GL-30
R01-NI-OC-01
R01-NI-OC-02
R01-NI-OC-03
R01-NI-OC-04
R01-NI-OC-05
R01-NI-OC-06
R01-NI-OC-07
R01-NI-OC-08
R01-NI-OC-09
R01-NI-OC-10
R01-NI-OC-11
R01-NI-OC-12
R01-NI-OC-13
R01-NI-OC-14
R01-NI-OC-15
R01-NI-OC-16
R01-NI-OC-17
R01-NI-OC-18
R01-NI-OC-19
R01-NI-OC-20
R01-NI-OC-21
R01-NI-OC-22
R01-NI-OC-23
R01-NI-OC-24
R01-NI-OC-25
R01-NI-OC-26
R01-NI-OC-27
R01-NI-OC-28
R01-NI-OC-29
R01-NI-OC-30
R01-SI-LC-01
R01-SI-LC-02
R01-SI-LC-03
R01-SI-LC-04
R01-SI-LC-05
R01-SI-LC-06
R01-SI-LC-07
R01-SI-LC-08
R01-SI-LC-09
R01-SI-LC-10
R01-SI-LC-11
R01-SI-LC-12
R01-SI-LC-13
R01-SI-LC-14
R01-SI-LC-15
R01-SI-LC-16
R01-SI-LC-17
R01-SI-LC-18
R01-SI-LC-19
R01-SI-LC-20
R01-SI-LC-21
R01-SI-LC-22
R01-SI-LC-23
R01-SI-LC-24
R01-SI-LC-25
R01-SI-LC-26
R01-SI-LC-27
R01-SI-LC-28
R01-SI-LC-29
R01-SI-LC-30
R01-SI-HR-01
R01-SI-HR-02
R01-SI-HR-03
R01-SI-HR-04
R01-SI-HR-05
R01-SI-HR-06
R01-SI-HR-07
R01-SI-HR-08
R01-SI-HR-09
R01-SI-HR-10
R01-SI-HR-11
R01-SI-HR-12"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D503_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D503_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_D503.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




