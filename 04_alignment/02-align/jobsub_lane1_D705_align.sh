#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D705_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D705_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-CW-MI-15
R02-CW-MI-16
R02-CW-MI-17
R02-CW-MI-18
R02-CW-MI-19
R02-CW-MI-20
R02-CW-MI-21
R02-CW-MI-22
R02-CW-MI-23
R02-CW-MI-24
R02-CW-MI-25
R02-CW-MI-26
R02-CW-MI-27
R02-CW-MI-28
R02-CW-MI-29
R02-CW-MI-30
R02-CW-KA-11
R02-CW-KA-12
R02-CW-KA-13
R02-CW-KA-14
R02-CW-KA-15
R02-CW-KA-16
R02-CW-KA-17
R02-CW-KA-18
R02-CW-KA-19
R02-CW-KA-20
R02-CW-KA-21
R02-CW-KA-22
R02-CW-KA-23
R02-CW-KA-24
R02-CW-KA-25
R02-CW-KA-26
R02-CW-KA-27
R02-CW-KA-28
R02-CW-KA-29
R02-CW-KA-30
R02-SC-IN-11
R02-SC-IN-12
R02-SC-IN-13
R02-SC-IN-14
R02-SC-IN-15
R02-SC-IN-16
R02-SC-IN-17
R02-SC-IN-18
R02-SC-IN-19
R02-SC-IN-20
R02-SC-IN-21
R02-SC-IN-22
R02-SC-IN-23
R02-SC-IN-24
R02-SC-IN-25
R02-SC-IN-26
R02-SC-IN-27
R02-SC-IN-28
R02-SC-IN-29
R02-SC-IN-30
R06-GH-DL-11
R06-GH-DL-12
R06-GH-DL-13
R06-GH-DL-14
R06-GH-DL-15
R06-GH-DL-16
R06-GH-DL-17
R06-GH-DL-18
R06-GH-DL-19
R06-GH-DL-20
R06-GH-DL-21
R06-GH-DL-22
R06-GH-DL-23
R06-GH-DL-24
R06-GH-DL-25
R06-GH-DL-26
R06-GH-DL-27
R06-GH-DL-28
R06-GH-DL-29
R06-GH-DL-30
R06-GI-EV-11
R06-GI-EV-12
R06-GI-EV-13
R06-GI-EV-14
R06-GI-EV-15
R06-GI-EV-16
R06-GI-EV-17
R06-GI-EV-18
R06-GI-EV-19
R06-GI-EV-20
R06-GI-EV-21
R06-GI-EV-22
R06-GI-EV-23
R06-GI-EV-24
R06-GI-EV-25
R06-GI-EV-26
R06-GI-EV-27
R06-GI-EV-28
R06-GI-EV-29
R06-GI-EV-30"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D705_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D705_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D705.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




