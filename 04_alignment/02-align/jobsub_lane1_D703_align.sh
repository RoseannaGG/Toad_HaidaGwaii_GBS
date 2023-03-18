#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D703_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D703_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R06-GH-PT-23
R06-GH-PT-24
R06-GH-PT-25
R06-GH-PT-26
R06-GH-PT-27
R06-GH-PT-28
R06-GH-PT-29
R06-GH-PT-30
R02-SS-EF-11
R02-SS-EF-12
R02-SS-EF-13
R02-SS-EF-14
R02-SS-EF-15
R02-SS-EF-16
R02-SS-EF-17
R02-SS-EF-18
R02-SS-EF-19
R02-SS-EF-20
R02-SS-EF-21
R02-SS-EF-22
R02-SS-EF-23
R02-SS-EF-24
R02-SS-EF-25
R02-SS-EF-26
R02-SS-EF-27
R02-SS-EF-28
R02-SS-EF-29
R02-SS-EF-30
R06-GI-CK-11
R06-GI-CK-12
R06-GI-CK-13
R06-GI-CK-14
R06-GI-CK-15
R06-GI-CK-16
R06-GI-CK-17
R06-GI-CK-18
R06-GI-CK-19
R06-GI-CK-20
R06-GI-CK-21
R06-GI-CK-22
R06-GI-CK-23
R06-GI-CK-24
R06-GI-CK-25
R06-GI-CK-26
R06-GI-CK-27
R06-GI-CK-28
R06-GI-CK-29
R06-GI-CK-30
R06-GI-GL-11
R06-GI-GL-12
R06-GI-GL-13
R06-GI-GL-14
R06-GI-GL-15
R06-GI-GL-16
R06-GI-GL-17
R06-GI-GL-18
R06-GI-GL-19
R06-GI-GL-20
R06-GI-GL-21
R06-GI-GL-22
R06-GI-GL-23
R06-GI-GL-24
R06-GI-GL-25
R06-GI-GL-26
R06-GI-GL-27
R06-GI-GL-28
R06-GI-GL-29
R06-GI-GL-30
R02-SS-LC-11
R02-SS-LC-12
R02-SS-LC-13
R02-SS-LC-14
R02-SS-LC-15
R02-SS-LC-16
R02-SS-LC-17
R02-SS-LC-18
R02-SS-LC-19
R02-SS-LC-20
R02-SS-LC-21
R02-SS-LC-22
R02-SS-LC-23
R02-SS-LC-24
R02-SS-LC-25
R02-SS-LC-26
R02-SS-LC-27
R02-SS-LC-28
R02-SS-LC-29
R02-SS-LC-30
R02-SS-LT-11
R02-SS-LT-12
R02-SS-LT-13
R02-SS-LT-14
R02-SS-LT-15
R02-SS-LT-16
R02-SS-LT-17
R02-SS-LT-18"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D703_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D703_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D703.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




