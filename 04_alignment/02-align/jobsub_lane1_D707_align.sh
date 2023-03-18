#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D707_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D707_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R06-GI-RR-02
R06-GI-RR-03
R06-GI-RR-04
R06-GI-RR-05
R06-GI-RR-06
R06-GI-RR-07
R06-GI-RR-08
R06-GI-RR-09
R06-GI-RR-10
R06-GI-RR-11
R06-GI-RR-12
R06-GI-RR-13
R06-GI-RR-14
R06-GI-RR-15
R06-GI-RR-16
R06-GI-RR-17
R06-GI-RR-18
R06-GI-RR-19
R06-GI-RR-20
R06-GI-RR-21
R06-GI-RR-22
R06-GI-RR-23
R06-GI-RR-24
R06-GI-RR-25
R06-GI-RR-26
R06-GI-RR-27
R06-GI-RR-28
R06-GI-RR-29
R06-GI-RR-30
R06-MI-MM-01
R06-MI-MM-02
R06-MI-MM-03
R06-MI-MM-04
R06-MI-MM-05
R06-MI-MM-06
R06-MI-MM-07
R06-MI-MM-08
R06-MI-MM-09
R06-MI-MM-10
R06-MI-MM-11
R06-MI-MM-12
R06-MI-MM-13
R06-MI-MM-14
R06-MI-MM-15
R06-MI-MM-16
R06-MI-MM-17
R06-MI-MM-18
R06-MI-MM-19
R06-MI-MM-20
R06-MI-MM-21
R06-MI-MM-22
R06-MI-MM-23
R06-MI-MM-24
R06-MI-MM-25
R06-MI-MM-26
R06-MI-MM-27
R06-MI-MM-28
R06-MI-MM-29
R06-MI-MM-30
R06-GI-LV-11
R06-GI-LV-12
R06-GI-LV-13
R06-GI-LV-14
R06-GI-LV-15
R06-GI-LV-16
R06-GI-LV-17
R06-GI-LV-18
R06-GI-LV-19
R06-GI-LV-20
R06-GI-LV-21
R06-GI-LV-22
R06-GI-LV-23
R06-GI-LV-24
R06-GI-LV-25
R06-GI-LV-26
R06-GI-LV-27
R06-GI-LV-28
R06-GI-LV-29
R06-GI-LV-30
R02-SS-FA-01
R02-SS-FA-02
R02-SS-FA-03
R02-SS-FA-04
R02-SS-FA-05
R02-SS-FA-06
R02-SS-FA-07
R02-SS-FA-08
R02-SS-FA-09
R02-SS-FA-10
R02-SS-FA-11
R02-SS-FA-12
R02-SS-FA-13
R02-SS-FA-14
R02-SS-FA-15
R02-SS-FA-16
R02-SS-FA-17"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D707_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D707_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D707.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




