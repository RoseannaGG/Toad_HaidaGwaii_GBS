#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D706_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D706_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R06-GH-PQ-01
R06-GH-PQ-02
R06-GH-PQ-03
R06-GH-PQ-04
R06-GH-PQ-05
R06-GH-PQ-06
R06-GH-PQ-07
R06-GH-PQ-08
R06-GH-PQ-09
R06-GH-PQ-10
R06-GH-PQ-11
R06-GH-PQ-12
R06-GH-PQ-13
R06-GH-PQ-14
R06-GH-PQ-15
R06-GH-PQ-16
R06-GH-PQ-17
R06-GH-PQ-18
R06-GH-PQ-19
R06-GH-PQ-20
R06-GH-PQ-21
R06-GH-PQ-22
R06-GH-PQ-23
R06-GH-PQ-24
R06-GH-PQ-25
R06-GH-GW-01
R06-GH-GW-02
R06-GH-GW-03
R06-GH-GW-04
R06-GH-GW-05
R06-GH-GW-06
R06-GH-GW-07
R06-GH-GW-08
R06-GH-GW-09
R06-GH-GW-10
R06-GH-GW-11
R06-GH-GW-12
R06-GH-GW-13
R06-GH-GW-14
R06-GH-GW-15
R06-GH-GW-16
R06-GH-GW-17
R06-GH-GW-18
R06-GH-GW-19
R06-GH-GW-20
R06-GH-GW-21
R06-GH-GW-22
R06-GH-GW-23
R06-GH-GW-24
R06-GH-GW-25
R06-GH-GW-26
R06-GH-GW-27
R06-GH-GW-28
R06-GH-GW-29
R06-GH-GW-30
R06-GH-LL-01
R06-GH-LL-02
R06-GH-LL-03
R06-GH-LL-04
R06-GH-LL-05
R06-GH-LL-06
R06-GH-LL-07
R06-GH-LL-08
R06-GH-LL-09
R06-GH-LL-10
R06-GH-LL-11
R06-GH-LL-12
R06-GH-LL-13
R06-GH-LL-14
R06-GH-LL-15
R06-GH-LL-16
R06-GH-LL-17
R06-GH-LL-18
R06-GH-LL-19
R06-GH-LL-20
R06-GH-LL-21
R06-GH-LL-22
R06-GH-LL-23
R06-GH-LL-24
R06-GH-LL-25
R06-GH-LL-26
R06-GH-LL-27
R06-GH-LL-28
R06-GH-LL-29
R06-GH-LL-30
R06-GI-LV-01
R06-GI-LV-02
R06-GI-LV-03
R06-GI-LV-04
R06-GI-LV-05
R06-GI-LV-06
R06-GI-LV-07
R06-GI-LV-08
R06-GI-LV-09
R06-GI-LV-10
R06-GI-RR-01"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D706_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D706_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D706.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




