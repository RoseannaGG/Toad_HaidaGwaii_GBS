#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane1_D702_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D702_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-SS-LT-07
R02-SS-LT-08
R02-SS-LT-09
R02-SS-LT-10
R06-GI-CN-01
R06-GI-CN-02
R06-GI-CN-03
R06-GI-CN-04
R06-GI-CN-05
R06-GI-CN-06
R06-GI-CN-07
R06-GI-CN-08
R06-GI-CN-09
R06-GI-CN-10
R06-GH-PT-01
R06-GH-PT-02
R06-GH-PT-03
R06-GH-PT-04
R06-GH-PT-05
R06-GH-PT-06
R06-GH-PT-07
R06-GH-PT-08
R06-GH-PT-09
R06-GH-PT-10
R02-SS-EF-01
R02-SS-EF-02
R02-SS-EF-03
R02-SS-EF-04
R02-SS-EF-05
R02-SS-EF-06
R02-SS-EF-07
R02-SS-EF-08
R02-SS-EF-09
R02-SS-EF-10
R06-GI-CK-01
R06-GI-CK-02
R06-GI-CK-03
R06-GI-CK-04
R06-GI-CK-05
R06-GI-CK-06
R06-GI-CK-07
R06-GI-CK-08
R06-GI-CK-09
R06-GI-CK-10
R02-SS-LC-01
R02-SS-LC-02
R02-SS-LC-03
R02-SS-LC-04
R02-SS-LC-05
R02-SS-LC-06
R02-SS-LC-07
R02-SS-LC-08
R02-SS-LC-09
R02-SS-LC-10
R06-GI-MY-01
R06-GI-MY-02
R06-GI-MY-03
R06-GI-MY-04
R06-GI-MY-05
R06-GI-MY-06
R06-GI-MY-07
R06-GI-MY-08
R06-GI-MY-09
R06-GI-MY-10
R06-GI-CN-11
R06-GI-CN-12
R06-GI-CN-13
R06-GI-CN-14
R06-GI-CN-15
R06-GI-CN-16
R06-GI-CN-17
R06-GI-CN-18
R06-GI-CN-19
R06-GI-CN-20
R06-GI-CN-21
R06-GI-CN-22
R06-GI-CN-23
R06-GI-CN-24
R06-GI-CN-25
R06-GI-CN-26
R06-GI-CN-27
R06-GI-CN-28
R06-GI-CN-29
R06-GI-CN-30
R06-GH-PT-11
R06-GH-PT-12
R06-GH-PT-13
R06-GH-PT-14
R06-GH-PT-15
R06-GH-PT-16
R06-GH-PT-17
R06-GH-PT-18
R06-GH-PT-19
R06-GH-PT-20
R06-GH-PT-21
R06-GH-PT-22"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D702_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D702_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D702.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




