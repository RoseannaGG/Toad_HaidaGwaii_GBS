#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_D502_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_D502_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R01-SI-DB-01
R01-SI-DB-02
R01-SI-DB-03
R01-SI-DB-04
R01-SI-DB-05
R01-SI-DB-06
R01-SI-DB-07
R01-SI-DB-08
R01-SI-DB-09
R01-SI-DB-10
R01-SI-DB-11
R01-SI-DB-12
R01-SI-DB-13
R01-SI-DB-14
R01-SI-DB-15
R01-SI-DB-16
R01-SI-DB-17
R01-SI-DB-18
R01-SI-DB-19
R01-SI-DB-20
R01-SI-DB-21
R01-SI-DB-22
R01-SI-DB-23
R01-SI-DB-24
R01-SI-DB-25
R01-SI-DB-26
R01-SI-DB-27
R01-SI-DB-28
R01-SI-DB-29
R01-SI-DB-30
R01-SI-KE-01
R01-SI-KE-02
R01-SI-KE-03
R01-SI-KE-04
R01-SI-KE-05
R01-SI-KE-06
R01-SI-KE-07
R01-SI-KE-08
R01-SI-KE-09
R01-SI-KE-10
R01-SI-KE-11
R01-SI-KE-12
R01-SI-KE-13
R01-SI-KE-14
R01-SI-KE-15
R01-SI-KE-16
R01-SI-KE-17
R01-SI-KE-18
R01-SI-KE-19
R01-SI-KE-20
R01-SI-KE-21
R01-SI-KE-22
R01-SI-KE-23
R01-SI-KE-24
R01-SI-KE-25
R01-SI-KE-26
R01-SI-KE-27
R01-SI-KE-28
R01-SI-KE-29
R01-SI-KE-30
R01-CR-RA-01
R01-CR-RA-02
R01-CR-RA-03
R01-CR-RA-04
R01-CR-RA-05
R01-CR-RA-06
R01-CR-RA-07
R01-CR-RA-08
R01-CR-RA-09
R01-CR-RA-10
R01-CR-RA-11
R01-CR-RA-12
R01-CR-RA-13
R01-CR-RA-14
R01-CR-RA-15
R01-CR-RA-16
R01-CR-RA-17
R01-CR-RA-18
R01-CR-RA-19
R01-CR-RA-20
R01-CR-RA-21
R01-CR-RA-22
R01-CR-RA-23
R01-CR-RA-24
R01-CR-RA-25
R01-CR-RA-26
R01-CR-RA-27
R01-CR-RA-28
R01-CR-RA-29
R01-CR-RA-30
R01-SI-GL-01
R01-SI-GL-02
R01-SI-GL-03
R01-SI-GL-04
R01-SI-GL-05
R01-SI-GL-06"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D502_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/D502_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_D502.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




