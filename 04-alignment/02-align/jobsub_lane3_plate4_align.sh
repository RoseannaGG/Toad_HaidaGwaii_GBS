#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate4_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate4_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R06-SK-FE-23
R06-SK-FE-24
R06-SK-FE-25
R06-SK-FE-26
R06-SK-FE-27
R06-SK-FE-28
R06-SK-FE-29
R06-SK-FE-30
R06-SK-ON-01
R06-SK-ON-02
R06-SK-ON-03
R06-SK-ON-04
R06-SK-ON-05
R06-SK-ON-06
R06-SK-ON-07
R06-SK-ON-08
R06-SK-ON-09
R06-SK-ON-10
R06-SK-ON-11
R06-SK-ON-12
R06-SK-ON-13
R06-SK-ON-14
R06-SK-ON-15
R06-SK-ON-16
R06-SK-ON-17
R06-SK-ON-18
R06-SK-ON-19
R06-SK-ON-20
R06-SK-ON-21
R06-SK-ON-22
R06-SK-ON-23
R06-SK-ON-24
R06-SK-ON-25
R06-SK-ON-26
R06-SK-ON-27
R06-SK-ON-28
R06-SK-ON-29
R06-SK-ON-30
R06-SK-ML-01
R06-SK-ML-02
R06-SK-ML-03
R06-SK-ML-04
R06-SK-ML-05
R06-SK-ML-06
R06-SK-ML-07
R06-SK-ML-08
R06-SK-ML-09
R06-SK-ML-10
R06-SK-ML-11
R06-SK-ML-12
R06-SK-ML-13
R06-SK-ML-14
R06-SK-ML-15
R06-SK-ML-16
R06-SK-ML-17
R06-SK-ML-18
R07-OM-FF-01
R07-OM-FF-02
R07-OM-FF-03
R07-OM-FF-04
R07-OM-FF-05
R07-OM-FF-06
R07-OM-FF-07
R07-OM-FF-08
R07-OM-FF-09
R07-OM-FF-10
R07-OM-FF-11
R07-OM-FF-12
R07-OM-FF-13
R07-OM-FF-14
R07-OM-FF-15
R07-OM-FF-16
R07-OM-FF-17
R07-OM-FF-18
R07-OM-FF-19
R07-OM-FF-20
R07-OM-FF-21
R07-OM-FF-22
R07-OM-FF-23
R07-OM-FF-24
R07-OM-FF-25
R07-OM-FF-26
R07-OM-FF-27
R07-OM-FF-28
R07-OM-FF-29
R07-OM-FF-30
R07-OM-EE-01
R07-OM-EE-02
R07-OM-EE-03
R07-OM-EE-04
R07-OM-EE-05
R07-OM-EE-06
R07-OM-EE-07
R07-OM-EE-08
R07-OM-EE-09
R07-OM-EE-10"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate4_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate4/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate4/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate4.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




