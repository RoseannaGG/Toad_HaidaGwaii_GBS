#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane2_B502_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_B502_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-CW-JO-26
R02-CW-JO-27
R02-CW-JO-28
R02-CW-JO-29
R02-CW-JO-30
R02-CW-AL-01
R02-CW-AL-02
R02-CW-AL-03
R02-CW-AL-04
R02-CW-AL-05
R02-CW-AL-06
R02-CW-AL-07
R02-CW-AL-08
R02-CW-AL-09
R02-CW-AL-10
R02-CW-AL-11
R02-CW-AL-12
R02-CW-AL-13
R02-CW-AL-14
R02-CW-AL-15
R02-CW-AL-16
R02-CW-AL-17
R02-CW-AL-18
R02-CW-AL-19
R02-CW-AL-20
R02-CW-AL-21
R02-CW-AL-22
R02-CW-AL-23
R02-CW-AL-24
R02-CW-AL-25
R02-CW-AL-26
R02-CW-AL-27
R02-CW-AL-28
R02-CW-AL-29
R02-CW-LA-01
R02-CW-LA-02
R02-CW-LA-03
R02-CW-LA-04
R02-CW-LA-05
R02-CW-LA-06
R02-CW-LA-07
R02-CW-LA-08
R02-CW-LA-09
R02-CW-LA-10
R02-CW-LA-11
R02-CW-LA-12
R02-CW-LA-13
R02-CW-LA-14
R02-CW-LA-15
R02-CW-LA-16
R02-CW-LA-17
R02-CW-LA-18
R02-CW-LA-19
R02-CW-LA-20
R02-CW-LA-21
R02-CW-LA-22
R02-CW-LA-23
R02-CW-LA-24
R02-CW-LA-25
R02-CW-LA-26
R02-CW-LA-27
R02-CW-LA-28
R02-CW-LA-29
R02-CW-LA-30
R02-CW-BE-01
R02-CW-BE-02
R02-CW-BE-03
R02-CW-BE-04
R02-CW-BE-05
R02-CW-BE-06
R02-CW-BE-07
R02-CW-BE-08
R02-CW-BE-09
R02-CW-BE-10
R02-CW-BE-11
R02-CW-BE-12
R02-CW-BE-13
R02-CW-BE-14
R02-CW-BE-15
R02-CW-BE-16
R02-CW-BE-17
R02-CW-BE-18
R02-CW-BE-19
R02-CW-BE-20
R02-CW-BE-21
R02-CW-BE-22
R02-CW-BE-23
R02-CW-BE-24
R02-CW-BE-25
R02-CW-BE-26
R02-CW-BE-27
R02-CW-BE-28
R02-CW-BE-29
R02-CW-BE-30"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/B502_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/B502_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_B502.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




