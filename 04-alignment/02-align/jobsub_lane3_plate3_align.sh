#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=80000
#SBATCH --job-name=jobsub_lane3_plate3_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane3_plate3_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R05-CB-HU-14
R05-CB-HU-15
R05-CB-HU-16
R05-CB-HU-17
R05-CB-HU-18
R05-CB-HU-19
R05-CB-HU-20
R05-CB-HU-21
R05-CB-HU-22
R05-CB-HU-23
R05-CB-HU-24
R05-CB-HU-25
R05-CB-HU-26
R05-CB-HU-27
R05-CB-HU-28
R05-CB-HU-29
R05-CB-HU-30
R05-CB-AL-01
R05-CB-AL-02
R05-CB-AL-03
R05-CB-AL-04
R05-CB-AL-05
R05-CB-AL-06
R05-CB-AL-07
R05-CB-AL-08
R05-CB-AL-09
R05-CB-AL-10
R05-CB-AL-11
R05-CB-AL-12
R05-CB-AL-13
R05-CB-AL-14
R05-CB-AL-15
R05-CB-AL-16
R05-CB-AL-17
R05-CB-AL-18
R05-CB-AL-19
R05-CB-AL-20
R05-CB-AL-21
R05-CB-AL-22
R05-CB-AL-23
R05-CB-AL-24
R05-CB-AL-25
R05-CB-AL-26
R05-CB-AL-27
R05-CB-AL-28
R05-CB-AL-29
R05-CB-AL-30
R06-SK-SI-01
R06-SK-SI-02
R06-SK-SI-03
R06-SK-SI-04
R06-SK-SI-05
R06-SK-RO-01
R06-SK-RO-02
R06-SK-RO-03
R06-SK-RO-04
R06-SK-RO-05
R06-SK-RO-06
R06-SK-RO-07
R06-SK-RO-08
R06-SK-RO-09
R06-SK-RO-10
R06-SK-RO-11
R06-SK-RO-12
R06-SK-RO-13
R06-SK-RO-14
R06-SK-RO-15
R06-SK-RO-16
R06-SK-RO-17
R06-SK-RO-18
R06-SK-RO-19
R06-SK-RO-20
R06-SK-RO-21
R06-SK-RO-22
R06-SK-FE-01
R06-SK-FE-02
R06-SK-FE-03
R06-SK-FE-04
R06-SK-FE-05
R06-SK-FE-06
R06-SK-FE-07
R06-SK-FE-08
R06-SK-FE-09
R06-SK-FE-10
R06-SK-FE-11
R06-SK-FE-12
R06-SK-FE-13
R06-SK-FE-14
R06-SK-FE-15
R06-SK-FE-16
R06-SK-FE-17
R06-SK-FE-18
R06-SK-FE-19
R06-SK-FE-20
R06-SK-FE-21
R06-SK-FE-22"

# rsync -zv /drives/f/GBS_data_03_02_21/Lane3_GBSdata_2023/jobsub_lane3_plate3_align.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/lane3_plate3/${sample}.1.fq.gz $src/Demultiplexing_stacks/lane3_plate3/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane3_plate3.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




