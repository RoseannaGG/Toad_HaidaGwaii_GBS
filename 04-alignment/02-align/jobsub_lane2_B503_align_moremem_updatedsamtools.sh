#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_lane2_B503_align_moremem_updatedsamtools
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane2_B503_align_moremem_updatedsamtools_%j.out


# transfer forward reads from B503 to /Demultiplexing_stacks/B503_norenz2_trimR2/


# rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/jobsub_lane2_B503_align_moremem_updatedsamtools.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

files="R02-SS-KH-01
R02-SS-KH-02
R02-SS-KH-03
R02-SS-KH-04
R02-SS-KH-05
R02-SS-KH-06
R02-SS-KH-07
R02-SS-KH-08
R02-SS-KH-09
R07-OM-WB-01
R07-OM-WB-02
R07-OM-WB-03
R07-OM-WB-04
R07-OM-WB-05
R07-OM-WB-06
R07-OM-WB-07
R07-OM-WB-08
R07-OM-WB-09
R07-OM-WB-10
R07-OM-WB-11
R07-OM-WB-12
R07-OM-WB-13
R07-OM-WB-14
R07-OM-WB-15
R07-OM-WB-16
R07-OM-WB-17
R07-OM-WB-18
R07-OM-WB-19
R07-OM-GG-01
R07-NE-AS-01
R07-NE-AS-02
R07-NE-AS-03
R07-NE-AS-04
R07-NE-AS-05
R07-NE-AS-06
R07-NE-AS-07
R07-NE-AS-08
R07-NE-AS-09
R07-NE-AS-10
R07-NE-AS-11
R07-NE-AS-12
R07-NE-AS-13
R07-NE-AS-14
R07-NE-AS-15
R07-NE-AS-16
R07-NE-AS-17
R07-NE-AS-18
R07-NE-AS-19"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.16.1


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/B503_norenz2_trimR2/${sample}.1.fq.gz $src/Demultiplexing_stacks/B503_norenz2_trimR2/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane2_B503__moremem_updatedsamtools.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignmentsmapq20_B503_moremem_updatedsamtools/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




