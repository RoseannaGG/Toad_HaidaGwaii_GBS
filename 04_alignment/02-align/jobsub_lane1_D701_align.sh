#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_lane1_D701_align
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_D701_align_%j.out

src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=48

files="R01-SI-LL-01
R01-SI-LL-02
R01-SI-LL-03
R01-SI-LL-04
R01-SI-LL-05
R01-SI-LL-06
R01-SI-LL-07
R01-SI-LL-08
R01-SI-LL-09
R01-SI-LL-10
R01-CR-CE-01
R01-CR-CE-02
R01-CR-CE-03
R01-CR-CE-04
R01-CR-CE-05
R01-CR-CE-06
R01-CR-CE-07
R01-CR-CE-08
R01-CR-CE-09
R01-CR-CE-10
R01-SI-CA-01
R01-SI-CA-02
R01-SI-CA-03
R01-SI-CA-04
R01-SI-CA-05
R01-SI-CA-06
R01-SI-CA-07
R01-SI-CA-08
R01-SI-CA-09
R01-SI-CA-10
R02-CW-MI-01
R02-CW-MI-02
R02-CW-MI-03
R02-CW-MI-04
R02-CW-MI-05
R02-CW-MI-06
R02-CW-MI-07
R02-CW-MI-08
R02-CW-MI-09
R02-CW-MI-10
R02-CW-KA-01
R02-CW-KA-02
R02-CW-KA-03
R02-CW-KA-04
R02-CW-KA-05
R02-CW-KA-06
R02-CW-KA-07
R02-CW-KA-08
R02-CW-KA-09
R02-CW-KA-10
R02-SC-IN-01
R02-SC-IN-02
R02-SC-IN-03
R02-SC-IN-04
R02-SC-IN-05
R02-SC-IN-06
R02-SC-IN-07
R02-SC-IN-08
R02-SC-IN-09
R02-SC-IN-10
R06-GH-DL-01
R06-GH-DL-02
R06-GH-DL-03
R06-GH-DL-04
R06-GH-DL-05
R06-GH-DL-06
R06-GH-DL-07
R06-GH-DL-08
R06-GH-DL-09
R06-GH-DL-10
R06-GI-EV-01
R06-GI-EV-02
R06-GI-EV-03
R06-GI-EV-04
R06-GI-EV-05
R06-GI-EV-06
R06-GI-EV-07
R06-GI-EV-08
R06-GI-EV-09
R06-GI-EV-10
R06-GI-GL-01
R06-GI-GL-02
R06-GI-GL-03
R06-GI-GL-04
R06-GI-GL-05
R06-GI-GL-06
R06-GI-GL-07
R06-GI-GL-08
R06-GI-GL-09
R06-GI-GL-10
R02-SS-LT-01
R02-SS-LT-02
R02-SS-LT-03
R02-SS-LT-04
R02-SS-LT-05
R02-SS-LT-06"



echo “Starting alignment at: `date`”


module load bwa/0.7.17
module load samtools/1.11


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/D701_norenz2_trimR2_PstI/${sample}.1.fq.gz $src/Demultiplexing_stacks/D701_norenz2_trimR2_PstI/${sample}.2.fq.gz 2> $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/logs/bwa_lane1_D701.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignmentsmapq20/${sample}.bam;
done

echo “Job alignment finished with exit code $? at: `date`”




