#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_HaidaGwaiifeb2023_ref_Aboreas
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_HaidaGwaiifeb2023_ref_Aboreas_%j.out 

#mkdir HaidaGwaii
#mkdir populations_maxhet0.6_writesinglesnp populations_writesinglesnp

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/jobsub_HaidaGwaiifeb2023_ref_Aboreas.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_385samples_June2021.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_385samples_June2021.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/info/ --progress

files="R06-GH-DL-01
R06-GH-DL-02
R06-GH-DL-03
R06-GH-DL-04
R06-GH-DL-05
R06-GH-DL-06
R06-GH-DL-07
R06-GH-DL-08
R06-GH-DL-09
R06-GH-DL-10
R06-GH-DL-11
R06-GH-DL-12
R06-GH-DL-13
R06-GH-DL-14
R06-GH-DL-15
R06-GH-DL-16
R06-GH-DL-17
R06-GH-DL-18
R06-GH-DL-19
R06-GH-DL-20
R06-GH-DL-21
R06-GH-DL-22
R06-GH-DL-23
R06-GH-DL-24
R06-GH-DL-25
R06-GH-DL-26
R06-GH-DL-27
R06-GH-DL-28
R06-GH-DL-29
R06-GH-DL-30
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
R06-GH-PQ-01
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
R06-GH-PT-22
R06-GH-PT-23
R06-GH-PT-24
R06-GH-PT-25
R06-GH-PT-26
R06-GH-PT-27
R06-GH-PT-28
R06-GH-PT-29
R06-GH-PT-30
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
R06-GI-EV-11
R06-GI-EV-12
R06-GI-EV-13
R06-GI-EV-14
R06-GI-EV-15
R06-GI-EV-16
R06-GI-EV-17
R06-GI-EV-18
R06-GI-EV-19
R06-GI-EV-20
R06-GI-EV-21
R06-GI-EV-22
R06-GI-EV-23
R06-GI-EV-24
R06-GI-EV-25
R06-GI-EV-26
R06-GI-EV-27
R06-GI-EV-28
R06-GI-EV-29
R06-GI-EV-30
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
R06-GI-MY-11
R06-GI-MY-12
R06-GI-MY-13
R06-GI-MY-14
R06-GI-MY-15
R06-GI-MY-16
R06-GI-MY-17
R06-GI-MY-18
R06-GI-MY-19
R06-GI-MY-20
R06-GI-MY-21
R06-GI-MY-22
R06-GI-MY-23
R06-GI-MY-24
R06-GI-MY-25
R06-GI-MY-26
R06-GI-MY-27
R06-GI-MY-28
R06-GI-MY-29
R06-GI-MY-30
R06-GI-RR-01
R06-GI-RR-02
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
R06-MI-MM-30"


src=/scratch/roseanna/
bwa_db=$src/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas
cpu=24

#files=/scratch/roseanna/info/list_test_samples_20.txt
#echo $files
#all_lines=`cat $files`


#turn txt file into array 
#myArray=(`cat "$files"`)

module load StdEnv/2020 stacks/2.62

module load bwa/0.7.17
module load samtools/1.16.1



echo “Starting alignment at: `date`”


#
# Align paired-end data with BWA, convert to BAM and SORT.
#
for sample in $files;
do
    bwa mem -M -t $cpu $bwa_db $src/Demultiplexing_stacks/all3lanessamples/${sample}.1.fq.gz $src/Demultiplexing_stacks/all3lanessamples/${sample}.2.fq.gz 2> $src/HaidaGwaii_Aboreas/logs/bwa_HG_ref_Aboreas.err|
      samtools view -b -q 20|
      samtools sort --threads $cpu > $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/${sample}.bam;
done

echo “Job finished with exit code $? at: `date`”


echo “Starting gstacks at: `date`”

#
# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.
#

gstacks -I $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/ -M $src/info/Pop_map_HaidaGwaii_385samples_June2021.tsv -O $src/HaidaGwaii_Aboreas/gstacks_minmapq20/ --min-mapq 20 -t $cpu

echo “Starting populations_maxhet0.6_writesinglesnp at: `date`”

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.

#
#populations -P $src/HaidaGwaii_Aboreas/gstacks_minmapq20/ -M $src/info/Pop_map_HaidaGwaii_385samples_June2021.tsv -O $src/HaidaGwaii_Aboreas/gstacks_minmapq20/populations_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --write-single-snp --fstats --hwe --vcf --plink --structure --genepop --treemix --log-fst-comp --verbose -t $cpu



#echo “Starting populations no filter except write single snp at: `date`”

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.

#
#populations -P $src/HaidaGwaii_Aboreas/gstacks_minmapq20/ -M $src/info/Pop_map_HaidaGwaii_385samples_June2021.tsv -O $src/HaidaGwaii_Aboreas/gstacks_minmapq20/populations_writesinglesnp/ --write-single-snp --fstats --hwe --vcf --plink --structure --genepop --treemix --log-fst-comp --verbose -t $cpu



echo “Job finished with exit code $? at: `date`”