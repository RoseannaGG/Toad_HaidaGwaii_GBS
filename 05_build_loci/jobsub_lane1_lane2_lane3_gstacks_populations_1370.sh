#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_lane1_lane2_lane3_gstacks_populations_1370
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_lane1_lane2_lane3_gstacks_populations_1370_%j.out




# rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/Pop_map_HGthesis_lane1_allexcptexperimnt_andnoLasqIs_lane2_LM_VanIs_lane3_NW_siteID.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/info/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/jobsub_lane1_lane2_lane3_gstacks_populations_1370.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress


module load StdEnv/2020
module load stacks/2.62

src=/scratch/roseanna/
cpu=24


# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.

echo “Starting gstacks at: `date`”

# increase the amount of files being able to open at one time - default is 1024, i have 1370 files
ulimit -n 4096


gstacks -I $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ --min-mapq 20 -t $cpu


# mkdir populations_ANBOref_maxhet0.6_writesinglesnp populations_ANBOref_writesinglesnp populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r10_R70pctoverall_p42_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r10_R70pctoverall_p47_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r10_R70pctoverall_p34_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r40_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r60_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r80_minmaf0.01_maxhet0.6_writesinglesnp populations_ANBOref_r60_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp



echo “Starting populations_ANBOref_r60_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu


echo “Starting populations_ANBOref_r60_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_r40_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r40_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-per-pop 0.4 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_r80_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r80_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-per-pop 0.8 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu


echo “Starting populations_ANBOref_r10_R70pctoverall_p34_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r10_R70pctoverall_p34_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 34 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu



echo “Starting populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 44 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu



echo “Starting populations_ANBOref_r10_R70pctoverall_p47_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r10_R70pctoverall_p47_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 47 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_maxhet0.6_writesinglesnp at: `date`”

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.


populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu





echo “Starting populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r10_R70pctoverall_p44_minmaf0.01_writesinglesnp/ --min-samples-overall 0.7 --min-samples-per-pop 0.1 --min-populations 44 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu



echo “Starting populations_ANBOref_writesinglesnp no filter except write single snp at: `date`”

#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.


populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_writesinglesnp/ --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu



echo “Job finished with exit code $? at: `date`”


## download populations outputs to: F:\GBS_data_03_02_21\Lane_1_2_3_feb2023\HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW