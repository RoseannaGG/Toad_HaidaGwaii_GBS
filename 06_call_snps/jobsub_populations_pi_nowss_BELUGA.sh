#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=168:00:00
#SBATCH --mem=512GB
#SBATCH --job-name=jobsub_populations_pi_nowss_BELUGA_%j
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_populations_pi_nowss_BELUGA_%j.out




# rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/Pop_map_HGthesis_lane1_allexcptexperimnt_andnoLasqIs_lane2_LM_VanIs_lane3_NW_siteID.tsv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/info/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/jobsub_populations_pi_nowss_BELUGA.sh roseanna@beluga.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/ --progress

# rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/jobsub_populations_pi_nowss_BELUGA_ /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/ --progress


module load StdEnv/2020
module load stacks/2.62

src=/scratch/roseanna/
cpu=32


# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.

#echo “Starting gstacks at: `date`”

# increase the amount of files being able to open at one time - default is 1024, i have 1370 files
ulimit -n 4096


#gstacks -I $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ --min-mapq 20 -t $cpu


# mkdir populations_ANBOref_r60_R60pctoverall_mm001_mh06_novcfall populations_ANBOref_r60_R60pctoverall_mm001_mh06



echo “Starting populations_ANBOref_r60_R60pctoverall_mm001_mh06 at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_r60_R60pctoverall_mm001_mh06_vcfall at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_vcfall/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu


echo “Job finished with exit code $? at: `date`”


## download populations outputs to: F:\GBS_data_03_02_21\Lane_1_2_3_feb2023\HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW