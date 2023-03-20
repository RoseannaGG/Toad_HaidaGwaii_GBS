#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --job-name=jobsub_HaidaGwaiifeb2023_ref_Aboreas_populations
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobsub_HaidaGwaiifeb2023_ref_Aboreas_populations_%j.out 

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/jobsub_HaidaGwaiifeb2023_ref_Aboreas_populations.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/jobsub_HaidaGwaiifeb2023_ref_Aboreas_populations.sh roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/ --progress

# rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/jobsub_HaidaGwaiifeb2023_ref_Aboreas_populations_35213589.out /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/ --progress


# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_385samples_June2021.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_385samples_June2021.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/info/ --progress

# rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/Pop_map_HaidaGwaii_385samples_June2021.tsv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/info/ --progress


src=/scratch/roseanna/
cpu=32

#files=/scratch/roseanna/info/list_test_samples_20.txt
#echo $files
#all_lines=`cat $files`


#turn txt file into array 
#myArray=(`cat "$files"`)

module load StdEnv/2020 stacks/2.62



# mkdir pop_ANBOref_r60_R60pctoverall_mm001_mh06




echo “Starting pop_ANBOref_r60_R60pctoverall_mm001_mh06 at: `date`”

populations -P $src/HaidaGwaii_Aboreas/gstacks_minmapq20/ -M $src/info/Pop_map_HaidaGwaii_385samples_June2021.tsv -O $src/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --fstats --hwe --smooth --smooth-fstats --vcf --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting pop_ANBOref_r60_R60pctoverall_mm001_mh06_vcfall at: `date`”

populations -P $src/HaidaGwaii_Aboreas/gstacks_minmapq20/ -M $src/info/Pop_map_HaidaGwaii_385samples_June2021.tsv -O $src/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_vcfall/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu


echo “Job finished with exit code $? at: `date`”


## download populations outputs to: F:\GBS_data_03_02_21\Lane_1_2_3_feb2023\HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW