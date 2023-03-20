#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --ntasks=65
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=jobSUB_structure_python_766_1368_k1to13_5reps_nodirforstrfile
#SBATCH --mail-user=roseanna.gamlen.greene@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=jobSUB_structure_python_766_1368_k1to13_5reps_nodirforstrfile_%j.out

module load StdEnv/2020 gcc/9.3.0 structure/2.3.4

cat /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/structure/job_structure_python_766_1368_k1to13_5reps_nodirforstrfile | parallel -j 65 --joblog structure_python_766_1368_k1to13_5reps_nodirforstrfile.log
