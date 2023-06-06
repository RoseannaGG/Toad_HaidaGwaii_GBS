1. create structure file in R ... set pops as regions -genind to strc function

2. input str file as table in R, and edit pop column to a number per pond, export the list of number per pond - and export these nunmbers for use later

3. open structure file in txt editor like sublime and see if pops are numbers per pond 

4. and delete "ind" and "pop" from header (or do in excel)  - and save as new file with "_nolab" at the end of name

5. and then make two tabs at start - and save should have two empty tabs at start of header 

6. should have one empty row after the last row with data - if not, make one, and save


7. save line endings as unix and encoding as UTF-8


8. file names CANNOT hvae dots in them


9. edit other files 

	a) main params file (loci 603 and inv 225 and popdata =1, label=1), 
	b) extra params file, 
	c) file with all the structure commands for k1 to k13, 
	d) and then the job script to run it on the cluster - remmeber to change number of nodes and parallel to number of lines in file c)

10. check if runs locally 

 copy str file, main parms and extra parms into folder in downloads folder 

cmd



cd C:\Users\Roseanna\Downloads\structure_march2023\225INDIV_603SNPS

C:\Users\Roseanna\Downloads\structure_windows_console\console\structure.exe -K 2 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o structure_225INDIV_603SNPS_SITEID_nolab_k2r1 -D 37847126

### SENDS OUTPUT TO LOG FILE
C:\Users\Roseanna\Downloads\structure_windows_console\console\structure.exe -K 2 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o structure_225INDIV_603SNPS_SITEID_nolab_k2r1 -D 37847122 > structure_225INDIV_603SNPS_SITEID_nolab_screenoutput_k2r1.log




11. upload files a-d and also the structure file to a folder on cluster

######## upload file _nodirforstrfile - i.e. has dir (file path) for all files except the structure file


structure_225INDIV_603SNPS_SITEID_nolab

rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/mainparams_225_603_SITEID_nolab /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/structure_225INDIV_603SNPS_SITEID_nolab.str /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/extraparams roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/ --progress




12. test run job on cluster

## make dir on cluster 
225INDIV_603SNPS/structure/

mkdir  outputs_k1to13_5rep joblogs_k1to13

cd /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/

module load nixpkgs/16.09 intel/2018.3 structure/2.3.4


structure_225INDIV_603SNPS_SITEID_nolab

structure -K 1 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122


"structure -K 1 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122


----------------------------------------------------
STRUCTURE by Pritchard, Stephens and Donnelly (2000)
     and Falush, Stephens and Pritchard (2003)
       Code by Pritchard, Falush and Hubisz
             Version 2.3.4 (Jul 2012)
----------------------------------------------------


Reading file "mainparams_225_603_SITEID_nolab".
datafile is
structure_225INDIV_603SNPS_SITEID_nolab.str
Reading file "extraparams".
Reading file "structure_225INDIV_603SNPS_SITEID_nolab.str".


Data file "structure_225INDIV_603SNPS_SITEID_nolab.str" (truncated) --

Ind:   Label Pop : Genotype_data . . . .
  1: R01-CR-CE-0    1  : -9   1  -9  -9   1  -9   1  . . . .   1
  1: R01-CR-CE-0    1  : -9   1  -9  -9   2  -9   1  . . . .   1
  2: R01-CR-CE-0    1  : -9   1   1   1  -9  -9   1  . . . .   1
  2: R01-CR-CE-0    1  : -9   1   1   1  -9  -9   1  . . . .   1
  3: R01-CR-CE-0    1  : -9   1   1  -9  -9  -9   1  . . . .   1
  3: R01-CR-CE-0    1  : -9   1   1  -9  -9  -9   2  . . . .   1
  4: R01-CR-CE-1    1  : -9   1   1   2  -9  -9  -9  . . . .   1
  4: R01-CR-CE-1    1  : -9   1   1   2  -9  -9  -9  . . . .   1

      *******

765: R06-SK-ML-1   44  :  1   1   1   1   1   1   1  . . . .   1
765: R06-SK-ML-1   44  :  1   1   2   1   1   1   1  . . . .   1
225: R06-SK-ML-1   44  :  1   1   1   1   1   1   1  . . . .  -9
225: R06-SK-ML-1   44  :  1   1   2   1   1   2   1  . . . .  -9

Number of alleles per locus: min= 2; ave=2.0; max= 2


--------------------------------------

Finished initialization; starting MCMC
250000 iterations + 75000 burnin

^C
"


########## testing whether the directory is the problem - it is!!! #######


### works
structure -K 1 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122


# works
structure -K 1 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o /outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122


## works
structure -K 1 -m mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122



# works
structure -K 1 -m /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/mainparams_225_603_SITEID_nolab -e extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122


# works
structure -K 1 -m /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/mainparams_225_603_SITEID_nolab -e /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122




### DOESN'T WORK - ie when add dir for struc file
structure -K 1 -m /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/mainparams_225_603_SITEID_nolab -e /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/extraparams -i /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/structure_225INDIV_603SNPS_SITEID_nolab.str -o /scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_k1r1 -D 37847122




12b. test with joblog

structure -K 1 -m /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/mainparams_225_603_SITEID_nolab -e /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/extraparams -i structure_225INDIV_603SNPS_SITEID_nolab.str -o /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/structure_225INDIV_603SNPS_SITEID_nolab_output_k1r1 -D 94109968 > /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/joblogs_k1to13/job_structure_python_225INDIV_603SNPS_k1to13_screenoutput_k1r1.log


#end job
CRTL + C

cd joblogs_k1to13

# yes there is a file there and there wasn't an error! I can submit the job!

13. SUBMIT JOB

#upload job script and pysthon script
rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/job_structure_python_225_603_k1to13_5reps_nodirforstrfile /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/jobSUB_structure_python_225_603_k1to13_5reps_nodirforstrfile.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/ --progress


## try indiv code

jobSUB_structure_python_nov18_2021parallel_225_603_k1to13_5reps_nodirforstrfile_INDIVCODEk1

rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/jobSUB_structure_python_nov18_2021parallel_225_603_k1to13_5reps_nodirforstrfile_INDIVCODEk* roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/ --progress

rsync -zv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/jobSUB_structure_python_nov18_2021parallel_225_603_k1to13_5reps_nodirforstrfile_INDIVCODEk4.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/ --progress


14. DOWNLOAD OUTPUTS ONCE DONE (3 DAYS TO A WEEK)

## download
rsync -zv roseanna@cedar.computecanada.ca:/scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/outputs_k1to13_5rep/* /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/outputs_k1to13_5rep/  --progress

rsync -zv roseanna@cedar.computecanada.ca:/scratch/roseanna/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_603SNPS/structure/jobSUB_structure_python_225_603_k1to13_5reps_nodirforstrfile_* /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/structure/  --progress


15. now compress outputs_k1to13_5rep folder

16.  upload compressed folder outputs_k1to13_5rep to structure harvester
http://taylor0.biology.ucla.edu/structureHarvester/




17. upload compressed folder outputs_k1to13_5rep to clumpak
http://clumpak.tau.ac.il/index.html

--- EDIT COLOURS

18. save pdf job pipeline summary and zip file output

