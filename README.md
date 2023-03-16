# Toad Genetics Pipeline  Feb 2023

## Data info
We have GBS double-digested (Sbfl and PstI), paired-end (2x101bp), Illumina Novaseq data with an individual level barcode on the forward read only. The data provided from the sequencing platform came back to us as demultiplexed plates with the Illumina adapter sequences and plate barcodes trimmed off (just leaving the individual barcode and cutsite on the forward read, and just the cutsite on the reverse read). There were three Novaseq lanes. This analysis started with 1370 samples, and 47 breeding sites - from the lower mainland, Vancouver Island, Haida Gwaii and northwest BC. DNA extractions done at UBC Hamelin lab, library prep done at Laval, sequencing done at Genome QC. Bioinformatics were done by RGG on the Compute Canada cluster. Most analyses were done in R. Research funded by National Geographic, BC Ministry of Forests and the UBC Public Scholar Intiative. RGGs stipend provided by Govt. Canada Vanier Scholarship and BC Ministry of Forests. 

Useful resources: https://www.biorxiv.org/content/10.1101/2021.11.02.466953v1.full.pdf

https://groups.google.com/g/stacks-users/

https://catchenlab.life.illinois.edu/stacks/

https://catchenlab.life.illinois.edu/stacks/manual/


Steps below, also see "to do.sh" file

## 1. Create directories on server

mkdir gstacks etc

## 2. Upload and check files on server 

A) upload
   - raw reads of the forward and reverse plates - fastq.gz files 
   - Anaxyrus boreas reference genome .fasta file
   - stacks barcode files 
   - pop.map .tsv  
   - job scripts 
    
      rsync -zv /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/*.txt /drives/f/GBS_data_03_02_21/Lane2_GBSdata_2022/Trim_demultiplex_Oct2022/*.sh roseanna@cedar.computecanada.ca:/home/roseanna/scratch/Demultiplexing_stacks/  --progress
   
   
   B) do md5sum check before downloading raw reads from nanuq and also upon download to laptop and after uploaded to server 


i) on hard drive (using windows commander)

CertUtil -hashfile "I:\Lane_3_GBS_data_13_02_23\NS.2053.003.B716---D502.Hamelin__20230109-Plate-2_R2.fastq.gz" MD5

ii) on server

md5sum NS.2053.003.B716---D502.Hamelin__20230109-Plate-2_R2.fastq.gz


## 3. FastQC plates

using interactive node for one plate, or do to all, create slurm script

salloc --time=01:00:00 --gres=gpu:1 --cpus-per-task=8 --mem=1000M --ntasks=1 --account=def-saitken

Fine except low quality score for reverse read cut site

module load fastqc

fastqc NS.1760.001.B711---D503.Hamelin_202110_plate2_R1.fastq.gz

investigate files
zcat NS.1760.001.B711---D503.Hamelin_202110_plate2_R1.fastq.gz | head -n 20


## 4. Trim reverse reads

using slurm

- Take out the restriction enzyme cut site on reverse read

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 4 -trimlog /scratch/roseanna/Trimmed_reverseplates/D709_P9_R2_trimmed_log.txt /project/def-saitken/roseanna/rawreads_md5checked_feb272021/NS.1470.002.D709.Hamelin_202010_plate__5_R2.fastq.gz /scratch/roseanna/Trimmed_reverseplates/NS.1470.002.D709.Hamelin_202010_plate__5_R2.trimmed.fastq.gz HEADCROP:3



## 5. Demultiplex/clean plates - stacks/process_radtags

using slurm

- discards reads with low quality (below 20) in sliding window of 25% of read length
- make sure stacks_barcode file is unix format
- can't use multiple cpus

process_radtags -1 /home/roseanna/scratch/rawreads_feb4th2023/NS.1760.001.B712---B501.Hamelin_202110_plate4_R1.fastq.gz -2 /scratch/roseanna/Trimmed_reverseplates/NS.1760.001.B712---B501.Hamelin_202110_plate4_R2.trimmed.fastq.gz -b /scratch/roseanna/Demultiplexing_stacks/stacks_barcode_B501.txt -o /scratch/roseanna/Demultiplexing_stacks/B501_norenz2_trimR2 -w 0.25 -s 20 -y gzfastq --inline_null --renz_1 pstI --quality --rescue --barcode-dist-1 1 --threads $cpu -D &> process_radtags_standoutputerror_B501_norenz2_trimR2.oe 

- download process radtags.oe files
-copy all .fq.gz files into one single folder (all the .rem.1, .1., .2, .rem.2 )


## 6. FastQC samples

 module load fastqc

fastqc file

## 7. Align the reads - BWA

using slurm
only keep reads that align with >20 quality score


a)
-defintely use multiple cpus - speeds it up a lot

bwa index -p /scratch/roseanna/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas $genome_fa &> /scratch/roseanna/Anaxyrus_boreas_genome/bwa/bwa_Anaxyrus_boreas_genome_index.oe

-download logs
-check if all bam files have data in them

### b) check alignment using samtools

module load samtools/1.16.1 StdEnv/2020

x) samtools view -q 30 -c R02-SC-IN-01.bam # counts reads mapped with over 30 quality
y) samtools flagstat R02-SC-IN-01.bam -O tsv # gives summary of total numner of reads, % mapped etc. 

Divide (x) by (y) to get % mapped reads over 30 quality

## 8. Build loci - stacks/gstacks

using slurm

gstacks -I $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/alignments_mapq20/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ --min-mapq 20 -t $cpu

-download logs and log distribs 
-see if can backup gstacks files on project server?

## 9. Call snps - stacks/populations

using slurm

- take first snp in locus (locus =  DNA between cut sites)
- run with pop map that does not include experiment pond samples R01-SS-EF - beacuse I have fawn lake in there as well
- run with various dif filtering criteria
- filter max obs het 0.6 across all samples

NB - run steps 7-9 on small test dataset first and play with values of M and n to see which to choose. As per Rochette and Catchan paper

Examples:
echo “Starting populations_ANBOref_r60_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu


echo “Starting populations_ANBOref_r60_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-samples-per-pop 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu

echo “Starting populations_ANBOref_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp at: `date`”

populations -P $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/ -M $src/info/Pop_map_HGthesis_lane1_allexcptexperimnt_lane2_LM_VanIs_lane3_NW_siteID.tsv -O $src/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_R60pctoverall_minmaf0.01_maxhet0.6_writesinglesnp/ --max-obs-het 0.6 --min-samples-overall 0.6 --min-maf 0.01 --write-single-snp --fstats --hwe --smooth --smooth-fstats --vcf --vcf-all --plink --structure --genepop --treemix --verbose -t $cpu



### 9b.  nucleotide diversity (pi)
  - run populations with some good filtering to get that has had some filtering
  - export sum stats summary file
  - manipulate in R to make plot
  - might also make P value? or new version of stacks gets a p value for this
  - include whatever samples (in R) for whatever regions you want
  - download vcfall file for this

## 10. Filtering - O'Leary paper
  ### Filter biallelic, mac=3, min depth =5, mean min depth =5, GQ =20, rm indels
  
  VCF tools
  
  ### Fitler missing SNPs and missing individuals iteratively until got dataset with SNPS missing in less than 10% of indiv, and indivs missing less than 25% of SNPs.

  VCF tools
  
  
  ## 10b. Export vcf file and missing individual files to desktop
  

   ##E filter out high obs het 

    - per region 0.6
   - R
  
    ### Fitler out siblings 0.2 full sibs and half sibs (if just full sibs 0.4 or 0.35)
  
  - plink and R
  - using 0.4 threshold

  ### Fitler SNPs out of HWE 0.01
  
  - Per region, averaged across each locus
  - R
  
    ### Fitler SNPs out of FIS 
  
  - Per pond, averaged across each locus
  - plot histogram to define cut offs - doesn't have to be symetrical 
  - R



  
  ## ANALYSIS ON WHOLE DATASET
  
  NB it is SUPER important to keep the popmap txt file with the samples in the same order as the stacks pipeline (popmap.tsv) used for all the analyses. This means whenever you need to amke changes to the pop map file (e.g. removing individuals) you MUST EDIT it in R and not excel so that you don't risk re-ordering the names.

## 11. Plot PCA of filtered snps

R

## 12. FST comparison 

R 

- between Haida Gwaii and Mainland & Vancouver Island
- within Haida Gwaii

## 13. Isolation-by-distance

R 

- make sure have the correct lat longs - there was a version that didn't have them correct
- make sure use the right dist calc - double check dist matric mataches real world

## 14. Structure

- export structure file from R 
- test structure locally
- upload to server and run
- upload output to structure harvester to get optimal K value 
      http://taylor0.biology.ucla.edu/structureHarvester/
- next could do CLUMPP and DISTRUCT locally to get plots (downlaod from standford)
- OR take strucutre outputs and put them in a webstie to get plots for a range of K values - and then clip those pdfs
       http://clumpak.tau.ac.il/index.html
       
 NB structure runs x10 faster if the directory in the command is written a certain way - see structure instructions

## 15. FIS

R

## 16. Expected het

R

  ## ANALYSIS ON SEPERATE DATASETS - E.G. Haida Gwaii and Vancouver Island
  
  ### Separate datasets in R - i.e. fully filtered dataset at the end of step 10.
  
  - seppop() function adegenet 

## run steps 11-16 for separate datasets (also the nucleotide diversity output from populations)

NB Bad apple 2021 Jose Cerca paper suggests running populations separately (with dif pop maps for region in my case) and then indentify individuals missing a lot of data first, then re-running populations for all regions at once but WITHOUT those bad apples in the pop map

## NE Estimator



