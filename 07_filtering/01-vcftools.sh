


########## CEDAR ###########


# http://www.ddocent.com/filtering/
salloc --time=01:00:00 --cpus-per-task=1 --mem=1000M --ntasks=1 --account=def-saitken




module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/
cd /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss




#############################################################
#############################  feb 24th 2023#########################

########## remove R02-SS-EF samples #############


############ EXPORT unfiltered data #########################


#rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/gapped0.9.recode.vcf /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/ --progress




### per site snp quality ###

vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --site-quality --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/sitequality


#################################################################################
## 1. prelim filtering before removing sibs
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/
cd /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss


 #compute imissing
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/


"
"

######### prelim fitlering
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --remove-indels --minDP 3 --min-alleles 2 --max-alleles 2 --minGQ 20 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic

"




After filtering, kept 1370 out of 1370 Individuals
Outputting VCF file...
After filtering, kept 165724 out of a possible 165724 Sites
Run Time = 516.00 seconds






"



##############################
####### INDIV > 30% #########
################################


 #compute imissing
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS

"
"

mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS.imiss | cut -f5 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF



                                       Histogram of % missing data per individual                      
                                                                                                       
  70 +-+-------------+---------------+---------------+--------------+---------------+---------------+-------------+-+
     +               +               +               +              +               +               +               +
     |                                                                                                      ******* |
  60 +-+                                         ***********************************                              +-+
     |                                           *                                 *                                |
     |                                           *                                 *                                |
     |                                           *                                 *                                |
  50 +-+                                         *                                 *                              +-+
     |                                           *                                 *                                |
     |                                           *                                 *                                |
  40 +-+                                         *                                 *                              +-+
     |                                           *                                 *                                |
     |                                           *                                 *                                |
  30 +-+                                         *                                 *                              +-+
     |                                           *                                 *                                |
     |                           **              *                                 *                                |
  20 +-+          *              **   **         *                              ** *                              +-+
     |     **    **  **          **   **   ***  **   *       *        **   **   ** *                                |
     |  *****************   ********* ****************   *******  **  ******** *****                                |
     |  ****************** ***************************************** ********* *****                                |
  10 +-+****************************************************************************                              +-+
     | *****************************************************************************                                |
     + *****************************************************************************+               +               +
   0 +-*****************************************************************************+---------------+-------------+-+
     0              0.2             0.4             0.6            0.8              1              1.2             1.4
                                                    % of missing data                                  
                            

# 80%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 20% of the samples
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP0.8



"
"


# 70%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.7 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP70%

"
"



# 60%
# remove snps with over 40% missing data across all individuals
#exclude snps that are not present in at least 40% of the samples
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.6 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%


"
After filtering, kept 1370 out of 1370 Individuals
Outputting VCF file...
After filtering, kept 10133 out of a possible 165724 Sites
Run Time = 61.00 seconds

"



 #compute imissing - * are indiv
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS

"
"

mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f5 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


                                      Histogram of % missing data per individual                      
                                                                                                       
  80 +-+-------------+---------------+---------------+--------------+---------------+---------------+-------------+-+
     +    **         +               +               +              +               +               +               +
     |    **                                                                                                ******* |
  70 +-+  **                                                                                                      +-+
     |    **                                                                                                        |
     |    **                                                                                                        |
  60 +-+  **                                                                                                      +-+
     |    **                                                                                                        |
     |    **                                                                                                        |
  50 +-+  ****                                                                     **                             +-+
     |    ****                                                                     **                               |
  40 +-+  ****                                                                     **                             +-+
     |    ****                                                                     **                               |
     |    **** **                                                                  **                               |
  30 +-+  *******                                                                  **                             +-+
     |    ********        **                                                       **                               |
     |    *********  **   **                                                       **                               |
  20 +-+  *************   ****                                                     **                             +-+
     |   *********************          **                                         **                               |
     | *************************  ***  ***          **                             **                               |
  10 +-**************************************** ****** ***  ******   **  **  **    **                             +-+
     |***************************************************** *********************  **                               |
     +*******************************************************************************               +               +
   0 +*******************************************************************************---------------+-------------+-+
     0              0.2             0.4             0.6            0.8              1              1.2             1.4
                                                    % of missing data                                  





## 30%  ## removing individuals with lss than 70% call rate, or that have more than 30% missing data
mawk '$5 > 0.3' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f1 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --remove  populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%

"



Excluding individuals in 'exclude' list
After filtering, kept 762 out of 1370 Individuals
Outputting VCF file...
After filtering, kept 10133 out of a possible 10133 Sites
Run Time = 20.00 seconds


"


## 40%  ## remove 369
mawk '$5 > 0.4' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f1 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --remove  populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv --recode --recode-INFO-all --out populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%

"

Excluding individuals in 'exclude' list
After filtering, kept 883 out of 1370 Individuals
Outputting VCF file...
After filtering, kept 10133 out of a possible 10133 Sites
Run Time = 23.00 seconds


"

 #compute imissing site (snp) - % missing indiv on a per SNP basis * is snp
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --missing-site --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing



"After filtering, kept 1370 out of 1370 Individuals
Outputting Site Missingness
After filtering, kept 10133 out of a possible 10133 Sites
Run Time = 2.00 seconds
"




mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing.lmiss | cut -f6 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing INDV per SNP"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


                                           Histogram of % missing INDV per SNP                        
                                                                                                       
  2000 +-+-----------+------------+-------------+-------------+------------+-------------+------------+-----------+-+
       +             +            +             +             +            +             +          ***             +
  1800 +-+                                                                                          * *     *******-+
       |                                                                                            * *             |
       |                                                                                         **** *             |
  1600 +-+                                                                                       *  * *           +-+
       |                                                                                         *  * *             |
  1400 +-+                                                                                       *  * *           +-+
       |                                                                                      ****  * *             |
  1200 +-+                                                                                    *  *  * *           +-+
       |                                                                                      *  *  * *             |
  1000 +-+                                                                                    *  *  * *           +-+
       |                                                                                   ****  *  * *             |
       |                                                                                   *  *  *  * *             |
   800 +-+                                                                               ***  *  *  * *           +-+
       |                                                                                 * *  *  *  * *             |
   600 +-+                                                                            **** *  *  *  * *           +-+
       |                                                                           ****  * *  *  *  * *             |
   400 +-+                                                                       ***  *  * *  *  *  * *           +-+
       |                                                                      **** *  *  * *  *  *  * *             |
       |                                                                *******  * *  *  * *  *  *  * *             |
   200 +-+                                                         ************************************           +-+
       +             +            +             +  *****************  * *  *  *  * *  *  * *  *  *  * *             +
     0 +-+*********************************************************************************************-----------+-+
      0.05          0.1          0.15          0.2           0.25         0.3           0.35         0.4           0.45
                                                     % of missing data                                 




# 70%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%.recode.vcf --max-missing 0.7 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%


"
After filtering, kept 883 out of 883 Individuals
Outputting VCF file...
After filtering, kept 10074 out of a possible 10133 Sites
Run Time = 20.00 seconds
"


# 80%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%.recode.vcf --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%


"
After filtering, kept 883 out of 883 Individuals
Outputting VCF file...
After filtering, kept 8450 out of a possible 10133 Sites
Run Time = 31.00 seconds

"



 #compute imissing - * are indiv
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS


mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f5 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

 
                                       Histogram of % missing data per individual                      
                                                                                                       
  90 +-+---------+------------+-----------+-----------+------------+-----------+-----------+------------+---------+-+
     +         ***            +           +           +            +           +           +            +           +
     |         * *                                                                                          ******* |
  80 +-+       * *                                                                                                +-+
     |         * *                                                                                                  |
  70 +-+       * *                                                                                                +-+
     |         * *                                                                                                  |
     |         * ****                                                                                               |
  60 +-+       * *  *                                                                                             +-+
     |         * *  *                                                                                               |
  50 +-+       * *  *                                                                                             +-+
     |         * *  *                                                                                               |
     |         * *  *                                                                                               |
  40 +-+       * *  * ******                                                                                      +-+
     |         * *  ***  * ****    ***                                                                              |
  30 +-+       * *  * *  * *  * **** *      ****                                                                  +-+
     |      **** *  * *  * *  * *  * *  *** *  *                                                                    |
     |      *  * *  * *  * *  * *  * *  * * *  *      ****              ***                                         |
  20 +-+  ***  * *  * *  * *  * *  * *  * * *  * ******  *    ***       * *  ***                                  +-+
     | **** *  * *  * *  * *  ***  * **** * *  * *  * *  * **** ********* *  * *    ***            ***              |
  10 +-*  * *  * *  * *  * *  * *  * *  * ***  ***  * *  ***  * *  * *  * **** ***  * ******       * ****         +-+
     | *  * *  * *  * *  * *  * *  * *  * * *  * *  * *  * *  * *  * *  * *  * * **** *  * ********* *  *           |
     + *  * *  * *  * *  * *  * *  * *  * * *  * *  * ***************************************************           +
   0 ****************************************************************************************************---------+-+
     0          0.05         0.1         0.15        0.2          0.25        0.3         0.35         0.4         0.45
                                                    % of missing data                                  
                                                                                                       



#allele freq per loci
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --freq2 --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq


populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq.frq




mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq.frq | cut -f6 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


 
                                        Histogram of % missing data per individual                     
                                                                                                       
  3000 +-+---------+-----------+-----------+-----------+------------+-----------+-----------+-----------+---------+-+
       +           +           +           +           +            +           +           +           +           +
       | ****                                                                                               ******* |
       | *  *                                                                                                       |
  2500 +-*  *                                                                                                     +-+
       | *  *                                                                                                       |
       | *  *                                                                                                       |
       | *  *                                                                                                       |
  2000 +-*  *                                                                                                     +-+
       | *  *                                                                                                       |
       | *  *                                                                                                       |
  1500 +-*  ***                                                                                                   +-+
       | *  * *                                                                                                     |
       | *  * *                                                                                                     |
       | *  * *                                                                                                     |
  1000 +-*  * *                                                                                                   +-+
       | *  * ****                                                                                                  |
       | *  * *  *                                                                                                  |
       ***  * *  ***                                                                                                |
   500 *-*  * *  * *                                                                                              +-+
       * *  * *  * ****                                                                                             |
       * *  * *  * *  *****                                                                                         |
       * *  * *  * *  * * ***********      +           +            +           +           +           +           +
     0 ****************************************************************************************************-------+-+
       0          0.05        0.1         0.15        0.2          0.25        0.3         0.35        0.4         0.45
                                                     % of missing data                                 
                     


### MAF
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_maf001


"
After filtering, kept 883 out of 883 Individuals
Outputting VCF file...
After filtering, kept 7794 out of a possible 8450 Sites
Run Time = 17.00 seconds
"




 #compute imissing - * are indiv
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS




## 40%  ## 
mawk '$5 > 0.4' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f1 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv

### ACTUALLY REMOVE INDIVI
vcftools --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --remove  populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv --recode --recode-INFO-all --out populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%



"Excluding individuals in 'exclude' list
After filtering, kept 877 out of 883 Individuals
Outputting VCF file...
After filtering, kept 8450 out of a possible 8450 Sites
Run Time = 33.00 seconds
"



 #compute imissing - * are indiv
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%_listindivmissingSNPS

## 30%
mawk '$5 > 0.3' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f1 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv

### ACTUALLY REMOVE INDIVI
vcftools --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --remove  populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%

"Excluding individuals in 'exclude' list
After filtering, kept 766 out of 883 Individuals
Outputting VCF file...
After filtering, kept 8450 out of a possible 8450 Sites
Run Time = 28.00 seconds
"


 #compute imissing - * are indiv
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --missing-indv --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS


"After filtering, kept 766 out of 766 Individuals
Outputting Individual Missingness
After filtering, kept 8450 out of a possible 8450 Sites
Run Time = 2.00 seconds
"


mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS.imiss | cut -f5 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF



                                       Histogram of % missing data per individual                      
                                                                                                       
  90 +-+---------------+------------------+------------------+-----------------+-----------------+----------------+-+
     +              ****                  +                  +                 +                 +                  +
     |              *  *                                                                                    ******* |
  80 +-+            *  *                                                                                          +-+
     |              *  *                                                                                            |
  70 +-+            *  *                                                                                          +-+
     |              *  *                                                                                            |
     |              *  *****                                                                                        |
  60 +-+            *  *   *                                                                                      +-+
     |              *  *   *                                                                                        |
  50 +-+            *  *   *                                                                                      +-+
     |              *  *   *                                                                                        |
     |              *  *   *                                                                                        |
  40 +-+            *  *   *   ********                                                                           +-+
     |              *  *   *****   *  *****      *****                                                              |
  30 +-+            *  *   *   *   *  *   *   ****   *          *****                                             +-+
     |          *****  *   *   *   *  *   *   *  *   *   *****  *   *                                               |
     |          *   *  *   *   *   *  *   *   *  *   *   *   *  *   *          *****                     *****      |
  20 +-+    *****   *  *   *   *   *  *   *   *  *   *   *  ***************************************************** +-+
     |   ****   *   *  *   *   *   *  *   *****  *   *****  **  *   *   *  *   *   *  *****   ************   *  *   |
  10 +-+ *  *   *   *  *   *   *   *  *   *   *  *   *   *  *****   *****  *   *   ****   *   *  *   *   *   **** +-+
     |   *  *   *   *  *   *   *   *  *   *   *  *   *   *  **  *   *   *  *   *   *  *   *   *  *   *   *   *  *   |
     +   *  *   *   *  *   *   *   *  *   *   *  *   *   *  **  *   *   *  *   *   *  *   *   *  *   *   *   *  *   +
   0 ************************************************************************************************************-+-+
     0                0.05               0.1                0.15              0.2               0.25               0.3
                                                    % of missing data                                  
                                                            


 #compute imissing site (snp) - % missing indiv on a per SNP basis * is snp
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --missing-site --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing.lmiss




mawk '!/IN/' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing.lmiss | cut -f6 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing_totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing INDV per SNP"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF



                                            Histogram of % missing INDV per SNP                        
                                                                                                       
  1600 +-+-------------------+---------------------+--------------------+---------------------+-------------------+-+
       +                     +                     +                    +                     +                     +
       |                                                                                                    ******* |
  1400 +-+                                                      *****                                             +-+
       |                                                        *   *                                               |
       |                                                        *   *                                               |
  1200 +-+                                                 ******   *****                                         +-+
       |                                                   *    *   *   *                                           |
       |                                                   *    *   *   *                                           |
  1000 +-+                                                 *    *   *   ******                                    +-+
       |                                                   *    *   *   *    *                                      |
   800 +-+                                             *****    *   *   *    *                                    +-+
       |                                               *   *    *   *   *    *                                      |
       |                                           *****   *    *   *   *    *                                      |
   600 +-+                                         *   *   *    *   *   *    *                                    +-+
       |                                           *   *   *    *   *   *    *****                                  |
       |                                      ******   *   *    *   *   *    *   *                                  |
   400 +-+                                    *    *   *   *    *   *   *    *   *                                +-+
       |                                      *    *   *   *    *   *   *    *   *                                  |
       |                                      *    *   *   *    *   *   *    *   *                                  |
   200 +-+                            *********    *   *   *    *   *   *    *   *****                            +-+
       |                         ******   *   *    *   *   *    *   *   *    *   *   ******                         |
       +            **************    *   *   *    *   *   *    *   *   *    *   *   *    *****                     +
     0 *********************************************************************************************--------------+-+
       0                    0.05                  0.1                  0.15                  0.2                   0.25
                                                     % of missing data                                 




## 60%  ## remove 369
mawk '$5 > 0.6' populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNPS.imiss | cut -f1 > populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv

cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%.recode.vcf --remove  populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv --recode --recode-INFO-all --out populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_rmINDIV60%


"Excluding individuals in 'exclude' list
After filtering, kept 883 out of 883 Individuals
Outputting VCF file...
After filtering, kept 10074 out of a possible 10074 Sites
Run Time = 20.00 seconds
"




### MAF
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001


"


After filtering, kept 762 out of 762 Individuals
Outputting VCF file...
After filtering, kept 9016 out of a possible 10133 Sites
Run Time = 16.00 seconds



"

# maf 0.01
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf001


"
After filtering, kept 766 out of 766 Individuals
Outputting VCF file...
After filtering, kept 7673 out of a possible 8450 Sites
Run Time = 26.00 seconds

"


# maf 0.03
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --maf 0.03 --remove-filtered-all --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003


"

After filtering, kept 766 out of 766 Individuals
Outputting VCF file...
After filtering, kept 3496 out of a possible 8450 Sites
Run Time = 13.00 seconds


"






### singletons gapped0.9_maxhet0.6_singleSNP_mac3_DP3_minGQ20_maf0.05_rmindels_biallelic_rmsibs_missingSNP0.6
vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --singletons --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003



"
After filtering, kept 766 out of 766 Individuals
Outputting Singleton Locations
After filtering, kept 3496 out of a possible 3496 Sites
Run Time = 0.00 seconds
"

#no singletons
cat populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.singletons
CHROM   POS     SINGLETON/DOUBLETON     ALLELE  INDV




#excludes singletons
#vcftools --vcf $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --remove-filtered-all --exclude-positions populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.singletons --recode --recode-INFO-all --out $src/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003

"

"






########### EXPORT  ###########

### populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001 ###


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001.recode.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/371INDIV_2392SNPS/ --progress



rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.log roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.sumstats_summary.tsv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.log.distribs roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/ --progress


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/catalog.calls roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/catalog.fa.gz /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/ --progress



# 883INDIV_7794SNPS
### populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_maf001


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_maf001.recode.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/883INDIV_7794SNPS/ --progress


# 766INDIV_7673SNPS
### populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf001


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf001.recode.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_7673SNPS/ --progress




# 766INDIV_3496SNPS
### populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf001


rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress












######### CONVERT TO PLINK ######
salloc --time=01:00:00 --cpus-per-task=1 --mem=1000M --ntasks=1 --account=def-saitken




module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/
cd /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss

module load plink/1.9b_6.21-x86_64 StdEnv/2020

plink --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --allow-extra-chr --recode --out 766INDIV_3496SNPS

"[roseanna@cedar5 populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss]$ plink --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --allow-extra-chr --recode --out 766INDIV_3496SNPS
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 766INDIV_3496SNPS.log.
Options in effect:
  --allow-extra-chr
  --out 766INDIV_3496SNPS
  --recode
  --vcf populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_DP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf

257860 MB RAM detected; reserving 128930 MB for main workspace.
Allocated 7259 MB successfully, after larger attempt(s) failed.
--vcf: 766INDIV_3496SNPS-temporary.bed + 766INDIV_3496SNPS-temporary.bim +
766INDIV_3496SNPS-temporary.fam written.
3496 variants loaded from .bim file.
766 people (0 males, 0 females, 766 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 766INDIV_3496SNPS.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 766 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.869599.
3496 variants and 766 people pass filters and QC.
Note: No phenotypes present.
--recode ped to 766INDIV_3496SNPS.ped + 766INDIV_3496SNPS.map ... done.
"

# upload list of ponds kept

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ponds_kept_766INDIV_3496SNPS.txt roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/ --progress



###########################################################################
####### remove sibs from filtered vcf file #############
#############################################################

module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/
cd /home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_1430SNPS_editlocinamesinR.vcf roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/listindiv_plink_sub_sib_rel_over_04_names_unique.indv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress


# remove
vcftools --vcf 766INDIV_1430SNPS_editlocinamesinR.vcf --remove listindiv_plink_sub_sib_rel_over_04_names_unique.indv --out 510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools

vcftools --vcf 766INDIV_1430SNPS_editlocinamesinR.vcf --remove listindiv_plink_sub_sib_rel_over_04_names_unique.indv --recode --recode-INFO-all --out 510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools


"VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf 766INDIV_1430SNPS_editlocinamesinR.vcf
        --remove listindiv_plink_sub_sib_rel_over_04_names_unique.indv
        --out 510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools.vcf

Excluding individuals in 'exclude' list
After filtering, kept 510 out of 766 Individuals
After filtering, kept 1430 out of a possible 1430 Sites
Run Time = 5.00 seconds
"

rsync -zv roseanna@cedar.computecanada.ca:/home/roseanna/scratch/ANBO_refassembly_HGthesis_lane1_lane2_lane3/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools.recode.vcf /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/populations_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress