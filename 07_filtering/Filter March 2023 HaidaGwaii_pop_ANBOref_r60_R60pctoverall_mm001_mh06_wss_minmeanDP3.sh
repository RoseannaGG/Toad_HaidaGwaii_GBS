


########## beluga ###########


# http://www.ddocent.com/filtering/
salloc --time=01:00:00 --cpus-per-task=1 --mem=1000M --ntasks=1 --account=def-saitken




module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/
cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss




#############################################################
#############################  feb 24th 2023#########################

########## remove R02-SS-EF samples #############


############ EXPORT unfiltered data #########################


#rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/gapped0.9.recode.vcf /drives/f/GBS_data_03_02_21/De_novo/stacks_gapped0.9.M3/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/ --progress


module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/
cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss


### per site snp quality ###

vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --site-quality --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/sitequality


#################################################################################
## 1. prelim filtering before removing sibs


 #compute imissing
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/


"
"

######### prelim fitlering

## mean min depth 3 - 
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/populations.snps.vcf --remove-indels --min-meanDP 3 --minDP 3 --min-alleles 2 --max-alleles 2 --minGQ 20 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic

"




After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 17529 out of a possible 53334 Sites
Run Time = 12.00 seconds









"




 #compute imissing
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.listindivmissingSNPS

"
"

mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

                                    Histogram of % missing data per individual                      
                                                                                                       
  10 +-+--------+----------+----------+----------+--**-------+----------+----------+----------+----------+--------+-+
     +          +          +          +          +  **       +          +          +          +          +          +
     |                                              **                                                      ******* |
   9 +-+                               ***          **                                                            +-+
     |                                 * *          **                                                              |
   8 +-+     **      ***             *** *          **     ***                                                    +-+
     |       **      ***             *** *          **     * *                                                      |
     |       **      ***             *** *          **     * *                                                      |
   7 +-+     **      *** ** **  **   *** *    ***   ***  *** *                                    **              +-+
     |       **      *** ** **  **   *** *    ***   ***  *** *                                    **                |
   6 +-+     ****    *** *********   *** *    ***   **** *** **    **         *** **         ***  **              +-+
     |       ****    *** ***** ***   *** *    ***   **** *** **    **         * * **         * *  **                |
     |       ****    *** ***** ***   *** *    ***   **** *** **    **         * * **         * *  **                |
   5 +-+     ******* *** ***** ***   *** * ** ********** *** **  ******       * * **         * *  *****  **       +-+
     |       ******* *** ***** ***   *** * ** **** ***** *** **  ******       * * **         * *  *** *  **         |
   4 +-+ *** ******* ********* ********* * ******* ********* **  ******       * * **        ** *  *** ** ****     +-+
     |   * * ******* ********* ********* * ******* ********* **  ******       * * **        ** *  *** ** ****       |
     |   * * ******* ********* ********* * ******* ********* **  ******       * * **        ** *  *** ** ****       |
   3 +-+ * * ***************** ********* ********* ****************************************************************-+
     |   * * ******* ********* ********* ********* ************* ******  ****** * **  ***** ** ****** *******     * |
   2 +-+ * * ******* ********* ********* ********* ********************  ****** *********** ** ****** ******* *** *-+
     |   * * ******* ********* ********* ********* ********************  ****** * ********* ** ****** ******* * * * |
     +   * * ******* ********* ********* ********* ******************** +****** * ********* **+****** ******* * * * +
   1 +-+-***-******************************************************************************************************-+
     0         0.1        0.2        0.3        0.4         0.5        0.6        0.7        0.8        0.9         1
                                                    % of missing data                                  
                          

# 80%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 20% of the samples
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP0.8



"

"


# 70%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.7 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP70%

"
"

##############################
####### SNP > 40% #########
################################

# 60%
# remove snps with over 40% missing data across all individuals
#exclude snps that are not present in at least 40% of the samples
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic.recode.vcf --max-missing 0.6 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%


"

After filtering, kept 385 out of 385 Individuals
Outputting VCF file...
After filtering, kept 3433 out of a possible 17529 Sites
Run Time = 3.00 seconds


"



 #compute imissing - * are indiv
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS

"
"

mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


                                       Histogram of % missing data per individual                      
                                                                                                       
  18 +-+-------+---------+---------+---------+---------+----------+---------+---------+---------+---------+-------+-+
     +    **   +         +         +         +         +          +         +         +         +         +         +
     |    **                                                                                                ******* |
  16 +-+  **                                                                                                      +-+
     |    **                                                                                                        |
  14 +-+  **                                                                                                      +-+
     |    ****    **          **                                                                                    |
     |    ****    **          **                                                                                    |
  12 +-+  ******  **          **                                                                                  +-+
     |    ******  **          ****                                                                                  |
  10 +-+  **********          ****                                                                                +-+
     |    **********     **   **** **                                                                               |
     |    **********     **   **** **                                                                               |
   8 +-+  **********     **   *******                                                                             +-+
     |   ***********     **   *******                                                                               |
   6 +-+ ******************** ********** **   **                                                                  +-+
     |   ******************** ********** **   **                                                                    |
     |  ********************* ********** **   **                                                                    |
   4 +-+********************* ******************** ***   ***  ***** ***            **       **                    +-+
     |  ********************* ************************   ***  *********     **     ************                     |
   2 +-+**************************************************** **********************************    ***            +-+
     |  **************************************************** **********************************    * *              |
     +  ***************************************************************************************************         +
   0 +-+***************************************************************************************************-------+-+
     0        0.1       0.2       0.3       0.4       0.5        0.6       0.7       0.8       0.9        1        1.1
                                                    % of missing data                                  
                                                                                                       


## 30%  ## removing individuals with lss than 70% call rate, or that have more than 30% missing data
mawk '$5 > 0.3' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%

"




"
##############################
####### INDIV > 40% #########
################################

## 40%  ## remove 369
mawk '$5 > 0.4' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%

"


Excluding individuals in 'exclude' list
After filtering, kept 252 out of 385 Individuals
Outputting VCF file...
After filtering, kept 3433 out of a possible 3433 Sites
Run Time = 4.00 seconds






"

 #compute imissing site (snp) - % missing indiv on a per SNP basis * is snp
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%.recode.vcf --missing-site --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing



""




mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing.lmiss | cut -f6 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listSITEmissing_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

  

# 70%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%.recode.vcf --max-missing 0.7 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%


"

"



##############################
####### SNP > 20% #########
################################


# 80%
# remove snps with over 30% missing data across all individuals
#exclude snps that are not present in at least 30% of the samples
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%.recode.vcf --max-missing 0.8 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%


"

After filtering, kept 252 out of 252 Individuals
Outputting VCF file...
After filtering, kept 2008 out of a possible 3433 Sites
Run Time = 1.00 seconds


"



 #compute imissing - * are indiv
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS


mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


                                       Histogram of % missing data per individual                      
                                                                                                       
  25 +-+-----------+-------------+-------------+-------------+------------+-------------+-------------+-----------+-+
     +             +             +             +             +            +             +             +             +
     |                                                                                                      ******* |
     |                                                                                                              |
     |       ****                                                                                                   |
  20 +-+     *  *                                                                                                 +-+
     |       *  *          ****                                                                                     |
     |       *  *          *  *                                                                                     |
     |       *  ****       *  *                                                                                     |
  15 +-+     *  *  *       *  *                                                                                   +-+
     |       *  *  *       *  *                                                                                     |
     |       *  *  *       *  *                                                                                     |
     |       *  *  *       *  *                                                                                     |
     |       *  *  *       *  ****                                                                                  |
  10 +-+     *  *  ****    *  *  *    ****                                                                        +-+
     |     ***  *  *  *    *  *  *    *  *                     ****             ****                                |
     |     * *  *  *  ***  *  *  **** *  *                     *  *********     *  *                                |
     |     * *  *  *  * ****  *  *  * *  ****  ***           ***  *  *  * *******  * ****                           |
   5 +-+   * *  *  *  * *  *  *  *  * *  *  **** *  ****     * *  *  *  * *  *  *  * *  ****  *****   ****        +-+
     |     * *  *  *  * *  *  *  *  * *  *  *  * ****  ****  * *  *  *  * *  *  *  * *  *  ****   *   *  *          |
     |     * *  *  *  * *  *  *  *  ***  *  *  * *  *  * ****************************************************       |
     |     * *  *  *  * *  *  *  *  * *  *  *  * *  *  * **  * *  *  *  * *  *  *  * *  *  *  *   *****  *  *       |
     +  **** *  *  *  * *  *  *  *  * *  *  *  * *  *  * **  * *  *  *  * *  *  *  * *  *  *  *   *   *  *  *       +
   0 +-+*****************************************************************************************************-----+-+
     0            0.05          0.1           0.15          0.2          0.25          0.3           0.35          0.4
                                                    % of missing data                                  

#allele freq per loci
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --freq2 --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq


pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq.frq




mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq.frq | cut -f6 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listSNPminorallelefreq_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


### MAF
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_maf001


"

"




 #compute imissing - * are indiv
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS




## 40%  ## 
mawk '$5 > 0.4' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv

### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_40%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%



""
##############################
####### INDIV > 30% #########
################################


 #compute imissing - * are indiv
#vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV40%_listindivmissingSNPS

## 30%
mawk '$5 > 0.3' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv

### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%

"

Excluding individuals in 'exclude' list
After filtering, kept 225 out of 252 Individuals
Outputting VCF file...
After filtering, kept 2008 out of a possible 2008 Sites
Run Time = 1.00 seconds



"



## 70% snpmissing

#compute missingness
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNPS

## 30%
mawk '$5 > 0.3' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_30%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_30%missINDV.indv

### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_30%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_rmINDIV30%


"
"


# maf 0.01
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_rmINDIV30%_maf001




 #compute imissing - * are indiv
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --missing-indv --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS



"


"


mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS.imiss | cut -f5 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listindivmissingSNPS_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF



 #compute imissing site (snp) - % missing indiv on a per SNP basis * is snp
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --missing-site --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing.lmiss




mawk '!/IN/' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing.lmiss | cut -f6 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing_totalmissing
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
plot 'pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_listSITEmissing_totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

 


## 60%  ## remove 369
mawk '$5 > 0.6' pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNPS.imiss | cut -f1 > pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv

cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv



### ACTUALLY REMOVE INDIVI
vcftools --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%.recode.vcf --remove  pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_listindivmissingSNP_60%missINDV.indv --recode --recode-INFO-all --out pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP70%_rmINDIV60%


""




### MAF
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV30%_maf001


"






"

# maf 0.01
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --maf 0.01 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf001


"




"

##############################
####### MAF 0.03 #########
################################


# maf 0.03
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%.recode.vcf --maf 0.03 --remove-filtered-all --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003


"


After filtering, kept 225 out of 225 Individuals
Outputting VCF file...
After filtering, kept 723 out of a possible 2008 Sites
Run Time = 1.00 seconds





"






### singletons gapped0.9_maxhet0.6_singleSNP_mac3_DP5_minGQ20_maf0.05_rmindels_biallelic_rmsibs_missingSNP0.6
vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --singletons --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003



"

"

#no singletons
cat pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.singletons
CHROM   POS     SINGLETON/DOUBLETON     ALLELE  INDV




#excludes singletons
#vcftools --vcf $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --remove-filtered-all --exclude-positions pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.singletons --recode --recode-INFO-all --out $src/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003

"

"






########### EXPORT  ###########


# 225INDIV_723SNPS
### pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003



rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_listindivmissingSNP_40%missINDV.indv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_listindivmissingSNP_30%missINDV.indv /drives/f/GBS_data_03_02_21/HaidaGwaii_ref_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/225INDIV_723SNPS/ --progress














######### CONVERT TO PLINK ######
salloc --time=01:00:00 --cpus-per-task=1 --mem=1000M --ntasks=1 --account=def-saitken




module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/
cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss

module load plink/1.9b_6.21-x86_64 StdEnv/2020

plink --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --allow-extra-chr --recode --out 766INDIV_3496SNPS

"[roseanna@beluga5 pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss]$ plink --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf --allow-extra-chr --recode --out 766INDIV_3496SNPS
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 766INDIV_3496SNPS.log.
Options in effect:
  --allow-extra-chr
  --out 766INDIV_3496SNPS
  --recode
  --vcf pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss_singleSNP_minDP3_minGQ20_rmindels_biallelic_missingSNP60%_rmINDIV40%_missingSNP80%_rmINDIV30%_maf003.recode.vcf

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

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ponds_kept_766INDIV_3496SNPS.txt roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/ --progress



###########################################################################
####### remove sibs from filtered vcf file #############
#############################################################

module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.14
src=/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/
cd /home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/766INDIV_1430SNPS_editlocinamesinR.vcf roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress

rsync -zv /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/Rmsibs/listindiv_plink_sub_sib_rel_over_04_names_unique.indv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress


# remove

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

rsync -zv roseanna@beluga.computecanada.ca:/home/roseanna/scratch/HaidaGwaii_Aboreas/gstacks_minmapq20/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/510INDIV_1430SNPS_editlocinamesinR_rmsibsvcftools.recode.vcf /drives/f/GBS_data_03_02_21/Lane_1_2_3_feb2023/gstacks_minmapq20_1370/pop_ANBOref_r60_R60pctoverall_mm001_mh06_wss/766INDIV_3496SNPS/ --progress