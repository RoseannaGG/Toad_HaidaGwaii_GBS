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
