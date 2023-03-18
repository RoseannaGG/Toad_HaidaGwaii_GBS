FastQC plates

using interactive node for one plate, or do to all, create slurm script

salloc --time=01:00:00 --gres=gpu:1 --cpus-per-task=8 --mem=1000M --ntasks=1 --account=def-saitken

Fine except low quality score for reverse read cut site

module load fastqc

fastqc NS.1760.001.B711---D503.Hamelin_202110_plate2_R1.fastq.gz

investigate files
zcat NS.1760.001.B711---D503.Hamelin_202110_plate2_R1.fastq.gz | head -n 20