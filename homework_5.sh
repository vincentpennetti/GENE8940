#!/bin/bash
#SBATCH --job-name=homework_5                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=6:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/vjp98982/homework_5"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# # curl the refseq E coli MG1655 genome sequence
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz | gunzip -c > ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa

module load kallisto/0.46.1-foss-2019b
kallisto index -i ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa.idx ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa

kallisto quant -t 6 -b 100 -i ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa.idx -o ${OUTDIR}/kallisto/SRR5344681 /work/gene8940/instructor_data/SRR5344681_1.fastq.gz /work/gene8940/instructor_data/SRR5344681_2.fastq.gz
kallisto quant -t 6 -b 100 -i ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa.idx -o ${OUTDIR}/kallisto/SRR5344682 /work/gene8940/instructor_data/SRR5344682_1.fastq.gz /work/gene8940/instructor_data/SRR5344682_2.fastq.gz
kallisto quant -t 6 -b 100 -i ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa.idx -o ${OUTDIR}/kallisto/SRR5344683 /work/gene8940/instructor_data/SRR5344683_1.fastq.gz /work/gene8940/instructor_data/SRR5344683_2.fastq.gz
kallisto quant -t 6 -b 100 -i ${OUTDIR}/GCF_000005845.2_ASM584v2_cds.fa.idx -o ${OUTDIR}/kallisto/SRR5344684 /work/gene8940/instructor_data/SRR5344684_1.fastq.gz /work/gene8940/instructor_data/SRR5344684_2.fastq.gz


conda activate R
R --no-save < $HOME/GENE8940/homework5.r
