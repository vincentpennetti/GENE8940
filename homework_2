#!/bin/bash
#SBATCH --job-name=homework_2                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=2:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/vjp98982/homework_2"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi


# curl the refseq genome sequence and GFF annotation
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > ${OUTDIR}/GCF_000005845.2_ASM584v2.fna
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz | gunzip -c > ${OUTDIR}/GCF_000005845.2_ASM584v2.gff

# load necessary modules for file processing
module load BEDOPS/2.4.39-foss-2019b
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load ucsc/359

# convert the GFF file to BED format (using convert2bed)
convert2bed --input=gff < ${OUTDIR}/GCF_000005845.2_ASM584v2.gff > ${OUTDIR}/GCF_000005845.2_ASM584v2.bed

# filter BED file to create a new BED file with only CDS regions (using grep)
grep "CDS" ${OUTDIR}/GCF_000005845.2_ASM584v2.bed > ${OUTDIR}/GCF_000005845.2_ASM584v2_CDS.bed

# create a "genome" index file for BEDtools
samtools faidx ${OUTDIR}/GCF_000005845.2_ASM584v2.fna
cut -f1,2 ${OUTDIR}/GCF_000005845.2_ASM584v2.fna.fai > ${OUTDIR}/GCF_000005845.2_ASM584v2.genome.txt

# use the CDS region BED file to create a complementary set of BED intervals for non-CDS regions (using BEDtools complement)
bedtools complement -i ${OUTDIR}/GCF_000005845.2_ASM584v2_CDS.bed -g ${OUTDIR}/GCF_000005845.2_ASM584v2.genome.txt > ${OUTDIR}/GCF_000005845.2_ASM584v2_intergenic.bed


# generates two files of fasta sequences for all CDS and all non-CDS regions respectively (using bedtools getfasta)
bedtools getfasta -fi ${OUTDIR}/GCF_000005845.2_ASM584v2.fna -bed ${OUTDIR}/GCF_000005845.2_ASM584v2_CDS.bed -fo ${OUTDIR}/GCF_000005845.2_ASM584v2_extractedCDS.fna
bedtools getfasta -fi ${OUTDIR}/GCF_000005845.2_ASM584v2.fna -bed ${OUTDIR}/GCF_000005845.2_ASM584v2_intergenic.bed -fo ${OUTDIR}/GCF_000005845.2_ASM584v2_extractedIntragenic.fna

# compute the GC content of  CDS and non-CDS intergenic regions (using faCount -summary)
faCount -summary ${OUTDIR}/GCF_000005845.2_ASM584v2_extractedCDS.fna
faCount -summary ${OUTDIR}/GCF_000005845.2_ASM584v2_extractedIntragenic.fna



# cut out the third column (specifying the nucleotide type), pipe to grep
# each line containing "CDS" is piped to wc
# wc counts the number of lines output from grep and stores the total in results.txt
cut -f3 ${OUTDIR}/ecoli_MG1655.gff | grep "CDS" | wc -l > ${OUTDIR}/results.txt
