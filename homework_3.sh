#!/bin/bash
#SBATCH --job-name=homework_3                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=6:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/vjp98982/homework_3"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load necessary modules for file processing
module load canu/1.9-GCCcore-8.3.0-Java-11
module load SPAdes/3.14.1-GCC-8.3.0-Python-3.7.4
module load QUAST/5.0.2-foss-2019b-Python-3.7.4
module load MUMmer/3.23_conda

# curl the refseq genome sequence
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > ${OUTDIR}/GCF_000005845.2_ASM584v2.fna

# download pacbio long reads
wget -O ${OUTDIR}/pacbio.fq.gz http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/ecoli_p6_25x.filtered.fastq.gz

# dowload illumina short reads
wget -O ${OUTDIR}/illumina_read1.fq.gz http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_1.fastq.gz
wget -O ${OUTDIR}/illumina_read2.fq.gz http://bergmanlab.genetics.uga.edu/data/downloads/gene8940/s_6_2.fastq.gz

# assemble the e coli MG1655 genome using PacBio long reads with Canu
gunzip -c ${OUTDIR}/pacbio.fq.gz > ${OUTDIR}/pacbio.fq
canu -p canuPacBio -d ${OUTDIR} genomeSize=4.8m useGrid=false -pacbio-raw ${OUTDIR}/pacbio.fq

# assemble the E coli MG1655 genome using Illumina short read data with SPADES
spades.py -t 6 -k 21,33,55,77 --isolate --memory 24 --pe1-1 ${OUTDIR}/illumina_read1.fq.gz --pe1-2 ${OUTDIR}/illumina_read2.fq.gz -o ${OUTDIR}

# run QUAST to generate assembly quality assessment statistics
quast.py -o $OUTDIR -t 6 -r ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/canuPacBio.contigs.fasta ${OUTDIR}/scaffolds.fasta

# mummer for pacbio
nucmer ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/canuPacBio.contigs.fasta -p PacBio_mummer
delta-filter -1 PacBio_mummer.delta > PacBio_mummer.1delta
mummerplot --size large -layout --color -f --png PacBio_mummer.1delta -p PacBio_mummer

# mummer for illumina_read1
nucmer ${OUTDIR}/GCF_000005845.2_ASM584v2.fna ${OUTDIR}/scaffolds.fasta -p illumina_mummer
delta-filter -1 illumina_mummer.delta > illumina_mummer.1delta
mummerplot --size large -layout --color -f --png illumina_mummer.1delta -p illumina_mummer
