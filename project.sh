#!/bin/bash
#SBATCH --job-name=project                           # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=8:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/vjp98982/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=vjp98982@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/work/gene8940/vjp98982/project"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

if [ ! -d ${OUTDIR}/fastqc_out ]
then
    mkdir -p ${OUTDIR}/fastqc_out
fi

if [ ! -d ${OUTDIR}/untrimmed_asm ]
then
    mkdir -p ${OUTDIR}/untrimmed_asm
fi


module load SRA-Toolkit/2.9.6-1-centos_linux64
module load FastQC/0.11.9-Java-11
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
module load SPAdes/3.14.1-GCC-8.3.0-Python-3.7.4



# dowload and extract paired end illumina reads from ncbi for Lactarius indigo SRA file for run SRR5804120
#prefetch -O ${OUTDIR} SRR5804120
#fastq-dump --split-files --gzip ${OUTDIR}/SRR5804120.sra -O ${OUTDIR}

# dowload and extract paired end illumina reads from ncbi for Lactarius indigo SRA file for run SRR5804121
#prefetch -O ${OUTDIR} SRR5804121
#fastq-dump --split-files --gzip ${OUTDIR}/SRR5804121.sra -O ${OUTDIR}

# perform fastqc to assess read quality
#fastqc ${OUTDIR}/SRR5804120_1.fastq.gz ${OUTDIR}/SRR5804120_2.fastq.gz ${OUTDIR}/SRR5804121_1.fastq.gz ${OUTDIR}/SRR5804121_2.fastq.gz --outdir ${OUTDIR}/fastqc_out

# use trimgalore to improve the quality of the reads prior to spades assembly. Trim out adapter sequences

##### this is the assembly command from HW3, use as a model for the project
# assemble the E coli MG1655 genome using Illumina short read data with SPADES
#spades.py -t 8 -k 21,33,55,77 --isolate --memory 32 --pe1-1 ${OUTDIR}/illumina_read1.fastq.gz --pe1-2 ${OUTDIR}/illumina_read2.fastq.gz -o ${OUTDIR}

spades.py -t 8 -k 21,33,55,77 --isolate --memory 32 --pe1-1 ${OUTDIR}/SRR5804120_1.fastq.gz --pe1-2 ${OUTDIR}/SRR5804120_2.fastq.gz --pe2-1 ${OUTDIR}/SRR5804121_1.fastq.gz --pe2-2 ${OUTDIR}/SRR5804121_2.fastq.gz-o ${OUTDIR}/untrimmed_asm

# quast to assess assembly contiguity
